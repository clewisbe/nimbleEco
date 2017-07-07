#BUGS MODULES#
library(nimble)
library(assertthat) ## to package dependencies


## nimble.lm, nim_lm, and lmPred Helper Functions ##

#Let's Deparse Run Over Multiple Lines
safeDeparse <- function(expr){
  ret <- paste(deparse(expr), collapse="")
  #rm whitespace
  gsub("[[:space:]][[:space:]]+", " ", ret)
}

#Return Number of Groups for Factor "A[1:Groups]"
no.groups <- function(term){
  splt <- as.numeric(strsplit(term, "[^0-9]+")[[1]])
  gp <- max(splt, na.rm = T) - min(splt, na.rm = T) + 1
  gp
}

#Removes Bracket Indexing for single [ ]
removeIndexing <- function(term) {
  return(gsub('\\[.*', '', term))
}

#Remove () and Split Random Effect Term by |
removeParen <- function(term) {
   out = gsub("[()]", "", term)
   return(unlist(strsplit(out, split='|',fixed = TRUE)))
}

#Make Bracket from Term with Indexing Running to Either No. Obs or No. Groups if Factor Variable
mk.bracket <- function(term, modframe){
  if(is.factor(modframe[[term]])){
    l = as.numeric(nlevels((modframe[[term]])))
  }
  else{
    l = as.numeric(length(modframe[[term]]))
  }
  substitute(TERM[1:L], list(TERM = as.name(term), L = l))
}


bracket <- function(term, mf2){
  if (grepl(":",term)==TRUE){
    terms.split <- unlist(strsplit(term, ":"))
    t.bracket = lapply(terms.split, mk.bracket, mf2)
    out <- prodTerm(t.bracket[[1]], t.bracket[[2]])

    ## append additional terms
    for(i in seq_along(t.bracket)[-c(1:2)]){
      out <- prodTerm(out, t.bracket[[i]])
    }
    out
  }
  else{
    out = mk.bracket(term, mf2)
    out
  }
}

#Same Bracket Function for Random Effects |
bracket.random <- function(term, mf){
  if(term[[2]]==1){
    lside=1}
  else{
    lside <- mk.bracket(term[[2]], mf)
  }
  rside <- mk.bracket(term[[3]], mf)
  out <- substitute((L|R), list(L=lside, R=rside))
  out
}

## These Functions Make the Name SomeFactor.X
make.effect.name <- function(term) {
  as.name(paste0(term, '.effect'))
}

make.coef.name <- function(term) {
  as.name(paste0('b.', term))
}

make.mean.name <- function(term) {
  as.name(paste0(term, '.mu'))
}

make.precision.name <- function(term) {
  as.name(paste0(term, '.tau'))
}

make.sigma.name <- function(term) {
  as.name(paste0(term, '.sigma'))
}

make.pred.name <- function(term) {
  as.name(paste0(term, '.pred'))
}


#Make Factor Variable
makeFac <- function(x,char.only=FALSE) {
  if (!is.factor(x) && (!char.only || is.character(x))) factor(x) else x
}


#Extracts Formula From Environment
RHSForm <- function(form,as.form=FALSE) {
  rhsf <- form[[length(form)]]
  if (as.form) reformulate(deparse(rhsf)) else rhsf
}

#Functions for Random Effects; Note These Functions Come From or are Modfied from lmer package

#Replaces '|' with '+'
subbars <- function(term)
{
  if (is.name(term) || !is.language(term)) return(term)
  if (length(term) == 2) {
    term[[2]] <- subbars(term[[2]])
    return(term)
  }
  stopifnot(length(term) >= 3)
  if (is.call(term) && term[[1]] == as.name('|'))
    term[[1]] <- as.name('+')
  for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
  term
}

#Find Random Effect Terms from Formula object
fb <- function(term)
{
  if (is.name(term) || !is.language(term)) return(NULL)
  if (identical(term[[1]], as.name("("))) return(fb(term[[2]]))
  stopifnot(is.call(term))
  if (term[[1]] == as.name('|')) return(term)
  if (length(term) == 2) return(fb(term[[2]]))
  c(fb(term[[2]]), fb(term[[3]]))
}

#Find Fixed Effect Terms
nobars <- function(term) {
  nb <- nobars_(term)  ## call recursive version
  if (is(term,"formula") && length(term)==3 && is.symbol(nb)) {
    ## called with two-sided RE-only formula:
    ##    construct response~1 formula
    nb <- reformulate("1",response=deparse(nb))
  }
  ## called with one-sided RE-only formula, or RHS alone
  if (is.null(nb)) {
    nb <- if (is(term,"formula")) ~1 else 1
  }
  nb
}

nobars_ <- function(term)
{
  if (!anyBars(term)) return(term)
  if (isBar(term)) return(NULL)
  if (isAnyArgBar(term)) return(NULL)
  if (length(term) == 2) {
    nb <- nobars_(term[[2]])
    if(is.null(nb)) return(NULL)
    term[[2]] <- nb
    return(term)
  }
  nb2 <- nobars_(term[[2]])
  nb3 <- nobars_(term[[3]])
  if (is.null(nb2)) return(nb3)
  if (is.null(nb3)) return(nb2)
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}

isBar <- function(term) {
  if(is.call(term)) {
    if((term[[1]] == as.name("|")) || (term[[1]] == as.name("||"))) {
      return(TRUE)
    }
  }
  FALSE
}

isAnyArgBar <- function(term) {
  if ((term[[1]] != as.name("~")) && (term[[1]] != as.name("("))) {
    for(i in seq_along(term)) {
      if(isBar(term[[i]])) return(TRUE)
    }
  }
  FALSE
}

anyBars <- function(term) {
  any(c('|','||') %in% all.names(term))
}


#Make Normal Diffuse Prior N ~ (0, 100,000), Uniform (0,1000), t(0,1,d), or point mass at 0
#We Could Also Allow an Argument to Specify the Shape and Scale Parameter?
make.prior <- function(term, prior,...){

  if (prior=="Normal"){
    substitute(TERM ~ dnorm(0, sd=1000), list(TERM=term))
  }

  else if (prior=="Uniform"){
    substitute(TERM ~ dunif(0, 1000), list(TERM=term))
  }

  else if (prior=="t"){
    substitute(TERM ~ dt(0, 1, d), list(TERM=term))
  }

  else if (prior=="Gamma"){
    substitute(TERM ~ dgamma(0.001, 0.001), list(TERM=term))
  }

  else{
    substitute(TERM <- 0, list(TERM=term))
  }
}

#Make Normal Prior for Random Effect
make.random.prior <- function(term, mean, sig){
  substitute(TERM ~ dnorm(MEAN, sd=SD), list(TERM=term, MEAN=mean, SD=sig))
}


#Adds Parameter Name to Covariate
make.param.name <- function(paramval, xvalue){
  substitute(B*X, list(B = paramval, X = xvalue))
}


## Takes Y (as a name) and constructs Y[i]
LHS2BUGSterm <- function(LHScode, indexName) {
  substitute(X[I], list(X = LHScode, I = indexName))
}

##Adds Term to a BUGS RHS.
addTerm <- function(currentTerms, newTerm) {
  substitute(A + B, list(A = currentTerms, B = newTerm))
}

##Multiply Two Terms
prodTerm <- function(currentTerms, newTerm) {
  substitute(A * B, list(A = currentTerms, B = newTerm))
}


## This will take a list of code lines and embed them in {}
embedLinesInCurlyBrackets <- function(lines) {
  as.call(c(list(quote(`{`)), lines))
}


## This will take a list of BUGS lines and embed them in a for loop
embedLinesInForLoop <- function(lines, indexName, start = 1, finish) {
  linesInBrackets <- embedLinesInCurlyBrackets(lines)
  rangeCall <- substitute(A:B, list(A = start, B = finish))
  ans <- substitute(for(INDEXNAME in RANGE) STUFF,
                    list(INDEXNAME = indexName,
                         RANGE = rangeCall,
                         STUFF = linesInBrackets))
  ans
}

#Make Starting Values for a Variable Passed To nimbleModel
make.init <- function(term, initialvalue, size){
  l = rep(initialvalue, size)
  init <- substitute(c(TERM,L), list(TERM = as.character(term), L = l))
  init
}


#Generate Initital Values for nimbleModel object
lmPred2init <- function(termbracket, factors, initial){

  term <- removeIndexing(termbracket)  #Remove RHS bracket from the Term

  if(term == "intercept"){
    return(make.init(make.coef.name(term), initial, 1))
  }

  if (length(strsplit(termbracket, '*', fixed = TRUE)[[1]]) > 2)  warning('Higher Than 2 Way Interactions Not Supported')

  if (length(strsplit(termbracket, '*', fixed = TRUE)[[1]]) == 2){
    interaction = unlist(lapply(strsplit(termbracket, '*', fixed = TRUE)[[1]], FUN = removeIndexing))

    if(interaction[1] %in% factors  | interaction[2] %in% factors){

      no.coef <- min(no.groups(unlist(strsplit(termbracket, '*', fixed = TRUE))[1]), no.groups(unlist(strsplit(termbracket, '*', fixed = TRUE))[2]))
      return(make.init(make.coef.name(paste0(interaction[1],interaction[2])), initial, no.coef))
    }

    if(all(interaction %in% factors)){
      no.coef <- no.groups(unlist(strsplit(termbracket, '*', fixed = TRUE))[1]) * no.groups(unlist(strsplit(termbracket, '*', fixed = TRUE))[2])
      return(make.init(make.coef.name(paste0(interaction[1],interaction[2])), initial, no.coef))
    }

    else{
    return(make.init(make.coef.name(paste0(interaction[1],interaction[2])), initial, 1))
    }
  }

  if(term %in% factors){
    return(make.init(make.coef.name(term), initial, no.groups(termbracket)))
    }
  else
    make.init(make.coef.name(term), initial, 1)
}

##Takes Term from lmPred and Converts it to a BUGS Term
lmPred2BUGSterm <- function(termbracket, factors, index.name) {
  ## Need to Account for Interaction Terms.  Represented by : in attr(terms)
  #Limiting to Two Way Interactions for Now.  Indexing/Checks Gets Nasty for >2 Way Interactions.  Maybe should create seperate ANOVA Module?

  term <- removeIndexing(termbracket)  #Remove RHS bracket from the Term

  if(term=="intercept"){
    return(make.coef.name(term))
  }

  if (length(strsplit(termbracket, '*', fixed = TRUE)[[1]]) > 2)  warning('Higher Than 2 Way Interactions Not Supported')

  if (length(strsplit(termbracket, '*', fixed = TRUE)[[1]]) == 2){
    interaction = unlist(lapply(strsplit(termbracket, '*', fixed = TRUE)[[1]], FUN = removeIndexing))

    if(interaction[1] %in% factors & !(interaction[2] %in% factors)){
      out = substitute(EFFECTNAME[TERM[INDEX]]*TERM2[INDEX],
                       list(EFFECTNAME = make.coef.name(paste0(interaction[1],interaction[2])), #Index Beta's on Covariates
                            TERM = as.name(interaction[1]),
                            TERM2 = as.name(interaction[2]),
                            INDEX = index.name))
      return(out)
    }

    if(interaction[2] %in% factors & !(interaction[1] %in% factors)){
      out = substitute(EFFECTNAME[TERM[INDEX]]*TERM2[INDEX],
                       list(EFFECTNAME = make.coef.name(paste0(interaction[1], interaction[2])),
                            TERM = as.name(interaction[2]),
                            TERM2 = as.name(interaction[1]),
                            INDEX = index.name))
      return(out)
    }

     if(all(interaction %in% factors)){
       out = substitute(EFFECTNAME[TERM[INDEX], TERM2[INDEX]],
                        list(EFFECTNAME = make.coef.name(paste0(interaction[1], interaction[2])),
                             TERM = as.name(interaction[2]),
                             TERM2 = as.name(interaction[1]),
                             INDEX = index.name))
       return(out)
     }

    else{
      out = substitute(TERM[INDEX]*TERM2[INDEX],
                       list(TERM = as.name(interaction[1]),
                            TERM2 = as.name(interaction[2]),
                            INDEX = index.name))}
    return(make.param.name(make.coef.name(paste0(interaction[1],interaction[2])),out))
  }

  if(term %in% factors){
    ## If term is a factor use BUGS indexing create b2[block[i]]
    out = substitute(EFFECTNAME[TERM[INDEX]],  #How to Keep double Brackets without Term Name in Front?
                     list(EFFECTNAME = make.coef.name(term) , #was (effect.name)
                          TERM = as.name(term),
                          INDEX = index.name))
    return(out)
  }
  else
    ## otherwise create BUGS term like x[i]
    out = substitute(TERM[INDEX],
                     list(TERM = as.name(term),
                          INDEX = index.name))

  make.param.name(make.coef.name(term),out)

}

## This Adds a Prior for Each Fixed Parameter on the RHS
## Still Need to Add Priors for Factor Interaction Variables; This will Require Double Indexing
lmPredTerm2PRIOR <- function(termbracket, RHSname, factors, index.name, prior, dropbase, int.check){

  term <- removeIndexing(termbracket)

  if(length(strsplit(as.character(RHSname), '*', fixed = TRUE)[[1]]) > 0){
  param.name <- as.name(gsub(" ", "", removeIndexing(strsplit(as.character(RHSname), '*', fixed=TRUE)[[1]][1]), fixed = TRUE))}

  else{
    param.name <- removeIndexing(RHSname)
  }

  if(term %in% factors) {
    ## if term is a factor, create priors for each level (note, dropbase=T will assign 0 prior to first level)

    out = list()  #Add Each Prior to List Object

    groups = no.groups(termbracket)  #Drop one group for model identification; Need to add Check if Intercept in Model Too

    if(dropbase==TRUE | int.check==1){
      out[[1]] <- make.prior(LHS2BUGSterm(param.name, 1), prior=0)
      st = 2
      idx = 2
    }

    else {st=1; idx=1}

    p = make.prior(LHS2BUGSterm(param.name, index.name), prior)
    out[[idx]] <- embedLinesInForLoop(p, index.name, start = st, finish = groups)
    return(embedLinesInCurlyBrackets(out))
  }

  #Priors for Interaction Terms with Factor Variable

  if (length(strsplit(termbracket, '*', fixed = TRUE)[[1]]) == 2){
    interaction = unlist(lapply(strsplit(termbracket, '*', fixed = TRUE)[[1]], FUN = removeIndexing))

  if(interaction[1] %in% factors | interaction[2] %in% factors){

      out = list()

      #Assume the Factor Group will Have a Lower Upper Index Number; Think About Better Way To Grab This
      groups = min(no.groups(strsplit(termbracket, '*', fixed=TRUE)[[1]][1]), no.groups(strsplit(termbracket, '*', fixed=TRUE)[[1]][2]))

      if(dropbase==TRUE){
        out[[1]] <- make.prior(LHS2BUGSterm(param.name, 1), prior=0)
        st = 2
        idx = 2
      }
      else {st = 1; idx = 1}

      p = make.prior(LHS2BUGSterm(param.name, index.name), prior)
      out[[idx]] <- embedLinesInForLoop(p, index.name, start=st, finish = groups)

      return(embedLinesInCurlyBrackets(out))
    }


    else {
      out = make.prior(make.coef.name(paste0(inter$V1,inter$V2)), prior)
    }
  }

  else{

    out = make.prior(param.name, prior)
  }

  return(out)
}


#lmPred to BUGS Expansion for Random Effect Term
lmPred2BUGSRandom <- function(termbracket, factors, index.name){

  term <- removeParen(termbracket)

  #Random Intercept
  if(term[[1]]==1 & removeIndexing(term[[2]]) %in% factors){

    out = substitute(EFFECTNAME[TERM[INDEX]],  #How to Keep double Brackets without Term Name in Front?
                     list(EFFECTNAME = make.coef.name(removeIndexing(term[[2]])), #was (effect.name)
                          TERM = as.name(removeIndexing(term[[2]])),
                          INDEX = index.name))
  }

  #Random Slope
  if (removeIndexing(term[[1]]) != 1 & removeIndexing(term[[2]]) %in% factors){
    out = substitute(EFFECTNAME[TERM[INDEX]]*TERM2,  #How to Keep double Brackets without Term Name in Front?
                     list(EFFECTNAME = make.effect.name(removeIndexing(term[[2]])), #was (effect.name)
                          TERM = as.name(removeIndexing(term[[2]])),
                          TERM2 = LHS2BUGSterm(as.name(removeIndexing(term[[1]])), index.name),
                          INDEX = index.name))

  }
  return(out)
}

#Prior for RANDOM EFFECTS Term
lmPredRandom2PRIOR<- function(termbracket, RHSname, factors, index.name, prior, dropbase, int.check){

  term <- removeParen(termbracket)
  groups = no.groups(term[2])
  param.name <- as.name(removeIndexing(RHSname))

  #Priors for Random Intercept
  if(length(strsplit(as.character(RHSname), '*', fixed = TRUE)) < 3){
    if(dropbase==TRUE | int.check==1){
      random <- make.random.prior(LHS2BUGSterm(param.name, index.name), mean = make.mean.name(param.name),
                                sig = make.sigma.name(param.name))
      prior.mu <- make.prior(make.mean.name(param.name), 0)
    }

    else{
      random <- make.random.prior(LHS2BUGSterm(param.name, index.name), mean = make.mean.name(param.name),
                                  sig = make.sigma.name(param.name))
      prior.mu <- make.prior(make.mean.name(param.name),"Normal")
    }

    random.loop <- embedLinesInForLoop(random, index.name, start = 1, finish = groups)

#Priors
    prior.sd <-  make.prior(make.sigma.name(param.name),"Uniform")
    out <- embedLinesInCurlyBrackets(list(random.loop, prior.sd, prior.mu))
    return(out)
}

  else{
    random <- make.random.prior(LHS2BUGSterm(param.name, index.name), mean = make.mean.name(param.name),
                              sig = make.sigma.name(param.name))
    prior.mu <- make.prior(make.mean.name(param.name), 0)
    #prior.mu <- make.prior(make.mean.name(param.name),"Normal")

  random.loop <- embedLinesInForLoop(random, index.name, start = 1, finish = groups)

  #Priors
  prior.sd <-  make.prior(make.sigma.name(param.name),"Uniform")
  out <- embedLinesInCurlyBrackets(list(random.loop, prior.sd, prior.mu))
  return(out)
  }
}





##### Main nimble.lm Function ####

## Future Options: shape/scale input for priors, additional family link functions, random effects
##dropbase will drop lowest group for cateogrical variables.
##Need to Add Nested Random Effects Too
nimble.lm <- function(mod, dat = NULL, family = c("Normal"), priors = c("Normal","t","Uniform"), dropbase = TRUE, niter = 10000, burnin = 1000,
                      initmcmc = 1, chains = 1, returncode = FALSE, returnsamp = FALSE){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)

  m <- match(c("mod","dat"), names(mf), 0L)

  #Construct Call
  mf <- mf[c(1L, m)]
  names(mf)[names(mf)=='dat'] <- 'data'
  names(mf)[names(mf)=='mod'] <- 'formula'
  mf$drop.unused.levels <- TRUE


  #Change Name of Function
  mf[[1L]] <- quote(stats::model.frame)

  formula <- as.formula(mf$formula)

  #Check Intercept in Model
  int.add <- attr(terms(formula), 'intercept')

  fr.form <- subbars(formula)  #turns "|" to "+" so we can use stats::model.frame

  mf$formula <- fr.form
  fr <- eval(mf, parent.frame()) # Gets us mf with data filled in from dat argument or from calling environment as needed

  #Convert Character Variables to Factors
  for (i in 1:dim(fr)[2]){
    fr[,i] <- makeFac(fr[,i], char.only = TRUE)
  }


  #Check Link and Prior Args
  family = match.arg(family)
  prior = match.arg(priors)

  #Make nim_lm
  index.name <- quote(i) #One Index for Now

  nobs = as.numeric(dim(fr)[1]) #number of rows/obs in data frame

  #Pull Out Fixed Terms
  fixedTm <- nobars(formula)
  terms.to.add <- attr(terms(fixedTm), 'term.labels')

  #Pull Out Random Terms
  randomTerms <- fb(formula)

  #Make nim_lm Formula
  LHS <- bracket(as.character(mf$formula[[2]]),fr)

  #Add RHS Terms
  if (int.add==1){
    RHS <- substitute(TERM[1:L], list(TERM=quote(intercept), L=nobs))
    #Append Additional Terms
    for(iTerm in seq_along(terms.to.add)) {
      RHS <- addTerm(RHS, bracket(as.character(terms.to.add[iTerm]), fr))
    }
  }

  if (int.add!=1  & length(terms.to.add>0)){
    RHS <- bracket(terms.to.add[1], fr)
    for(iTerm in seq_along(terms.to.add[-1]) + 1){
      RHS <- addTerm(RHS, bracket(as.character(terms.to.add[iTerm]), fr))
    }
  }

  #Add Random Terms
  if (length(randomTerms)>0 & exists("RHS")){
    for (iTerm in seq_along(randomTerms)){
      RHS <- addTerm(RHS, bracket.random(as.character(randomTerms[[iTerm]]), fr))
    }
  }

  if(length(randomTerms)>0 & !exists("RHS")){
    RHS <- bracket.random(randomTerms[[1]], fr)
    for (iTerm in seq_along(randomTerms[-1]) + 1){
      RHS <- addTerm(RHS, bracket.random(as.character(randomTerms[[iTerm]]), fr))
    }
  }


  RHS.LM <- call("nim_lm",RHS)
  if(length(names(grep("factor", sapply(fr, class), value=TRUE))) > 0) {RHS.LM[[3]] <- names(grep("factor", sapply(fr, class), value=TRUE))}
      else{RHS.LM[[3]] <- "None"}
  RHS.LM[[4]] <- cl$priors
  RHS.LM[[5]] <- cl$dropbase
  names(RHS.LM)[3:5] <- c("factors", "priors","dropbase")


  #Expand Code Using BUGS Modules Option
  lm.expand <- nim_lm$process(LHS, RHS.LM)
  pred.expand <- lmPred$process(lm.expand$LHS, lm.expand$RHS)

  #Full BUGS Code Expansion
  full.expand <- embedLinesInCurlyBrackets(lines = list(lm.expand$probmod, pred.expand$code))

  #Make nimbleModel Object and Run MCMC
  #Set Initial Values
  terms.init <- strsplit(safeDeparse(lm.expand$RHS$RHSmodel), split='+', fixed=TRUE)[[1]]
  terms.init.list <- lapply(terms.init, FUN = function(x) {gsub(" ", "", x , fixed = TRUE)})

  nimInitial <- list() #Need to Add for Random Terms Too

  for(iTerm in seq_along(terms.init.list)){
    if(grepl('|', terms.init.list[[iTerm]], fixed=TRUE) == FALSE){
    nimInitial[[iTerm]] <- lmPred2init(terms.init.list[[iTerm]], factors = RHS.LM$factors, initial = initmcmc)
    }
  }


  values <- lapply(nimInitial, `[[`, 3)
  names(values) <- lapply(nimInitial, `[[`, 2)

  #Add Prior for Sigma
  values$sigma <- initmcmc


  #Define Data for nimModel function
  response.var = list(fr[,1])
  names(response.var) = names(fr)[1]

  #Make Factor Variables Numeric
  constants <- fr[,c(-1), drop=FALSE]
  indx <- sapply(constants, is.factor)
  constants[indx] <- lapply(constants[indx], function(x) as.numeric(as.factor(x)))

  browser()

  #Make Nimble Model Object, Compile, and run MCMC
  nimMod.obj <- nimbleModel(code = full.expand, inits = values, constants = constants, data = response.var)
  nimMod.C <- compileNimble(nimMod.obj)

  #Set Up MCMC
  config.Mod <- configureMCMC(nimMod.C, print = TRUE)
  mod.MCMC <- buildMCMC(config.Mod)
  C.mod.MCMC <- compileNimble(mod.MCMC)
  samplesList <- runMCMC(C.mod.MCMC, niter = niter, nburnin = burnin, nchains = chains)

  #Extract Samples
  samples <- as.matrix(C.mod.MCMC$mvSamples)

  out <- list("Summary" = summary(samples))
  if(returncode==TRUE){out[["BUGScode"]] = full.expand}
  if(returnsamp==TRUE){out[["Samples"]] = samplesList}

  return(out)
}


## BUGS Modules Functions ##
makeBUGSmodule <- function(fun) {

  ans <- structure(list(process = fun), class = "BUGSmodule")
  ans
}

nim_lm <- makeBUGSmodule(
  function(LHS, RHS) {
  RHSargs <- match.call(function(mod, factors, priors, dropbase){}, RHS)

  #Set Up LHS
  index.name <- quote(i)
  l.var <- LHS2BUGSterm(LHS[[2]], index.name)
  pred.d <- substitute(dnorm(MEAN, sd=SIGMA), list(MEAN=LHS2BUGSterm(make.pred.name(LHS[[2]]), index.name), SIGMA=quote(sigma)))
  meanfctn <-  substitute(LHS ~ RHS, list(LHS = l.var, RHS = pred.d))
  Y <- embedLinesInForLoop(meanfctn, index.name, start=LHS[[3]][[2]], finish=LHS[[3]][[3]])
  sd.prior <- quote(sigma ~ dunif(0, 1000))

  #Set up RHS
  r.var <- call("lmPred", RHSargs$mod)
  r.var[[3]] <- RHSargs$factors
  r.var[[4]] <- RHSargs$priors
  r.var[[5]] <- "Normal"
  r.var[[6]] <- "Identity"
  r.var[[7]] <- RHSargs$dropbase
  names(r.var)[2:7] <- c("RHSmodel", "factors", "priors", "family", "link", "dropbase")
  X <- substitute(LHS ~ RHS, list(LHS = LHS2BUGSterm(make.pred.name(LHS[[2]]), index.name) , RHS = r.var))

  #Return Code
  newCode <- embedLinesInCurlyBrackets(lines = list(sd.prior, Y, X))
  top.expand <- embedLinesInCurlyBrackets(lines = list(sd.prior, Y))
  return(list(code = newCode, LHS = LHS2BUGSterm(make.pred.name(LHS[[2]]), index.name), RHS = r.var, probmod = top.expand))
})

## Make lmPred function ##
lmPred <- makeBUGSmodule(
  function(LHS, RHS) {

  RHSargs <- match.call(function(RHSmodel, factors, priors, family, link, dropbase){}, RHS)
  link <- match.arg(RHSargs[['link']], c('Identity','log', 'Poisson', 'Bernoulli'))  #Can Add More Link Functions

  #Generate Code
  if(link == 'Identity') {
    LHS <- LHS
  }

  terms <- strsplit(safeDeparse(RHSargs$RHSmodel), split='+', fixed=TRUE)[[1]]
  terms.list <- lapply(terms, FUN = function(x) {gsub(" ", "", x , fixed = TRUE)})

  RHS.fix.tm <- list()
  RHS.random.tm <- list()

  for(iTerm in seq_along(terms.list)){
    if(grepl('|', terms.list[[iTerm]], fixed=TRUE) == FALSE){
      RHS.fix.tm[[iTerm]] <- terms.list[[iTerm]]}
    else{
      RHS.random.tm[[iTerm]] <- terms.list[[iTerm]]
    }
  }

  #Drop Null Slots.  Look into Vectorizing
  RHS.fix.tm <- RHS.fix.tm[!sapply(RHS.fix.tm,is.null)]
  RHS.random.tm <- RHS.random.tm[!sapply(RHS.random.tm,is.null)]


  index.name <- quote(i)
  int.check <- sum(!(is.na(sapply(terms.list, function(x) pmatch("intercept", x)))))


  #Generate Terms and Priors for Fixed Terms
  RHS.fix <- list()
  Priors.fix <- list()
  for(iTerm in seq_along(RHS.fix.tm)){
        RHS.fix[[iTerm]] <- lmPred2BUGSterm(RHS.fix.tm[[iTerm]], RHSargs$factors, index.name)
        Priors.fix[[iTerm]] <- lmPredTerm2PRIOR(RHS.fix.tm[[iTerm]], RHS.fix[iTerm], factors = RHSargs$factors, index.name, prior = RHSargs$priors,
                                            dropbase = RHSargs$dropbase, int.check = int.check)
    }

  #Generate Terms and Priors for Random Terms
  RHS.random <- list()
  Priors.random <- list()
  for(iTerm in seq_along(RHS.random.tm)){
      RHS.random[[iTerm]] <- lmPred2BUGSRandom(RHS.random.tm[[iTerm]], RHSargs$factors, index.name)
      Priors.random[[iTerm]] <- lmPredRandom2PRIOR(RHS.random.tm[[iTerm]], RHS.random[iTerm], factors = RHSargs$factors, index.name, prior = RHSargs$priors,
                                          dropbase = RHSargs$dropbase, int.check = int.check)
  }


  Priors.Combo <- append(Priors.fix, Priors.random)
  Terms.Combo <- append(RHS.fix, RHS.random)

  RHS2 <- Terms.Combo[[1]]
  for (i in seq_along(Terms.Combo)[-1]){
    RHS2 <- addTerm(RHS2, Terms.Combo[[i]])
  }

  #Make Model Form and Embed in Loop
  fullLine <- substitute(LHS <- RHS, list(LHS = LHS, RHS = RHS2))
  forLoop <- embedLinesInForLoop(fullLine, index.name, start = 1, finish = 10)
  Priors.All <- embedLinesInCurlyBrackets(Priors.Combo)
  newCode <- embedLinesInCurlyBrackets(list(forLoop, Priors.All))

  return(list(code = newCode))
})






