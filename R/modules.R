#BUGS MODULES#

#E.G., z[1:z_length] ~ nim_lm( A[1:A_length] + (xm[1:n] | A[1:A_length]), factors = 'A', priors = 'Normal')

library(nimble)
library(assertthat) ## to package dependencies


## nimble.lm, nim_lm, and lmPred Helper Functions ##

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


#Extracts formula From Environment
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


#Adds Parameter Symbol to Covariate
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


#Converts Term from randomEffects List to BUGS term and makes prior; Could Run This Function over All Random Terms in List
#term[[2]] is left of | term[[3]] is right of |
formulaRandom2BUGS <- function(term, fr, index.name, intcpt, fixedterms){

  out <- list()  #First Element of List Random Effect Terms, 2nd Priors

  #Random Intercept
  if(term[[2]]==1 & is.factor(fr[[as.character(term[[3]])]])){

    out[[1]] = substitute(EFFECTNAME[TERM[INDEX]],  #How to Keep double Brackets without Term Name in Front?
                          list(EFFECTNAME = make.effect.name(term[[3]]), #was (effect.name)
                               TERM = as.name(term[[3]]),
                               INDEX = index.name))

    groups <-  nlevels(as.factor((fr[[as.character(term[[3]])]])))

    if(intcpt==1){
      random <- make.random.prior(LHS2BUGSterm(make.effect.name(term[[3]]), index.name), mean = make.mean.name(term[[3]]),
                                  sig = make.sigma.name(term[[3]]))
      prior.mu <- make.prior(make.mean.name(term[[3]]), 0)
    }
    else{
      random <- make.random.prior(LHS2BUGSterm(make.effect.name(term[[3]]), index.name), mean = make.mean.name(term[[3]]),
                                  sig = make.sigma.name(term[[3]]))
      prior.mu <- make.prior(make.mean.name(term[[3]]),"Normal")
    }
    random.loop <- embedLinesInForLoop(random, index.name, start = 1, finish = as.numeric(groups))

    #Priors
    prior.sd <-  make.prior(make.sigma.name(term[[3]]),"Uniform")
    out[[2]] <- embedLinesInCurlyBrackets(list(random.loop, prior.sd, prior.mu))
  }

  #Random Slope
  if (as.character(term[[2]]) %in% names(fr)){
    out[[1]] = substitute(EFFECTNAME[TERM[INDEX]]*TERM2,  #How to Keep double Brackets without Term Name in Front?
                          list(EFFECTNAME = make.effect.name(term[[2]]), #was (effect.name)
                               TERM = as.name(term[[3]]),
                               TERM2 = LHS2BUGSterm(term[[2]], index.name),
                               INDEX = index.name))

    groups <-  nlevels(as.factor((fr[[as.character(term[[3]])]])))

    if(as.character(term[[2]]) %in% fixedterms){
      random <- make.random.prior(LHS2BUGSterm(make.effect.name(term[[2]]), index.name), mean = make.mean.name(term[[2]]),
                                  sig = make.sigma.name(term[[2]]))
      prior.mu <- make.prior(make.mean.name(term[[2]]), 0)
    }
    else{
      random <- make.random.prior(LHS2BUGSterm(make.effect.name(term[[2]]), index.name), mean = make.mean.name(term[[2]]),
                                  sig = make.sigma.name(term[[2]]))
      prior.mu <- make.prior(make.mean.name(term[[2]]),"Normal")
    }

    random.loop <- embedLinesInForLoop(random, index.name, start = 1, finish = as.numeric(groups))

    #Priors
    prior.sd <-  make.prior(make.sigma.name(term[[2]]),"Uniform")
    out[[2]] <- embedLinesInCurlyBrackets(list(random.loop, prior.sd, prior.mu))

  }

  return(out)
}


##Converts a Term from lmPred and Converts it to a BUGS Term
lmPred2BUGSterm <- function(termbracket, factors, index.name) {
  ## Need to Account for Interaction Terms.  Represented by : in attr(terms)
  #Limiting to Two Way Interactions for Now.  Indexing/Checks Gets Nasty for >2 Way Interactions.  Maybe should create ANOVA Module?

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

lmPredTerm2PRIOR <- function(termbracket, RHSname, factors, index.name, dropbase, prior){

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

    if(dropbase==TRUE){
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


##### Main nimble.lm Function ####

## Future Options: shape/scale input for priors, additional family link functions, random effects
##dropbase will drop lowest group for cateogrical variables.
##Need to Add Nested Random Effects Too

nimble.lm <- function(mod, dat=NULL, family=c("Normal"), priors=c("Normal","t","Uniform"), dropbase=TRUE){

  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)

  m <- match(c("mod","dat"), names(mf), 0L)

  ## construct a call
  mf <- mf[c(1L, m)]
  names(mf)[names(mf)=='dat'] <- 'data'
  names(mf)[names(mf)=='mod'] <- 'formula'
  mf$drop.unused.levels <- TRUE


  #Change name of function
  mf[[1L]] <- quote(stats::model.frame)

  formula <- as.formula(mf$formula) #RHSForm(mf$formula)

  #Intercept in Model
  int.add <- attr(terms(formula), 'intercept')

  fr.form <- subbars(formula)  #turns "|" to "+" so we can use stats::model.frame


  mf$formula <- fr.form
  fr <- eval(mf, parent.frame()) # Gets us mf with data filled in from dat argument or from calling environment as needed


  #Convert Character Variables to Factors
  for (i in 1:dim(fr)[2]){
    fr[,i] <- makeFac(fr[,i],char.only = TRUE)
  }


  #Check link and prior arguments are supported
  family = match.arg(family)
  prior = match.arg(priors)

  #Make nim_lm
  index.name <- quote(i)

  nobs = as.numeric(dim(fr)[1]) #number of rows/obs in data frame

  #Pull out Fixed Terms
  fixedTm <- nobars(formula)
  terms.to.add <- attr(terms(fixedTm), 'term.labels')

  #Pull Out Random Terms
  randomTerms <- fb(formula)

  #Make nim_lm formula
  LHS <- bracket(as.character(mf$formula[[2]]),fr)
  #Add RHS Terms
  if (int.add==1){
    RHS <- substitute(TERM[1:L], list(TERM=quote(intercept), L=nobs))
    ## append additional terms
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
      RHS <- addTerm(RHS, bracket.random(as.character(randomTerms[[iTerm]]),fr))
    }
  }

  if(length(randomTerms)>0 & !exists("RHS")){
    RHS <- bracket.random(randomTerms[[1]],fr)
    for (iTerm in seq_along(randomTerms[-1]) + 1){
      RHS <- addTerm(RHS, bracket.random(as.character(randomTerms[[iTerm]]),fr))
    }
  }

  RHS.LM <- call("nim_lm",RHS)
  if(length(names(grep("factor", sapply(fr, class), value=TRUE))) > 0) {RHS.LM[[3]] <- names(grep("factor", sapply(fr, class), value=TRUE))}
      else{RHS.LM[[3]] <- "None"}
  RHS.LM[[4]] <- cl$priors
  names(RHS.LM)[3:4] <- c("factors", "priors")

  return(list(LHS = LHS, RHS = RHS.LM))
}


## Main nim_lm Function ##
## Removing the MakeBUGS Modele Piece for Now

nim_lm <- function(LHS, RHS){
  RHSargs <- match.call(function(mod, factors, priors){}, RHS)

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
  names(r.var)[2:6] <- c("RHSmodel", "factors", "priors", "family", "link")
  X <- substitute(LHS ~ RHS, list(LHS = LHS2BUGSterm(make.pred.name(LHS[[2]]), index.name) , RHS = r.var))

  #Return Code
  newCode <- embedLinesInCurlyBrackets(lines = list(sd.prior, Y, X))
  return(list(code = newCode, LHS = LHS2BUGSterm(make.pred.name(LHS[[2]]), index.name) ,RHS = r.var))
}

## Make lmPred function ##
## Removing the MakeBUGS Modele Piece for Now
lmPred <- function(LHS, RHS){

  RHSargs <- match.call(function(RHSmodel, factors, priors, family, link){}, RHS)
  link <- match.arg(RHSargs[['link']], c('Identity','log', 'Poisson', 'Bernoulli'))  #Can Add More Link Functions

  #Generate Code
  if(link == 'Identity') {
    LHS <- LHS
  }

  terms <- strsplit(safeDeparse(RHSargs$RHSmodel), split='+', fixed=TRUE)[[1]]
  terms.list <- lapply(terms, FUN = function(x) {gsub(" ", "", x , fixed = TRUE)})

  index.name <- quote(i)

  RHS <- list()
  Priors <- list()

  #Generate Terms and Priors for Fixed Terms
  for(iTerm in seq_along(terms.list)){
    if(grepl('|', terms.list[[iTerm]], fixed=TRUE)==FALSE){
        RHS[[iTerm]] <- lmPred2BUGSterm(terms.list[[iTerm]], RHSargs$factors, index.name)
        Priors[[iTerm]] <- lmPredTerm2PRIOR(terms.list[[iTerm]], RHS[iTerm], factors=RHSargs$factors, index.name, dropbase=TRUE, prior=RHSargs$priors)
    }
  }

  RHS2 <- RHS[[1]]
  for (i in seq_along(RHS)[-1]){
    RHS2 <- addTerm(RHS2, RHS[[i]])
  }

  #Make Model Form and Embed in Loop
  fullLine <- substitute(LHS <- RHS, list(LHS = LHS, RHS = RHS2))
  forLoop <- embedLinesInForLoop(fullLine, index.name, start = 1, finish = 10)
  Priors.All <- embedLinesInCurlyBrackets(Priors)
  newCode <- embedLinesInCurlyBrackets(list(forLoop, Priors.All))

  return(list(code = newCode))
}

##Test Function

sillyData <- list(z = rpois(10,2), B = c(rep("colin",10)), q=rnorm(10), xm = rnorm(10), ym=rgamma(10,3), A = c(rep('Red',7), rep('Orange',1), rep('Blue', 2)))

#debug(nimble.lm)

test <- nimble.lm(mod = z ~ 1 + xm*A, dat = sillyData, priors="Normal", dropbase=TRUE)
test2 <- nim_lm(test$LHS, test$RHS)
test3 <- lmPred(test2$LHS, test2$RHS)


#Need to Think About how to Combine All 3 Functions to Return One Block of Expanded Code to Send to nimbleModel()





