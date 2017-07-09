### NIMBLE ABUNDANCE N Mixture Models for Repeated count Data ###
library(nimble)
#library(devtools)
#install_github("nimble-dev/nimble", ref = "devel", subdir = "packages/nimble")
library(tidyverse)
library(coda)
nimbleOptions(enableBUGSmodules = TRUE)

makeBUGSmodule <- function( fun ) {
  ans <- structure(list(process = fun), class = "BUGSmodule")
  ans
}

#Helper function to gather data
gatherfctn <- function(mydata, key.col, val.col, gather.cols) {
  new.data <- gather_(data = mydata,
                      key_col = key.col,
                      value_col = val.col,
                      gather_cols = colnames(mydata)[gather.cols])
  return(new.data)
}


## Make Tidy Data From 3 Abundance Data Sets ##
## y is matrix, sitevars is matrix, obsvars is named list
maketidy <- function(y, sitevars = NULL, obsvars = NULL){
  final.df <- list()

  if(is.null(y)) {stop('Count Data is Required')}
  y <- as.data.frame(y)
  L = dim(y)[2]
  VisitID <- paste("Visit", 1:L, sep="")
  colnames(y) <- VisitID
  y$Site <- as.numeric(rownames(y))
  long.y <- gather(y, Visit, Count, 1:L)


  #Site Variables
  if (!is.null(sitevars)){
    sitevars <- as.data.frame(sitevars)
    sitevars$Site <- as.numeric(rownames(sitevars))

  #Check Dimensions
  if(dim(y)[1] != dim(sitevars)[1]) {stop('Number of Rows of Y and Site Variables Do Not Match')}

  #Join Response with Site Level Variables
  final.df[[1]] <- left_join(long.y, sitevars, by = "Site")
  }

  else{
    final.df[[1]] <- long.y
  }

  #Combine Observation Level Data
  if (!is.null(obsvars)){
    varnames.obs <- names(obsvars)

    for (i in seq_along(obsvars)){
      mf = as.data.frame(obsvars[[i]])
      if(L != dim(mf)[2]){stop('Number of Visits of Y and Obs Covariate Variables Do Not Match')}
      colnames(mf) <- VisitID
      mf$Site <- as.numeric(rownames(mf))
      final.df[[i + 1]] <- gatherfctn(mf, "Visit", varnames.obs[i], 1:L)
      }
  }
  return(Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Site", "Visit")), final.df))
}



## BUGS Code Generation Helper Functions ##

#Takes Y (as a name) and constructs Y[i]
LHS2BUGSterm <- function(LHScode, indexName) {
  substitute(X[I], list(X = LHScode, I = indexName))
}

#Make the Name X.effect/mu/tau
make.effect.name <- function(term, level) {
  if (level == site){
  as.name(paste0(term, '.effect.s'))
  }
  else{
  as.name(paste0(term, '.effect.o'))
  }
}

#To avoid repeated names, seperate naming for Site and Obs level parameters
make.coef.name <- function(term, level) {
  if (as.character(level) == "site"){
    as.name(paste0('s.', term))
  }
  else{
  as.name(paste0('o.', term))
  }
}

make.coef.obs <- function(term) {
  as.name(paste0('o.', term))
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

make.param.name <- function(paramval, xvalue){
  substitute(B*X, list(B = paramval, X = xvalue))
}

#Make Normal Prior for Random Effect
make.random.prior <- function(term, mean, sig){
  substitute(TERM ~ dnorm(MEAN, sd=SD), list(TERM=term, MEAN=mean, SD=sig))
}

#Function to Make Many Priors
make.prior <- function(term, prior,...){

  if (prior=="Normal"){
    substitute(TERM ~ dnorm(0, sd=1000), list(TERM=term))
  }

  else if (prior=="Uniform"){
    substitute(TERM ~ dunif(0, 1), list(TERM=term))
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


#Lets Deparse Run Over Multiple Lines
safeDeparse <- function(expr){
  ret <- paste(deparse(expr), collapse="")
  #rm whitespace
  gsub("[[:space:]][[:space:]]+", " ", ret)
}

#Create A[1:L] where L is number of levels for factor
factor.bracket <- function(term, df){
  if(is.factor(df[[term]])){
    l = as.numeric(nlevels((df[[term]])))
    return(substitute(TERM[1:L], list(TERM = as.name(term), L = l)))
  }
}

#Pulls out formula from model frame
RHSForm <- function(form, as.form=FALSE){
  rhsf <- form[[length(form)]]
  if (as.form) reformulate(deparse(rhsf)) else rhsf
}

#Remove Parentheses from term
removeParen <- function(term) {
  out = gsub("[()]", "", term)
  return(unlist(strsplit(out, split='|',fixed = TRUE)))
}

#Adds Term to a BUGS RHS.
addTerm <- function(currentTerms, newTerm) {
  substitute(A + B, list(A = currentTerms, B = newTerm))
}

prodTerm <- function(currentTerms, newTerm) {
  substitute(A * B, list(A = currentTerms, B = newTerm))
}

#Make Factor Variable
makeFac <- function(x,char.only=FALSE) {
  if (!is.factor(x) && (!char.only || is.character(x))) factor(x) else x
}

#Put Block Indexing Around Term
bracket <- function(term, index){
  substitute(TERM[INDEX], list(TERM = term, INDEX = index))
}

#Make Double Bracket for Factor A[i,j]
doublebracket <- function(term, index1, index2){
  substitute(TERM[INDEX1, INDEX2], list(TERM = term, INDEX1 = index1, INDEX2 = index2))
}


#Removes Bracket Indexing for single [ ]
removeIndexing <- function(term) {
  return(gsub('\\[.*', '', term))
}

make.interaction <- function(interaction, factors, index.name){
  if(interaction %in% factors){
    out = substitute(TERM[INDEX],
                     list(TERM = as.name(interaction),
                          INDEX = index.name))
    out <- substitute((TERM), list(TERM = out))

  }
  else{
    out = substitute(TERM[INDEX],
                     list(TERM = as.name(interaction),
                          INDEX = index.name))}
  return(out)
}


no.groups <- function(term){
  term <- as.character(term)
  splt <- as.numeric(strsplit(term, "[^0-9]+")[[1]])
  return(max(splt, na.rm = T) - min(splt, na.rm = T) + 1)
}


# This will take a list of code lines and embed them in {}
embedLinesInCurlyBrackets <- function(lines) {
  as.call(c(list(quote(`{`)), lines))
}


#Put Lines in a For Loop
embedLinesInForLoop <- function(lines, indexName, start = 1, finish) {
  linesInBrackets <- embedLinesInCurlyBrackets(lines)
  rangeCall <- substitute(A:B, list(A = start, B = finish))
  ans <- substitute(for(INDEXNAME in RANGE) STUFF,
                    list(INDEXNAME = indexName,
                         RANGE = rangeCall,
                         STUFF = linesInBrackets))
  ans
}


#Find Fixed Effect Terms (from LME4 Package)
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

#Make nim_glm object to pass to BUGS Module for Processing
make.glm <- function(RHS, factors = NULL, cl, mixture, link, N = NULL, Site = NULL, level = NULL){
  RHS.out <- call("nim_glm", RHSForm(RHS))
  RHS.out[[3]] <- factors
  RHS.out[[4]] <- as.name(cl$priors)
  RHS.out[[5]] <- cl$dropbase
  RHS.out[[6]] <- mixture
  RHS.out[[7]] <- link
  RHS.out[[8]] <- level
  names(RHS.out)[3:8] <- c("factors", "priors","dropbase", "family", "link", "level")
  if (mixture %in% c("Binomial", "Bernoulli")){
    RHS.out[[9]] <- list(N, Site)
    names(RHS.out)[9] <- c("args")
  }
  return(RHS.out)
}

## nim_glm BUGS Module ##
#Returns Full Code, LHS, RHS, and the Probability Model
nim_glm <- makeBUGSmodule(
  function(LHS, RHS) {
    RHSargs <- match.call(function(mod, factors, priors, family, link, dropbase, level, args){}, RHS)

    #Set Up LHS
    index.name <- quote(i)

    if (RHSargs$family %in% c("Poisson", "NB", "ZIP", "Bernoulli")) {distn <- quote(dpois); link = quote(log); param = quote(lambda);
        pred.d <- substitute(DIST(PARAM), list(DIST = distn, PARAM=LHS2BUGSterm(param, index.name)))
    }

    if (RHSargs$family =="Binomial") {distn <- quote(dbinom); link = quote(logit); param = quote(p);
        pred.d <- substitute(DIST(PARAM, N), list(DIST = distn, PARAM=LHS2BUGSterm(param, index.name),
          N = bracket(RHSargs$args[[1]],bracket(RHSargs$args[[2]],quote(i)))))
    }

    meanfctn <-  substitute(LHS ~ RHS, list(LHS = LHS2BUGSterm(LHS[[2]], index.name), RHS = pred.d))
    meanfctn.loop <- embedLinesInForLoop(meanfctn, index.name, start = 1, finish=LHS[[3]][[3]])

    if (RHSargs$family =="NB"){
      dispersion <- substitute(LAMBDA <- MU*DISP, list(LAMBDA = LHS2BUGSterm(param, index.name), DISP = LHS2BUGSterm(quote(rho), index.name),
                                                       MU = LHS2BUGSterm(quote(mu), index.name)))
      rho.dist <- substitute(X ~ dgamma(alpha, alpha), list(X = LHS2BUGSterm(quote(rho), index.name)))
      meanfctn.no.prior <- embedLinesInForLoop(list(meanfctn, dispersion, rho.dist), index.name, start = 1, finish=LHS[[3]][[3]])

      rho.param <- substitute(ALPHA <- exp(LOGALPHA), list(ALPHA = quote(alpha), LOGALPHA = quote(logalpha)))
      rho.prior <- make.prior(quote(logalpha), "Normal")

      meanfctn.loop <- embedLinesInCurlyBrackets(list(rho.prior, rho.param, meanfctn.no.prior))
    }

    if (RHSargs$family =="ZIP"){
      zip <- substitute(S ~ dbern(THETA), list(S = LHS2BUGSterm(quote(s), index.name), THETA = quote(theta)))
      zip.prior <- make.prior(quote(theta), "Uniform")
      distn <- quote(dpois); link = quote(log); param = quote(lambda);
      pred.d <- substitute(DIST(PARAM*ZIP), list(DIST = distn, PARAM=LHS2BUGSterm(param, index.name), ZIP = LHS2BUGSterm(quote(s), index.name)))
      meanfctn <-  substitute(LHS ~ RHS, list(LHS = LHS2BUGSterm(LHS[[2]], index.name), RHS = pred.d))
      meanfctn.no.prior <- embedLinesInForLoop(list(meanfctn, zip), index.name, start = 1, finish=LHS[[3]][[3]])
      meanfctn.loop <- embedLinesInCurlyBrackets(list(zip.prior, meanfctn.no.prior))
    }

    if (RHSargs$family =="Bernoulli"){
      if(RHSargs$level == quote(site)){
        param <- quote(psi)
        presence <- substitute(dbern(PSI), list(PSI = LHS2BUGSterm(param, index.name)))
      }
      if(RHSargs$level == quote(obs)){
        param <- quote(p)
        presence <- substitute(dbern(Z*P), list(Z = bracket(quote(z), bracket(RHSargs$args[[2]],quote(i))), P = LHS2BUGSterm(param, index.name)))
      }
      meanfctn <-  substitute(LHS ~ RHS, list(LHS = LHS2BUGSterm(LHS[[2]], index.name), RHS = presence))
      meanfctn.loop <- embedLinesInForLoop(meanfctn, index.name, start = 1, finish=LHS[[3]][[3]])
    }

    #Set up RHS
    lpred <- call("lmPred", RHSargs$mod)
    lpred[[3]] <- RHSargs$factors
    lpred[[4]] <- RHSargs$priors
    lpred[[5]] <- as.name(RHSargs$link)
    lpred[[6]] <- RHSargs$dropbase
    lpred[[7]] <- RHSargs$level
    names(lpred)[3:7] <- c("factors", "priors", "link", "dropbase", "level")

    if (RHSargs$family =="NB"){
    lpred.full <- substitute(LHS <- RHS, list(LHS = bracket(quote(mu), LHS[[3]]), RHS = lpred))
    LHS.Param <- bracket(quote(mu), LHS[[3]])
    }

    else{
    lpred.full <- substitute(LHS <- RHS, list(LHS = bracket(param, LHS[[3]]), RHS = lpred))
    LHS.Param <- bracket(param, LHS[[3]])
    }

    #Return Code
    newCode <- embedLinesInCurlyBrackets(lines = list(meanfctn.loop, lpred.full))
    return(list(code = newCode, LHS = LHS.Param, RHS = lpred, prob.mod = meanfctn.loop))
  })


#Takes Term from lmPred and Converts it to a BUGS Term
#Note: : is not supported for interaction term, just *
lmPred2BUGSterm <- function(term, factors, index.name, level) {

  factors.name <- lapply(factors, removeIndexing)
  factors.name <- unlist(lapply(factors.name, '[[',2))
  if(term=="1"){
    return(make.coef.name(quote(intercept), level))
  }

  if (length(strsplit(term, '*', fixed = TRUE)[[1]]) == 2){
    interaction = unlist(strsplit(term, '*', fixed=TRUE))

    if(interaction[1] %in% factors.name & !(interaction[2] %in% factors.name)){
      out = substitute(EFFECTNAME[TERM[INDEX]]*TERM2[INDEX],
                       list(EFFECTNAME = make.effect.name(paste0(interaction[1],interaction[2]), level), #Index Beta's on Covariates
                            TERM = as.name(interaction[1]),
                            TERM2 = as.name(interaction[2]),
                            INDEX = index.name))
      return(out)
    }

    if(interaction[2] %in% factors.name & !(interaction[1] %in% factors.name)){
      out = substitute(EFFECTNAME[TERM[INDEX]]*TERM2[INDEX],
                       list(EFFECTNAME = make.effect.name(paste0(interaction[1], interaction[2]), level),
                            TERM = as.name(interaction[2]),
                            TERM2 = as.name(interaction[1]),
                            INDEX = index.name))
      return(out)
    }

    if(all(interaction %in% factors.name)){
      out = substitute(EFFECTNAME[TERM[INDEX], TERM2[INDEX]],
                       list(EFFECTNAME = make.effect.name(paste0(interaction[1], interaction[2]), level),
                            TERM = as.name(interaction[2]),
                            TERM2 = as.name(interaction[1]),
                            INDEX = index.name))
      return(out)
    }
    else{
      out = substitute(EFFECTNAME*TERM[INDEX]*TERM2[INDEX],
                       list(EFFECTNAME = make.coef.name(paste0(interaction[1], interaction[2]), level),
                            TERM = as.name(interaction[1]),
                            TERM2 = as.name(interaction[2]),
                            INDEX = index.name))
      return(out)
    }
  }

    if (length(strsplit(term, '*', fixed = TRUE)[[1]])  > 2){
      interaction = unlist(strsplit(term, '*', fixed=TRUE))

      tm.name <- make.coef.name(paste(interaction, collapse = ""), level)

      if(interaction[1] %in% factors){
        out = substitute(EFFECTNAME[TERM[INDEX]],
                       list(EFFECTNAME = tm.name, #Index Beta's on Covariates
                            TERM = as.name(interaction[1]),
                            INDEX = index.name))}
      else{
        out = substitute(EFFECTNAME*TERM[INDEX],
                         list(EFFECTNAME = tm.name, #Index Beta's on Covariates
                              TERM = as.name(interaction[1]),
                              INDEX = index.name))
      }

    for(i in seq_along(interaction[-1])) {
      out <- prodTerm(out, make.interaction(interaction[[i + 1]], factors, index.name))
    }
      return(out)
  }

  if(term %in% factors.name){
    ## If term is a factor use BUGS indexing create b2[block[i]]
    out = substitute(EFFECTNAME[TERM[INDEX]],  #How to Keep double Brackets without Term Name in Front?
                     list(EFFECTNAME = make.effect.name(term, level),
                          TERM = as.name(term),
                          INDEX = index.name))
    return(out)
  }
  else{
    ## otherwise create BUGS term like x[i]
    out = substitute(TERM[INDEX],
                     list(TERM = as.name(term),
                          INDEX = index.name))

  make.param.name(make.coef.name(term, level),out)}

}

#lmPred to BUGS Expansion for Random Effect Term
lmPred2BUGSRandom <- function(termparen, factors, index.name, level){

  term <- removeParen(termparen)

  factors.name <- lapply(factors, removeIndexing)
  factors.name <- unlist(lapply(factors.name, '[[',2))

  #Random Intercept
  if(term[[1]]==1 & removeIndexing(term[[2]]) %in% factors.name){

    out = substitute(EFFECTNAME[TERM[INDEX]],  #How to Keep double Brackets without Term Name in Front?
                     list(EFFECTNAME = make.coef.name(removeIndexing(term[[2]]), level), #was (effect.name)
                          TERM = as.name(removeIndexing(term[[2]])),
                          INDEX = index.name))
  }

  #Random Slope
  if (removeIndexing(term[[1]]) != 1 & removeIndexing(term[[2]]) %in% factors.name){
    out = substitute(EFFECTNAME[TERM[INDEX]]*TERM2,
                     list(EFFECTNAME = make.effect.name(paste0(term[1], term[2]), level),
                          TERM = as.name(removeIndexing(term[[2]])),
                          TERM2 = LHS2BUGSterm(as.name(removeIndexing(term[[1]])), index.name),
                          INDEX = index.name))

  }

  if(!(term[[2]] %in% factors.name)){stop('Conditioned Term is Not a Factor')}
  return(out)
}

#Make Prior for fixed Term
lmPredTerm2PRIOR <- function(term, RHSname, factors, index.name, prior, dropbase, int.check){
  priors = list()

  factors.name <- lapply(factors, removeIndexing)
  factors.name <- unlist(lapply(factors.name, '[[',2))


  if(length(strsplit(as.character(RHSname), '*', fixed = TRUE)[[1]]) > 0){
    param.name <- as.name(removeIndexing((gsub(" ", "", strsplit(as.character(RHSname), '*', fixed = TRUE)[[1]][[1]]))))
    }

  else{
    param.name <- as.name(removeIndexing(RHSname))
  }

  if(term %in% factors.name) {
    ## if term is a factor, create priors for each level (note, dropbase=T will assign 0 prior to first level)

    fct.index = factors[which(as.character(term) %in% factors.name)]
    groups = no.groups(fct.index)  #Drop one group for model identification; Need to add Check if Intercept in Model Too

    if(dropbase==TRUE | int.check==1){
      priors[[1]] <- make.prior(LHS2BUGSterm(param.name, 1), prior = 0)
      st = 2
      idx = 2
    }

    else {st=1; idx=1}

    p = make.prior(LHS2BUGSterm(param.name, index.name), prior)
    priors[[idx]] <- embedLinesInForLoop(p, index.name, start = st, finish = groups)
    initial <- list(groups)
    names(initial) <- quote(param.name)
    #initial[[1]] <- cbind(as.character(param.name), groups)
    return(c(embedLinesInCurlyBrackets(priors), initial))
  }

  #Priors for Interaction Terms with Factor Variable
  if(length(strsplit(as.character(term), '*', fixed = TRUE)[[1]]) == 2){
      interaction = unlist(strsplit(as.character(term), '*', fixed=TRUE))
#Taking this Piece Our for Now.  Determining Sum to 0 constraints for 2 and 3 Level Anova gets messy quick.
#     if(interaction[1] %in% factors.name & interaction[2] %in% factors.name){  #Do We want to Allow Users to Have Interaction Terms for Multiple Factors?
#         out <- list()
#
#         group1 <- no.groups(factors[which(interaction[1] %in% factors.name)])
#         group2 <- no.groups(factors[which(interaction[2] %in% factors.name)])
#
#         if(dropbase==TRUE | int.check==1){  #Identifiability Constraints
#           out[[1]] <- make.prior(doublebracket(param.name, 1, 1), prior = 0)
#           out[[2]] <- embedLinesInForLoop(make.prior(doublebracket(param.name, 1, index.name), prior = 0), index.name, start = 2, finish = group1)
#           out[[3]] <- embedLinesInForLoop(make.prior(doublebracket(param.name, index.name, 1), prior = 0), index.name, start = 2, finish = group2)
#           out[[4]] <- embedLinesInForLoop(embedLinesInForLoop(make.prior(doublebracket(param.name, index.name, quote(j)), prior = "Normal"), index.name, start = 2, finish = group2),
#                                            quote(j), start = 2, finish = group1)
#         }
#         return(out)
#         }


    if(interaction[1] %in% factors.name | interaction[2] %in% factors.name){

      #Assume the Factor Group will Have a Lower Upper Index Number; Think About Better Way To Grab This
      fct.index = factors[which(interaction %in% factors.name)]
      groups = no.groups(fct.index)

      if(dropbase==TRUE){
        priors[[1]] <- make.prior(LHS2BUGSterm(param.name, 1), prior=0)
        st = 2
        idx = 2
      }
      else {st = 1; idx = 1}

      p = make.prior(LHS2BUGSterm(param.name, index.name), prior)
      priors[[idx]] <- embedLinesInForLoop(p, index.name, start=st, finish = groups)
      initial <- list(groups)
      names(initial) <- quote(param.name)
      return(c(embedLinesInCurlyBrackets(priors), initial))
    }

      else{
        priors = make.prior(param.name, prior)
        initial = list(1)
        names(initial) <- quote(param.name)
        #initial[[1]] = cbind(as.character(param.name), 1)
        return(c(priors, initial))
      }
  }

    else{
      priors = make.prior(param.name, prior)
      initial = list(1)
      names(initial) <- quote(param.name)
      #initial[[1]] = cbind(as.character(param.name), 1)
      return(c(priors, initial))
      }
}

#Prior for RANDOM EFFECTS Term
lmPredRandom2PRIOR<- function(termparen, RHSname, factors, index.name, prior, dropbase, int.check){

  priors = list()

  factors.name <- lapply(factors, removeIndexing)
  factors.name <- unlist(lapply(factors.name, '[[',2))

  term <- removeParen(termparen)
  groups <- no.groups(factors[which(term[2] %in% factors.name)])
  param.name <- as.name(removeIndexing(RHSname))

  #Priors for Random Intercept
  if(term[1] == 1){
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
    priors <- embedLinesInCurlyBrackets(list(random.loop, prior.sd, prior.mu))
    initial <- list(1,1)
    initial <- setNames(initial, c(as.character(make.sigma.name(param.name)), as.character(make.mean.name(param.name))))
    #initial[[1]] <-cbind(as.character(make.sigma.name(param.name)), 1, as.character(make.mean.name(param.name)), 1 )
    return(c(priors, initial))
  }

  else{
    random <- make.random.prior(LHS2BUGSterm(param.name, index.name), mean = make.mean.name(param.name),
                                sig = make.sigma.name(param.name))
    prior.mu <- make.prior(make.mean.name(param.name), 0)
    #prior.mu <- make.prior(make.mean.name(param.name),"Normal")

    random.loop <- embedLinesInForLoop(random, index.name, start = 1, finish = groups)

    #Priors
    prior.sd <-  make.prior(make.sigma.name(param.name),"Uniform")
    priors <- embedLinesInCurlyBrackets(list(random.loop, prior.sd, prior.mu))
    initial <- list(1,1)
    initial <- setNames(initial, c(as.character(make.sigma.name(param.name)), as.character(make.mean.name(param.name))))
    #initial[[1]] <- c(as.character(make.sigma.name(param.name)), 1, as.character(make.mean.name(param.name)), 1)
    return(c(priors, initial))
  }
}


## Make lmPredictor function ##
# Returns Full BUGS code and List of Parameters and Initital Values for MCMC
lmPredictor <- makeBUGSmodule(
  function(LHS, RHS) {

    RHSargs <- match.call(function(RHSmodel, factors, priors, family, link, level, dropbase){}, RHS)

    link <- match.arg(as.character(RHSargs$link), c("identity", "log", "logit", "Poisson", "Bernoulli"))

    LHS.Index <- bracket(as.name(removeIndexing(LHS)[[2]]), quote(i))


    #LHS with Link
    if(link == 'identity') {
      LHS.out <- LHS.Index
    }

    if(link == 'log') {
      LHS.out <- substitute(LINK(LHS), list(LINK = quote(log), LHS = LHS.Index))
    }

    if(link == 'logit') {
      LHS.out <- substitute(LINK(LHS), list(LINK = quote(logit), LHS = LHS.Index))
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
    int.check <- sum(!(is.na(sapply(terms.list, function(x) pmatch("1", x)))))


    #Generate Terms and Priors for Fixed Terms
    RHS.fix <- list()
    Priors.fix <- list()

    for(iTerm in seq_along(RHS.fix.tm)){
      RHS.fix[[iTerm]] <- lmPred2BUGSterm(RHS.fix.tm[[iTerm]], RHSargs$factors, index.name, RHSargs$level)
      Priors.fix[[iTerm]] <- lmPredTerm2PRIOR(RHS.fix.tm[[iTerm]], RHS.fix[iTerm], factors = RHSargs$factors, index.name, prior = RHSargs$priors,
                                             dropbase = RHSargs$dropbase, int.check = int.check)
    }

    #Generate Terms and Priors for Random Terms
    RHS.random <- list()
    Priors.random <- list()
    for(iTerm in seq_along(RHS.random.tm)){
      RHS.random[[iTerm]] <- lmPred2BUGSRandom(RHS.random.tm[[iTerm]], RHSargs$factors, index.name, RHSargs$level)
      Priors.random[[iTerm]] <- lmPredRandom2PRIOR(RHS.random.tm[[iTerm]], RHS.random[iTerm], factors = RHSargs$factors, index.name, prior = RHSargs$priors,
                                                 dropbase = RHSargs$dropbase, int.check = int.check)
    }


    Priors.Combo <- append(lapply(Priors.fix, '[[', 1), lapply(Priors.random, '[[', 1))
    Terms.Combo <- append(RHS.fix, RHS.random)

    #Initial.Values <- append(lapply(Priors.fix, '[[', 2), lapply(Priors.random, '[[', 2))
    Inits.Fix <- sapply(Priors.fix, function(x) x[2])
    Inits.Random <- append(sapply(Priors.random, function(x) x[2]), sapply(Priors.random, function(x) x[3]))
    Inits.All <- append(Inits.Fix, Inits.Random)

    RHS2 <- Terms.Combo[[1]]
    for (i in seq_along(Terms.Combo)[-1]){
      RHS2 <- addTerm(RHS2, Terms.Combo[[i]])
    }

    #Make Model Form and Embed in Loop
    fullLine <- substitute(LHS <- RHS, list(LHS = LHS.out, RHS = RHS2))
    forLoop <- embedLinesInForLoop(fullLine, index.name, start = 1, finish = LHS[[3]][[3]])
    Priors.All <- embedLinesInCurlyBrackets(Priors.Combo)

    newCode <- embedLinesInCurlyBrackets(list(forLoop, Priors.All))

    return(list(code = newCode, inits = Inits.All))
  })


## Functions to Get Code, Data, and Initial Values ready for NIMBLE ##

#Make Starting Values for a Variable Passed To nimbleModel
make.init <- function(term, initialvalue, size){
  l = rep(initialvalue, size)
  init <- substitute(c(TERM,L), list(TERM = as.character(term), L = l))
  init
}


lmPred2init <- function(term, value){
    out = rep(value, term)
    term <- out
    return(term)
  }

## Functions that Interact with NIMBLE ##

#Call Nimble and Pass the BUGS code and data
callNim <- function(modobj, niter, burn, chains){
  nimMod.C <- compileNimble(modobj)

  #Set Up MCMC
  config.Mod <- configureMCMC(nimMod.C, print = TRUE)
  mod.MCMC <- buildMCMC(config.Mod)
  C.mod.MCMC <- compileNimble(mod.MCMC)
  samplesList <- runMCMC(C.mod.MCMC, niter = niter, nburnin = burn, nchains = chains, returnCodaMCMC = TRUE)
  return(samplesList)
}



#Take the MCMC samples and return summary statistics
makeoutput <- function(codaobj, sitevars, obsvars){
  quantile <- summary(codaobj)$quantiles
  meansd <- summary(codaobj)$statistics[,1:2]

  #Back Transform Site and Obs Estimates
  #out.names <- names(quantile[,1])
  #site.match <- exp(quantile[match(sitevars, out.names), drop = FALSE,])
  #obs.match <-  exp(quantile[match(obsvars, out.names), drop = FALSE,])/(1 + exp(quantile[match(obsvars, out.names), drop = FALSE,]))
  #rownames(site.match) <- paste0("exp.", sitevars)
  #rownames(obs.match) <- paste0("inv.logit.", obsvars)
  #quantile.adj <- rbind(obs.match, site.match)


  rhat <- tryCatch(gelman.diag(codaobj)$psrf, error = function(e) NA)
  Effsamp <- effectiveSize(codaobj)
  output  <- as.data.frame(cbind(meansd, quantile, rhat, Effsamp))
  names(output)[names(output) == 'Point est.'] <- 'Rhat'
  output$`Upper C.I.`<-NULL
  return(list(output))
}


## Main Function for Abundance N Mixture Models ##
nimble.abund <- function(siteformula = NULL, obsformula = NULL, y = NULL, sitevars = NULL, obsvars = NULL, mixture = c("Poisson", "ZIP", "NB"),
                         priors = c("Normal","t","Uniform"), dropbase = TRUE, niter = 10000, burnin = 1000, initmcmc = 1, chains = 1,
                         returncode = FALSE, returnsamp = FALSE){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)

  obsform <- as.formula(mf$obsformula)
  stateform <- as.formula(mf$siteformula)

  df <- maketidy(y, sitevars, obsvars) #turn 3 data pieces into one tiny data frame


  #Length of Long DF
  L <- as.numeric(dim(df)[1])

  #No. Sites
  S <- as.numeric(dim(y)[1])

  #Convert Character Variables to Factors
  for (i in 1:dim(df)[2]){
    df[,i] <- makeFac(df[,i], char.only = TRUE)
  }

  # Make Site Level Formula #
  LHS.Site <- substitute(N[1:L], list(L = S))
  factors <- names(Filter(is.factor, df))
  factors.size <- lapply(factors, factor.bracket, df)
  RHS.Site <- make.glm(mf$siteformula, factors.size,  cl = cl, mixture, "log", level = quote(site))

  LHS.Obs <- substitute(y[1:S], list(S=L))
  RHS.Obs <- make.glm(mf$obsformula, factors.size,  cl= cl, "Binomial", "log", N = quote(N), Site = quote(Site), level = quote(obs))

  glm.expand.Site <- nim_glm$process(LHS.Site, RHS.Site)
  glm.expand.Obs <- nim_glm$process(LHS.Obs, RHS.Obs)


  lm.expand.Site <- lmPredictor$process(glm.expand.Site$LHS, glm.expand.Site$RHS)
  lm.expand.Obs <- lmPredictor$process(glm.expand.Obs$LHS, glm.expand.Obs$RHS)

  #Full BUGS Code Expansion
  full.expand <- embedLinesInCurlyBrackets(lines = list(glm.expand.Site$prob.mod, lm.expand.Site$code, glm.expand.Obs$prob.mod, lm.expand.Obs$code))

  #Get Data and Model Readyfor nimbleModel and MCMC

  #Make Factor Variables Numeric
  indx <- sapply(df, is.factor)
  df[indx] <- lapply(df[indx], function(x) as.numeric(as.factor(x)))

  start.values <- lapply(append(lm.expand.Site$inits, lm.expand.Obs$inits), lmPred2init, initmcmc)

  #Inititalize N, the latent counts
  N <- apply(y, 1, max);  N[is.na(N)] <- 0 ; N = N + 1
  start.values <- append(start.values, list(N = N))

  if(mixture == "NB"){
    start.values <- append(start.values, list(logalpha = 1))
  }

  if(mixture == "ZIP"){
    start.values <- append(start.values, list(theta = .5))
  }

  nimMod.obj <- nimbleModel(code = full.expand, inits = start.values, constants = as.list(df), data = list(y = df$Count))
  browser()
  mod.output <- callNim(nimMod.obj, niter = niter, burn = burnin, chains = chains)

  if(is.list(mod.output)){mod.out = do.call(rbind, mod.output)} else{mod.out = mod.output}

  results <- list("Summary" = makeoutput(mod.output, names(lm.expand.Site$inits), names(lm.expand.Obs$inits)))   #Q: Do we want to backtransform parameters from log and logit scale?


  if(returncode==TRUE){results[["BUGScode"]] = full.expand}
  if(returnsamp==TRUE){results[["Samples"]] = mod.out}
  return(results)
}


## Main Function for Occupancy Models ##
nimble.occup <- function(siteformula = NULL, obsformula = NULL, y = NULL, sitevars = NULL, obsvars = NULL, mixture = c("Bernoulli"),
                         priors = c("Normal","t","Uniform"), dropbase = TRUE, niter = 10000, burnin = 1000, initmcmc = 1, chains = 1,
                         returncode = FALSE, returnsamp = FALSE){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)

  obsform <- as.formula(mf$obsformula)
  stateform <- as.formula(mf$siteformula)

  df <- maketidy(y, sitevars, obsvars) #turn 3 data pieces into one tiny data frame


  #Length of Long DF
  L <- as.numeric(dim(df)[1])

  #No. Sites
  S <- as.numeric(dim(y)[1])

  #Convert Character Variables to Factors
  for (i in 1:dim(df)[2]){
    df[,i] <- makeFac(df[,i], char.only = TRUE)
  }

  # Make Site Level Formula #
  LHS.Site <- substitute(z[1:L], list(L = S))
  factors <- names(Filter(is.factor, df))
  factors.size <- lapply(factors, factor.bracket, df)
  RHS.Site <- make.glm(mf$siteformula, factors.size,  cl, mixture, "logit", level = quote(site))

  LHS.Obs <- substitute(y[1:S], list(S=L))
  RHS.Obs <- make.glm(mf$obsformula, factors.size,  cl, mixture, "logit", N = quote(N), Site = quote(Site), level = quote(obs))
  browser()

  glm.expand.Site <- nim_glm$process(LHS.Site, RHS.Site)
  glm.expand.Obs <- nim_glm$process(LHS.Obs, RHS.Obs)

  lm.expand.Site <- lmPredictor$process(glm.expand.Site$LHS, glm.expand.Site$RHS)
  lm.expand.Obs <- lmPredictor$process(glm.expand.Obs$LHS, glm.expand.Obs$RHS)

  #Full BUGS Code Expansion
  full.expand <- embedLinesInCurlyBrackets(lines = list(glm.expand.Site$prob.mod, lm.expand.Site$code, glm.expand.Obs$prob.mod, lm.expand.Obs$code))

  #Get Data and Model Readyfor nimbleModel and MCMC

  #Make Factor Variables Numeric
  indx <- sapply(df, is.factor)
  df[indx] <- lapply(df[indx], function(x) as.numeric(as.factor(x)))

  start.values <- lapply(append(lm.expand.Site$inits, lm.expand.Obs$inits), lmPred2init, initmcmc)

  #Inititalize N, the latent counts
  z <- apply(y, 1, max);  z[is.na(z)] <- 0
  start.values <- append(start.values, list(z = z))

  nimMod.obj <- nimbleModel(code = full.expand, inits = start.values, constants = as.list(df), data = list(y = df$Count))
  browser()
  mod.output <- callNim(nimMod.obj, niter = niter, burn = burnin, chains = chains)

  if(is.list(mod.output)){mod.out = do.call(rbind, mod.output)} else{mod.out = mod.output}

  results <- list("Summary" = makeoutput(mod.output, names(lm.expand.Site$inits), names(lm.expand.Obs$inits)))   #Q: Do we want to backtransform parameters from log and logit scale?


  if(returncode==TRUE){results[["BUGScode"]] = full.expand}
  if(returnsamp==TRUE){results[["Samples"]] = mod.out}
  return(results)
}




