## Nimble LM Model##
## BUGS Code Generation for nimbleCode Object for Linear Model ##
library(nimble)
library(assertthat) ## to package dependencies

## Helper Functions ##

## These Functions Make the Name someFactor.X
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


#Functions for Random Effects

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
    if((term[[1]] == as.name("|"))) {
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


## This will take a list of code lines and embed them in {}
## Would be nice to have function that embeds lines without any syntax. BUGS should ignore extra {} so may be OK for now.
embedLinesInCurlyBrackets <- function(lines) {
  as.call(c(list(quote(`{`)), lines))
}


## This will take a list of BUGS lines and embed them in a for loop
## In future we may want to call this repeatedly to create nested loops
## or modify it to create nested loops in one go
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



##Converts a Term from a Formula to a BUGS Term
formulaTerm2BUGSterm <- function(term, mf, index.name) {
  ## term will be a string like "x"
  ## mf will be the model.frame
  ## index.name will be the name of the for-loop index

  ## Need to Account for Interaction Terms.  Represented by : in attr(terms)
  ## Converts Interaction term to BUGS code, e.g. x[i]*y[i]; Need to Expand to take greater than 2 way interactions
  if (length(plyr::ldply(strsplit(term, split = ":"))) > 1){
    inter = plyr::ldply(strsplit(term, split = ":"))
        if(is.factor(mf[[inter[,1]]])){
          out = substitute(EFFECTNAME[TERM[INDEX]]*TERM2[INDEX],
                           list(EFFECTNAME = make.coef.name(paste0(inter$V1,inter$V2)), #Index Beta's on Covariates
                                TERM = as.name(inter[,1]),
                                TERM2 = as.name(inter[,2]),
                                INDEX = index.name))
          return(out)
        }
        if(is.factor(mf[[inter[,2]]])){
        out = substitute(EFFECTNAME[TERM[INDEX]]*TERM2[INDEX],
                       list(EFFECTNAME = make.coef.name(paste0(inter$V2,inter$V1)),
                            TERM = as.name(inter[,2]),
                            TERM2 = as.name(inter[,1]),
                            INDEX = index.name))
      return(out)
      }

      else{
              out = substitute(TERM[INDEX]*TERM2[INDEX],
                     list(TERM = as.name(inter[,1]),
                     TERM2 = as.name(inter[,2]),
                     INDEX = index.name))}
    return(make.param.name(make.coef.name(paste0(inter$V1,inter$V2)),out))
  }

  if(is.factor(mf[[term]])){
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
formulaTerm2PRIOR <- function(term, mf, index.name, dropbase, prior) {
  ##Would be nice (maybe easier?) to Parse apart the RHS language term and define priors from that object?

  if(is.factor(mf[[term]])){
    ## if term is a factor, create priors for each level (note, dropbase=T will assign 0 prior to first level)
    out = list()  #Add Each Prior to List Object

    groups = nlevels(as.factor(mf[[term]]))  #Drop one group for model identification; Need to add Check if Intercept in Model Too
    #beta.no = as.name(paste0(quote(b),id))

    if(dropbase==TRUE){
      out[[1]] <- make.prior(LHS2BUGSterm(make.coef.name(term), 1), prior=0)
      st = 2
      idx = 2
      }
    else {st=1; idx=1}

    p = make.prior(LHS2BUGSterm(make.coef.name(term), index.name), prior)
    out[[idx]] <- embedLinesInForLoop(p, index.name, start=st, finish=as.numeric(groups))
    return(embedLinesInCurlyBrackets(out))
  }

   #Priors for Interaction Terms with Factor Variable
   if (length(plyr::ldply(strsplit(term, split = ":"))) > 1){
     inter = plyr::ldply(strsplit(term, split = ":"))
     #beta.no = as.name((paste0(quote(b),id)))
     if(is.factor(mf[[inter[,1]]]) | is.factor(mf[[inter[,2]]])){

       out = list()  #Add Each Prior to List Object

       groups = max((nlevels(mf[[inter[,1]]])), nlevels(mf[[inter[,2]]]))  #Allow More Than 3  Way Interaction?
       #beta.no = as.name(paste0(quote(b),id))

       if(dropbase==TRUE){
         out[[1]] <- make.prior(LHS2BUGSterm(make.coef.name(paste0(inter$V1,inter$V2)), 1), prior=0)
         st = 2
         idx = 2
       }
       else {st = 1; idx = 1}

       p = make.prior(LHS2BUGSterm(make.coef.name(paste0(inter$V1,inter$V2)), index.name), prior)
       out[[idx]] <- embedLinesInForLoop(p, index.name, start=st, finish=as.numeric(groups))

       return(embedLinesInCurlyBrackets(out))
      }


   else {
     out = make.prior(make.coef.name(paste0(inter$V1,inter$V2)), prior)
        }
   }
   else{
     out = make.prior(make.coef.name(term), prior)
   }

  return(out)
}


## Main nimble.lm Function ##
## Future Options: shape/scale input for priors, additional family link functions, random effects
##  dropbase will drop lowest group for cateogrical variables.

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
  #environment(fr.form) <- environment(formula)


  mf$formula <- fr.form
  fr <- eval(mf, parent.frame()) # Gets us mf with data filled in from dat argument or from calling environment as needed


  #Convert Character Variables to Factors
  for (i in 1:dim(fr)[2]){
    fr[,i] <- makeFac(fr[,i],char.only = TRUE)
  }


  #Check link and prior arguments are supported
  family = match.arg(family)
  prior = match.arg(priors)


  index.name <- quote(i)


  #Likelihood

  #Generate mean function mu[i]
  if (family=="Normal"){
    meanfcn <- substitute(dnorm(MEAN, sd=SIGMA), list(MEAN=LHS2BUGSterm(make.mean.name(RHSForm(mf$formula[2])), index.name),
                SIGMA=quote(sigma)))}  #could change to other disributions later on
    dep <- LHS2BUGSterm(make.pred.name(mf$formula[[2]]), index.name)
    meanfctn <-  substitute(LHS ~ RHS, list(LHS = dep, RHS = meanfcn))

    #Pull out Fixed Terms
    fixedTm <- nobars(formula)
    terms.to.add <- attr(terms(fixedTm), 'term.labels')

    randomTerms <- fb(formula)  #Puts Random Effect Terms in a List

    #Store Random Effect Terms and Priors
    RandomParms <- list()
    RandomPriors <- list()


    #Generate Random Terms and Priors First
    if(!is.null(randomTerms)){
      for (i in 1:length(randomTerms)){
        RandomParms[[i]] <- formulaRandom2BUGS(randomTerms[[i]], fr, index.name, intcpt = int.add, fixedterms = terms.to.add)[[1]]
        RandomPriors[[i]] <- formulaRandom2BUGS(randomTerms[[i]], fr, index.name, intcpt = int.add, fixedterms = terms.to.add)[[2]]
      }
    }

  #Add Fixed Effect Terms
  assert_that(length(terms.to.add) > 0 | length(randomTerms) > 0)

  ## Check for Overal Model Intercept
  if (int.add==1){
      RHS <- quote(b0)
      ## append additional terms
      for(iTerm in seq_along(terms.to.add)) {
        RHS <- addTerm(RHS, formulaTerm2BUGSterm(terms.to.add[iTerm], fr, index.name))
    }
  }

  if(int.add==0 & length(terms.to.add) > 0){
    RHS <- formulaTerm2BUGSterm(terms.to.add[1], fr, index.name)
    for(iTerm in seq_along(terms.to.add[-1]) + 1){
      RHS <- addTerm(RHS, formulaTerm2BUGSterm(terms.to.add[iTerm], fr, index.name))
    }
  }

  #Add Back in Random Terms; Need Check if There are No Fixed Terms Or Intercept Though to Define RHS
  if(!is.null(randomTerms) & exists("RHS")){
    for (i in 1:length(RandomParms)){
      RHS <- addTerm(RHS, RandomParms[[i]])
    }
  }

  if(!is.null(randomTerms) & !exists("RHS")){
    RHS <- RandomParms[[1]]
    for (i in seq_along(RandomParms[-1]) + 1){
      RHS <- addTerm(RHS, RandomParms[[i]])
    }
  }

  #Make Term for mean function, mu
  LHS <- LHS2BUGSterm(make.mean.name(RHSForm(mf$formula[2])), index.name)
  fullLine <- substitute(LHS <- RHS, list(LHS = LHS, RHS = RHS))
  lines <- list(meanfctn, fullLine) ## could have potentially more
  forLoop <- embedLinesInForLoop(lines, index.name, start = 1, finish = as.numeric(dim(fr)[1])) ## num.data will have to be figured out


  #Specify Priors


  priorsList <- list()

  #Prior for Variance Term, could allow options for this later
  priorsList[[1]] <-  quote(sigma ~ dunif(0, 1000))

  #Priors for Beta Terms
  #Have a List of Index IDs to loop Over.  Is using same index for loops a good idea?
  for(iTerm in seq_along(terms.to.add)) {
    priorsList[[iTerm + 1]] <- formulaTerm2PRIOR(terms.to.add[iTerm], fr, index.name, dropbase, prior = priors)
  }

  #Add in Random Effect Priors
  if (length(RandomPriors) > 0){
    priorsList[[length(priorsList) + 1]] <- embedLinesInCurlyBrackets(RandomPriors)
  }

  #Add Prior for Intercept
  if(int.add==1){
    priorsList[[length(priorsList) + 1]] <- make.prior(quote(b0), priors)
  }

  Priors <- embedLinesInCurlyBrackets(priorsList)
  bugs.out <- embedLinesInCurlyBrackets(list(forLoop, Priors))

  bugs.out
}


##Test Function

sillyData <- list(z = rpois(10,2), q=rnorm(10), xm = rnorm(10), ym=rgamma(10,3), A = c(rep('Red',7), rep('Orange',1), rep('Blue', 2)))

#debug(nimble.lm)
test <- nimble.lm(mod = z ~ xm + A, dat = sillyData, priors="Normal", dropbase=TRUE)

test

##Run in NIMBLE to Check if There are NAMESPACE ISSUES

inits <- list(sigma=1, b.xm = 1, b.A=c(1,1,1), b0 = 1)

dat <- list(z.pred=sillyData$z)

#Make A a Numeric
A.num <-as.numeric(as.factor(sillyData$A))

cont <- list(xm=sillyData$xm, A=A.num)


silly <- nimbleModel(code=test, name="mod1", inits=inits, constants=cont, data=dat)
silly$plotGraph()
Csilly <- compileNimble(silly)

#Set Up MCMC
sillyConf <- configureMCMC(silly, print = TRUE)
sillyMCMC <- buildMCMC(sillyConf)
CsillyMCMC <- compileNimble(sillyMCMC, project = silly)
niter <- 10000
set.seed(1)
CsillyMCMC$run(niter)

#Look at Samples
samples <- as.matrix(CsillyMCMC$mvSamples)
plot(samples[ , "b.A[3]"], type = "l", xlab = "iteration",
     ylab = expression(sigma))

median(samples[,2])

#Check with lm function
mod <- lm(z ~ xm + A, data=sillyData)
