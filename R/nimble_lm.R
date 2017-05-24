## Nimble LM Model##
## BUGS Code Generation for nimbleCode Object for Linear Model ##

library(assertthat) ## to package dependencies

## Helper Functions ##

## This Function Makes the Name someFactor.effects
make.effect.name <- function(term) {
  paste0(term, '.effect')
}

#Make Factor Variable
makeFac <- function(x,char.only=FALSE) {
  if (!is.factor(x) && (!char.only || is.character(x))) factor(x) else x
}




#Make Normal Diffuse Prior N ~ (0, 100,000), Uniform (0,1000), t(0,1,d), or point mass at 0
#We Could Also Allow an Argument to Specify the Shape and Scale Parameter?
make.prior <- function(term,prior,...){

  if (prior=="Normal"){
  substitute(TERM ~ dnorm(0,0.0001), list(TERM=term))
  }

  else if (prior=="Uniform"){
    substitute(TERM ~ dunif(0,1000), list(TERM=term))
  }

  else if (prior=="t"){
    substitute(TERM ~ dt(0,1,d), list(TERM=term))
  }

  else{
    substitute(TERM <- 0, list(TERM=term))
    }
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



##Converts a Term from a Formula to a BUGS Term
formulaTerm2BUGSterm <- function(term, mf, index.name,id) {
  ## term will be a string like "x"
  ## mf will be the model.frame
  ## index.name will be the name of the for-loop index

  ## Need to Account for Interaction Terms.  Represented by : in attr(terms)
  ## Converts Interaction term to BUGS code, e.g. x[i]*y[i]; Need to Expand to take greater than 2 way interactions
  if (length(plyr::ldply(strsplit(term, split = ":"))) > 1){
    inter = plyr::ldply(strsplit(term, split = ":"))
    beta.no = as.name((paste0(quote(b),id)))
        if(is.factor(mf[[inter[,1]]])){
          beta.no = as.name((paste0(quote(b),id)))
          out = substitute(EFFECTNAME[TERM[INDEX]]*TERM2[INDEX],
                           list(EFFECTNAME = as.name(beta.no), #Index Beta's on Covariates
                                TERM = as.name(inter[,1]),
                                TERM2 = as.name(inter[,2]),
                                INDEX = index.name))
          return(out)
        }
        if(is.factor(mf[[inter[,2]]])){
        beta.no = as.name((paste0(quote(b),id)))
        out = substitute(EFFECTNAME[TERM[INDEX]]*TERM2[INDEX],
                       list(EFFECTNAME = as.name(beta.no),
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
    return(make.param.name(beta.no,out))
  }

  if(is.factor(mf[[term]])){
    ## If term is a factor use BUGS indexing create b2[block[i]]
    beta.no = as.name((paste0(quote(b),id)))
    out = substitute(EFFECTNAME[TERM[INDEX]],  #How to Keep double Brackets without Term Name in Front?
               list(EFFECTNAME = as.name(beta.no), #was (effect.name)
                    TERM = as.name(term),
                    INDEX = index.name))
    return(out)
  }
  else
    ## otherwise create BUGS term like x[i]
    out = substitute(TERM[INDEX],
               list(TERM = as.name(term),
                    INDEX = index.name))

    beta.no = as.name((paste0(quote(b),id)))
    make.param.name(beta.no,out)

}


## This Adds a Prior for Each Parameter on the RHS

formulaTerm2PRIOR <- function(term, mf, index.name, id, dropbase, prior) {
  ## term will be a string like "x"
  ## mf will be the model.frame
  ## index.name will be the name of the for-loop index
  ##Would be nice (maybe easier?) to Parse apart the RHS language term and define priors from that object?

  if(is.factor(mf[[term]])){
    ## if term is a factor, create priors for each level (note, dropbase=T will assign 0 prior to first level)
    out = list()  #Add Each Prior to List Object

    groups = nlevels(as.factor(mf[[term]]))  #Drop one group for model identification; Need to add Check if Intercept in Model Too
    beta.no = as.name(paste0(quote(b),id))

    if(dropbase==TRUE){
      out[[1]] <- make.prior(substitute(beta.no[ID],list(ID=1,beta.no=beta.no)),prior=0)
      st = 2
      idx = 2
      }
    else {st=1; idx=1}

    p = make.prior(LHS2BUGSterm(beta.no,index.name), prior)
    out[[idx]] <- embedLinesInForLoop(p, index.name, start=st, finish=as.numeric(groups))
    return(embedLinesInCurlyBrackets(out))
  }



   #Priors for Interaction Terms with Factor Variable
   if (length(plyr::ldply(strsplit(term, split = ":"))) > 1){
     inter = plyr::ldply(strsplit(term, split = ":"))
     beta.no = as.name((paste0(quote(b),id)))
     if(is.factor(mf[[inter[,1]]]) | is.factor(mf[[inter[,2]]])){

       out = list()  #Add Each Prior to List Object

       groups = max((nlevels(mf[[inter[,1]]])), nlevels(mf[[inter[,2]]]))  #Allow More Than 3  Way Interaction?
       beta.no = as.name(paste0(quote(b),id))

       if(dropbase==TRUE){
         out[[1]] <- make.prior(substitute(beta.no[ID],list(ID=1,beta.no=beta.no)),prior=0)
         st = 2
         idx = 2
       }
       else {st=1; idx=1}

       p = make.prior(LHS2BUGSterm(beta.no,index.name), prior)
       out[[idx]] <- embedLinesInForLoop(p, index.name, start=st, finish=as.numeric(groups))

       return(embedLinesInCurlyBrackets(out))

     }

   }

   else
     beta.no = as.name((paste0(quote(b),id)))
     out = make.prior(beta.no, prior)
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

  ## change name of function
  mf[[1L]] <- quote(stats::model.frame)

  ## and eval it in calling environment
  mf <- eval(mf, parent.frame())
  ## This gets us mf with data filled in from dat argument or from calling environment as needed

  #Convert Character Variables to Factors
  for (i in 1:dim(mf)[2]){
    mf[,i] <- makeFac(mf[,i],char.only = TRUE)
  }


  #Check link and prior arguments are supported
  family = match.arg(family)
  prior = match.arg(priors)

  mod = as.formula(mod)

  index.name <- quote(i)

  #Likelihood
  if (family=="Normal"){
    meanfcn <- quote(dnorm(mu[i],tau))}  #could change to other disributions later on
    dep <- LHS2BUGSterm(mod[[2]], index.name)
    meanfctn <-  substitute(LHS ~ RHS, list(LHS = dep, RHS = meanfcn))


    #Generate mean function mu[i]
    #Intercept
    int.to.add <- attr(terms(mod), 'intercept')

    terms.to.add <- attr(terms(mod), 'term.labels')
    assert_that(length(terms.to.add) > 0)


  ## Assume Fixed Intercept for Now
  if (int.to.add==1){
      RHS <- quote(b0)
      ## append additional terms
      for(iTerm in seq_along(terms.to.add)) {
        RHS <- addTerm(RHS, formulaTerm2BUGSterm(terms.to.add[iTerm], mf, index.name,id = iTerm))
    }
  }

  if(int.to.add==0){
    RHS <- formulaTerm2BUGSterm(terms.to.add[1], mf, index.name, id=1)
    for(iTerm in seq_along(terms.to.add[-1])+1) {
      RHS <- addTerm(RHS, formulaTerm2BUGSterm(terms.to.add[iTerm], mf, index.name,id = iTerm))
    }
  }


  #Make Term for mean function, mu
  LHS <- LHS2BUGSterm(quote(mu), index.name)
  fullLine <- substitute(LHS ~ RHS, list(LHS = LHS, RHS = RHS))
  lines <- list(meanfctn, fullLine) ## could have potentially more
  forLoop <- embedLinesInForLoop(lines, index.name, start = 1, finish = as.numeric(dim(mf)[1])) ## num.data will have to be figured out


  #Specify Priors


  priorsList <- list()

  #Prior for Variance Term, could allow options for this later
  priorsList[[1]] <-  quote(tau ~ dgamma (0.001, 0.001))
  priorsList[[2]] <-  quote(sigma2 <- 1/tau)


  #Priors for Beta Terms
  #Have a List of Index IDs to loop Over.  Is using same index for loops a good idea?
  for(iTerm in seq_along(terms.to.add)) {
    priorsList[[iTerm+2]] <- formulaTerm2PRIOR(terms.to.add[iTerm], mf, index.name,id = iTerm, dropbase, prior = priors)
  }

  browser()

  #Add Prior for Intercept
  if(int.to.add==1){
    priorsList[[length(priorsList)+1]] <- make.prior(quote(b0),priors)
  }


  Priors <- embedLinesInCurlyBrackets(priorsList)
  bugs.out <- embedLinesInCurlyBrackets(list(forLoop,Priors))

  bugs.out
}


##Test Function

sillyData <- list(z = rpois(10,2), q=rnorm(10), xm = rnorm(10), ym=rgamma(10,3), A = c(rep('Red',7), rep('Orange',1), rep('Blue', 2)))

#debug(nimble.lm)

test <- nimble.lm(mod = z ~ xm + ym*A, dat = sillyData, priors="Uniform", dropbase=TRUE)

test

## suggest writing unit tests with testthat hand-in-hand with code development

