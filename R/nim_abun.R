### NIMBLE ABUNDANCE ###
library(nimble)
nimbleOptions(enableBUGSmodules = TRUE)

#Make Bracket from Term with Indexing Running to Either No. Obs or No. Groups if Factor Variable
mk.bracket <- function(term, sitemf, obslist, M, J){
  if (term %in% names(sitemf)==TRUE){

      if(is.factor(sitemf[[term]])){
      l = as.numeric(nlevels((sitemf[[term]])))
    }
    else{
      l = as.numeric(length(sitemf[[term]]))
    }

    return(substitute(TERM[1:L], list(TERM = as.name(term), L = l)))
  }

  if (term %in% names(obslist)==TRUE){
    if(is.factor(obslist[[term]])){
      l = as.numeric(nlevels((obslist[[term]])))  #Need to Figure Out how to Determine Factors from a List of Variables
    }
    else{
      l = M
    }
    return(substitute(TERM[1:M, 1:J], list(TERM = as.name(term), M = l, J = J)))
  }

  else{stop('A Term was Supplied that is Not in the Data')}
}


bracket <- function(term, sitemf, obslist, M, J){
  if (grepl(":",term)==TRUE){
    terms.split <- unlist(strsplit(term, ":"))
    t.bracket = lapply(terms.split, mk.bracket, sitemf, obslist, M, J)
    out <- prodTerm(t.bracket[[1]], t.bracket[[2]])

    ## append additional terms
    for(i in seq_along(t.bracket)[-c(1:2)]){
      out <- prodTerm(out, t.bracket[[i]])
    }
    out
  }
  else{
    out = mk.bracket(term, sitemf, obslist, M, J)
    out
  }
}

#Same Bracket Function for Random Effects |
bracket.random <- function(term, sitemf, obslist, M, J){
  if(term[[2]]==1){
    lside = 1}
  else{
    lside <- mk.bracket(term[[2]], sitemf, obslist, M, J)
  }
  rside <- mk.bracket(as.character(term[[3]]), sitemf, obslist, M, J)
  out <- substitute((L|R), list(L=lside, R=rside))
  out
}

#Make nim_glm object to pass to BUGS Module
make.glm.rhs <- function(RHS, sitemf, obslist, cl, mixture, link){
  RHS.State <- call("nim_glm", RHS)
  if(length(names(grep("factor", sapply(sitemf, class), value=TRUE))) > 0)
    {RHS.State[[3]] <- names(grep("factor", sapply(sitemf, class), value=TRUE))}
  else{
    RHS.State[[3]] <- "None"
  }
  RHS.State[[4]] <- cl$priors
  RHS.State[[5]] <- cl$dropbase
  RHS.State[[6]] <- mixture
  RHS.State[[7]] <- link
  names(RHS.State)[3:7] <- c("factors", "priors","dropbase", "family", "link")
  return(RHS.State)
}

#Make RHS with Brackets for Fixed Terms
make.RHS.bracket <- function(formula, terms, siteCov, ObsCov, M, J){
  if(length(terms)==0 & attr(terms(formula), 'intercept')==0){return()}
  else{
    if (attr(terms(formula), 'intercept')==1){
    RHS <- substitute(TERM[1:M], list(TERM = quote(intercept), M = M))
    #Append Additional Terms
    for(iTerm in seq_along(terms)) {
      RHS <- addTerm(RHS, bracket(as.character(terms[iTerm]), siteCov, ObsCov, M, J ))
    }
  }

  if (attr(terms(formula), 'intercept')!=1  & length(terms > 0)){
    RHS <- bracket(terms[1], siteCov, ObsCov, M, J)
      for(iTerm in seq_along(terms[-1]) + 1){
        RHS <- addTerm(RHS, bracket(as.character(terms[iTerm]), siteCov, ObsCov, M, J))
    }
  }
  }
    return(RHS)
}


#Add Random Terms
make.RHS.bracket.random <- function(terms, siteCov, ObsCov, M, J, RHS){
if (length(terms) > 0 & !is.null(RHS)){
  for (iTerm in seq_along(terms)){
    RHS <- addTerm(RHS, bracket.random(as.character(terms[[iTerm]]), siteCov, ObsCov, M, J))
  }
}

if(length(terms) > 0 & is.null(RHS)){
  RHS <- bracket.random(terms[[1]], siteCov, ObsCov, M, J)
  for (iTerm in seq_along(terms[-1]) + 1){
    RHS <- addTerm(RHS, bracket.random(as.character(terms[[iTerm]]), siteCov, ObsCov, M, J))
  }
}
  return(RHS)
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




nimble.abund <- function(formula, y = NULL, siteCov = NULL, ObsCov = NUll, mixture = c("Poisson", "ZIP", "NB"), priors = c("Normal","t","Uniform"),
                         dropbase = TRUE, niter = 10000, burnin = 1000, initmcmc = 1, chains = 1, returncode = FALSE, returnsamp = FALSE){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)

  obsformula <- as.formula(mf$formula[[2]])
  stateformula <- as.formula(paste("~", mf$formula[3], sep=""))

  #Get Number of Sites (M) and Number of Visits (J)
  M = as.numeric(dim(y)[1])
  J = as.numeric(dim(y)[2])

  # Make Site Level Formula #
  LHS.Site = substitute(N[1:L], list(L = M))
  fixedTm <- nobars(stateformula)
  terms.to.add <- attr(terms(fixedTm), 'term.labels')

  #Extract Random Terms
  randomTm <- fb(stateformula)

  #Add RHS Terms and Create nim_glm object for site level
  RHS.fix <- make.RHS.bracket(fixedTm, terms.to.add, siteCov, ObsCov, M, J)
  RHS.all <- make.RHS.bracket.random(randomTm, siteCov, ObsCov, M, J, RHS.fix)
  RHS.Site <- make.glm.rhs(RHS.all, siteCov, ObsCov, cl,"Poisson", "log")

  # Make Obs Level Formula #
  LHS.Obs = substitute(y[1:M, 1:J], list(M=M, J=J))
  fixedTmObs <- nobars(obsformula)
  terms.to.add <- attr(terms(fixedTmObs), 'term.labels')

  #Extract Random Terms
  randomTmObs <- fb(obsformula)

  #Add RHS Terms and Create nim_glm object for obs level
  RHS.fix.obs <- make.RHS.bracket(fixedTmObs, terms.to.add, siteCov, ObsCov, M, J)
  RHS.all.obs <- make.RHS.bracket.random(randomTmObs, siteCov, ObsCov, M, J, RHS.fix.obs)

  RHS.Obs <- make.glm.rhs(RHS.all.obs, siteCov, ObsCov, cl, "Binomial", "logit")

  return(list(RHS.Obs = RHS.Obs, RHS.Site = RHS.Site))
  #Next Step is Module Processing of These Two Pieces
}





