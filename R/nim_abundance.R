## Main Function for Abundance N Mixture Models ##
nimble.abund <- function(siteformula = ~1, obsformula = ~1, y = NULL, sitevars = NULL, obsvars = NULL, mixture = c("Poisson", "ZIP", "NB"),
                         priors = c("Normal","t","Uniform"), dropbase = TRUE, niter = 10000, burnin = 1000, initmcmc = 1, chains = 1,
                         returncode = FALSE, returnsamp = FALSE){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)

  obsform <- as.formula(mf$obsformula)
  stateform <- as.formula(mf$siteformula)

  if(length(obsform)==0 | length(stateform)==0){stop('A State and Observation Formula Must be Provided')}

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
  RHS.Obs <- make.glm(mf$obsformula, factors.size,  cl= cl, "Binomial", "logit", N = quote(N), Site = quote(Site), level = quote(obs))

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
  mod.output <- callNim(nimMod.obj, niter = niter, burn = burnin, chains = chains)

  if(is.list(mod.output)){mod.out = do.call(rbind, mod.output)} else{mod.out = mod.output}

  results <- list("Summary" = makeoutput(mod.output, names(lm.expand.Site$inits), names(lm.expand.Obs$inits)))   #Q: Do we want to backtransform parameters from log and logit scale?


  if(returncode==TRUE){results[["BUGScode"]] = full.expand}
  if(returnsamp==TRUE){results[["Samples"]] = mod.out}
  return(results)
}

