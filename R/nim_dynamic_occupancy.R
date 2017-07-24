#Main NIMBLE Function for Dynamic Occupancy Models

nimble.dynamic.occ <- function(psiformula = NULL, phiformula = NULL, gammaformula = NULL, pformula = NULL, y = NULL, site = NULL, site.season = NULL,
                               site.season.survey = NULL, priors = c("Normal","t","Uniform"), mixture = "Bernoulli", dropbase = TRUE, niter = 10000,
                               burnin = 1000, initmcmc = .1, chains = 4, returncode = TRUE, returnsamp = TRUE){

  cl <- match.call.defaults()  #Gets All Arguments, Including Defaults
  mf <- match.call(expand.dots = FALSE)

  #Pull of Model Formulas
  psiform <- as.formula(mf$psiformula)
  phiform <- as.formula(mf$phiformula)
  gammaform <- as.formula(mf$gammaformula)
  pform <- as.formula(mf$pformula)

  nsite <- as.numeric(dim(y)[1])
  nsurvey <- as.numeric(dim(y)[2])
  nseason <- as.numeric(dim(y)[3])

  index.name <- quote(i)

  #Make Input Data Tidy
  tidy.data <- tidy.dynam.occ(y, site, site.season, site.season.survey)

  tidy.data.sort = arrange(tidy.data, Survey)   #Allows Looping Over Seasons w/o Skipping Rows with Multiple Surveys Per Season

  tidy.data.sort$SiteSeason <- rep(1:(nsite*nseason), times = nsurvey, each = 1)  #for nested indexing over latent states

  #Identify Factor Variables
  factors <- names(Filter(is.factor, tidy.data))
  factors.size <- lapply(factors, factor.bracket, tidy.data)


  #Ecological SubModels

  #Model for Site
  LHS.Site <- substitute(z[1:L], list(L = nsite))
  RHS.Site <- make.glm(psiform, factors= factors.size , cl, mixture, "logit", level = quote(site))
  glm.expand.Site <- nim_glm$process(LHS.Site, RHS.Site)
  lm.expand.Site <- lmPredictor$process(glm.expand.Site$LHS, glm.expand.Site$RHS)

  #Models for Phi and Gamma
  LHS.Phi <- substitute(z[1:L], list(L = nsite*nseason))
  RHS.Phi <- make.glm(phiform, factors = factors.size, cl, mixture, "logit", level = quote(phi))
  glm.expand.Phi <- nim_glm$process(LHS.Phi, RHS.Phi)
  lm.expand.Phi <- lmPredictor$process(glm.expand.Phi$LHS, glm.expand.Phi$RHS)  #See if We Need LHS here.

  LHS.Gamma <- substitute(z[1:L], list(L = nsite*nseason))
  RHS.Gamma <- make.glm(gammaform, factors = factors.size, cl, mixture, "logit", level = quote(gam))
  glm.expand.Gamma <- nim_glm$process(LHS.Gamma, RHS.Gamma)
  lm.expand.Gamma <- lmPredictor$process(glm.expand.Gamma$LHS, glm.expand.Gamma$RHS)  #See if We Need LHS here.

  #Mean of Latent States, muZ
  latent.mean <- latent.mean.loop(nsite, nseason)

  #Latent State z
  latent.phi <- embedLinesInForLoop(substitute(S ~ dbern(THETA), list(S = LHS2BUGSterm(quote(z), index.name), THETA = LHS2BUGSterm(quote(psi), index.name))), index.name, start = 1, finish = nsite) #latent state for first season
  latent.z <- latent.state.loop(nsite, nseason, index.name) #make latent state for season s + 1 to nseason

  #Observation Model
  LHS.P <- substitute(z[1:L], list(L = nsite*nseason*nsurvey))
  RHS.P <- make.glm(pform, factors = factors.size, cl, mixture, "logit", level = quote(p))
  glm.expand.P <- nim_glm$process(LHS.P, RHS.P)
  lm.expand.P <- lmPredictor$process(glm.expand.P$LHS, glm.expand.P$RHS)  #See if We Need LHS here.

  #Model for y
  Obs.Mod <- make.obs.loop(nsite, nseason, nsurvey)

  #Full BUGS Code Expansion
  full.expand <- embedLinesInCurlyBrackets(lines = list(lm.expand.Site$code, lm.expand.Phi$code, lm.expand.Gamma$code, latent.mean, latent.phi, latent.z,
                                                        lm.expand.P$code, Obs.Mod))

  #Make Factor Variables Numeric
  indx <- sapply(tidy.data.sort, is.factor)
  tidy.data.sort[indx] <- lapply(tidy.data.sort[indx], function(x) as.numeric(as.factor(x)))


  combine.inits <- list(lm.expand.Site$inits, lm.expand.Phi$inits, lm.expand.Gamma$inits, lm.expand.P$inits)
  start.values <- lapply(unlist(combine.inits), lmPred2init, initmcmc)

  #Inititalize z, latent occupancy
  z <- c(apply(y, c(1,3), max))
  start.values <- append(start.values, list(z = z))

  browser()

  nimMod.obj <- nimbleModel(code = full.expand, inits = start.values, constants = as.list(tidy.data.sort))

  if (burnin >= niter){stop('Burnin Needs to Be Less Than the Number of Iterations')}
  mod.output <- callNim(nimMod.obj, niter = niter, burn = burnin, chains = chains)


  if(is.list(mod.output)){mod.out = do.call(rbind, mod.output)} else{mod.out = mod.output}

  results <- list("Summary" = makeoutput(mod.output))   #Q: Do we want to backtransform parameters from log and logit scale?


  if(returncode==TRUE){results[["BUGScode"]] = full.expand}
  if(returnsamp==TRUE){results[["Samples"]] = mod.out}
  return(results)
}



#Still Doing Testing#


look <- nimble.dynamic.occ(psiformula = ~ 1, phiformula = ~ Season, gammaformula = ~ Season, pformula = ~ Season, y = z, site = NULL, site.season = NULL,
                  site.season.survey = NULL, priors = "Normal", dropbase = TRUE, burnin = 100, niter = 1000)








