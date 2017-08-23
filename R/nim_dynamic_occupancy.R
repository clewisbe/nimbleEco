#'Estimation of Dynamic (MultiSeason) Single-Species Occupancy Models
#'
#'\code{nimble.dynamic.occ} fits single species multi-season occupancy
#'models using Bayesian estimation.  The function returns model summary
#'statistics, posterior draws, and the BUGS code used to estimate the model.
#'
#'This is a function for fitting Multiseason Single Species Site Occupancy
#'Models using a Bayeisan framework.  The model structure is: \deqn{z_{1} ~
#'Bernoulli(\psi_{1})} \deqn{z_{i,k} ~ Bernoulli(\psi_{i,k})}
#'\deqn{z_{k+1}|z_{k} ~ Bernoulli(z_{k} * \phi_{i,k} + (1-z{k}) * \gamma_{i,k})}
#'\deqn{y_{i,j,k} | z_{i,k} ~ Bernoulli(z_{i,k} * p_{i,j,k})}
#'
#'Where i indexing site, j indexes survey, and k indexes the season. z is a
#'binary latent variable for the true observed state (occupied or not occupied).
#'The parameter \eqn{\psi_{1}} (first year occupancy) can be modeled as a linear
#'function of site covariates; \eqn{\phi_{i,k}} (probability of local survial),
#'and \eqn{\gamma_{i,k}} (colonization probabiliy) can be modeled as linear
#'functions of site and site.season covariates; and \eqn{p_{i,j,k}} (detection
#'probability) can be modeled as a linear function of site, site.season, or
#'site.season.survey covariates.  A logit link is used for all 3 models.    \cr
#'\cr Estimation is performed using NIMBLE.  In addition to fitting the model,
#'generated BUGS code for the full model, as well as posterior draws, are
#'returned.
#'
#'
#'@param psiformula  A linear model formula for first year occupancy:
#'  \eqn{logit(\psi_{1})}.  Formula must start with a tilde (~).  Add -1 to
#'  remove intercept. Random effects are specified with a pipe (see
#'  \href{https://www.rdocumentation.org/packages/lme4/versions/1.1-13/topics/lmer}{lmer}).
#'  If model contains an intercept random effects will be centered at 0.
#'  Interactions can be specified with * (e.g., x*y).
#'@param phiformula A linear model formula for survival probability:
#'  \eqn{logit(\phi_{i,k})}
#'@param gammaformula A linear model formula for the probability of
#'  colonization: \eqn{logit(\gamma_{i,k})}.
#'@param pformula A linear model formula for the probability of detection:
#'  \eqn{logit(p{i,j,k})}.
#'@param y A 3 dimensinoal array (site, survey, season) with binary values: 1
#'  for detection, 0 otherwise.
#'@param site A data frame or matrix.  Rows are sites and columns (with variable
#'  names) are independent variables.
#'@param site.season A named list (rows site, columns season) with each list
#'  element corresponding to an independent variable. Variables specified in the
#'  model must match the variable names in the list.
#'@param site.season.survey  A named list of arrays (site, survey, season) with
#'  each name corresponding to an independent variable. Variables specified in
#'  the model must match the variable names in the list.
#'@param priors One prior is specified for all model parameters in the
#'  siteformula and obsformula model statements.  The default prior is a
#'  Normal(mu = 0, sd = 1).  Other options are t(mu = 0, tau = 1000, df = 5),
#'  Uniform(0, 1000), and Gamma(shape = .001, rate = .001).
#'@param dropbase  Drops the first level of a factor variable.  Default is TRUE.
#'@param droplast Drops the last level of a factor variable.  Default is TRUE.
#'@param niter Number of iterations to run the MCMC.  Default is 10,000
#'  iterations.
#'@param burnin Number of burnin iterations for MCMC.  Default is 1000
#'  iterations.
#'@param initmcmc Vector of initial values for MCMC.  Length of vector must
#'  equal the number of chains.  The same initial values are used for all
#'  parameters. Default is c(1, 5, 10, 15)
#'@param chains Number of MCMC chains to run.  Default is 4.
#'@param returncode Option to return the BUGS code passed to NIMBLE.  Default is
#'  TRUE.
#'@param returnsamp Option to return all posterior draws.  Default is TRUE
#'@param ... Additional arguments can be passed to specify different parameters
#'  for the priors.  Options are location, scale, lb, ub, df, shape, and rate.
#'
#'
#'
#'
#'@return Output is a named list with the following elements: Summary, BUGScode,
#'  and Samples.  In the Summary statement, variables from the psiformula will
#'  start with "s.", variables in the phiformula will start with "phi.",
#'  variables from the gammaformula will start with "gam.", and variables from
#'  the pformula will start with "p.".  Also, note parameters are on the
#'  transformed scale (logit). In addition to quantiles, the effective sample
#'  size \link[coda]{effectiveSize} and Gelman Rubin diagnoistic
#'  \link[coda]{gelman.diag} are provided from the coda package.
#'
#'
#'
#'
#'@author Colin Lewis-Beck
#'
#' @examples
#' # Simualte Data R <- 100 #Number of Sites J <- 3 #Number of Surveys K <- 10
#' #Number of Seasons psi1 <- 0.4 #First year occupancy probability psi <-
#' rep(NA, K) muZ <- z <- array(dim = c(R, K)) y <- array(NA, dim = c(R, J, K))
#' psi[1] <- psi1 p <- runif(n = K, min = 0.20, max = 0.40) phi <- runif(n =
#' K-1, min = 0.60, max = 0.80) gamma <- runif(n = K-1, min = 0, max = 0.10)
#'
#' #Generate Occurance states z[,1] <- rbinom(R, 1, psi[1]) for (i in 1:R){ for
#' (k in 2:K){ muZ[k] <- z[i, k-1]*phi[k-1] + (1-z[i,k-1])*gamma[k-1] z[i,k] <-
#' rbinom(1,1,muZ[k]) } }
#'
#' #Generate detection/nondetection data for (i in 1:R){ for(k in 1:K){ prob <-
#' z[i,k]*p[k] for (j in 1:J){ y[i,j,k] <- rbinom(1, 1, prob) } } }
#'
#'
#' #Fit Model
#' model <- nimble.dynamic.occ(psiformula = ~ 1, phiformula = ~ -1 +
#' Season, gammaformula = ~ -1 + Season, pformula = ~ -1 + Season, y = y,
#' initmcmc = 1, chains = 1, dropbase = FALSE)
#'
#'
#'@export
nimble.dynamic.occ <- function(psiformula = NULL, phiformula = NULL, gammaformula = NULL, pformula = NULL, y = NULL, site = NULL, site.season = NULL,
                               site.season.survey = NULL, priors = c("Normal", "t", "Uniform", "Gamma"), dropbase = TRUE,
                               droplast = TRUE, niter = 10000, burnin = 1000, initmcmc = c(.5, 1, 5, 10), chains = 4, returncode = TRUE, returnsamp = TRUE, ...) {

  cl <- match.call.defaults()  #Gets All Arguments, Including Defaults
  mf <- match.call(expand.dots = FALSE)

  # Pull of Model Formulas
  psiform <- as.formula(mf$psiformula)
  phiform <- as.formula(mf$phiformula)
  gammaform <- as.formula(mf$gammaformula)
  pform <- as.formula(mf$pformula)

  nsite <- as.numeric(dim(y)[1])
  nsurvey <- as.numeric(dim(y)[2])
  nseason <- as.numeric(dim(y)[3])

  index.name <- quote(i)

  cl$priors <- match.arg(priors, c("Normal", "t", "Uniform", "Gamma"))
  prior.params = list(...)
  argnames <- names(prior.params)

  # Default Prior Set to ~ N(0,1)
  if(length(prior.params)==0){
    prior.params <- list(scale = 1)
  }

  # Make Input Data Tidy
  tidy.data <- tidy.dynam.occ(y, site, site.season, site.season.survey)
  tidy.data.sort <- tidy.data[order(tidy.data$Survey, decreasing = FALSE), ] #Allows Looping Over Seasons w/o Skipping Rows with Multiple Surveys Per Season
  #tidy.data.sort = dplyr::arrange(tidy.data, Survey)

  tidy.data.sort$SiteSeason <- rep(1:(nsite * nseason), times = nsurvey, each = 1)  #for nested indexing over latent states

  # Identify Factor Variables
  factors <- names(Filter(is.factor, tidy.data))
  factors.size <- lapply(factors, factor.bracket, tidy.data)


  # Ecological SubModels

  # Model for Site
  LHS.Site <- substitute(z[1:L], list(L = nsite))
  RHS.Site <- make.glm(psiform, factors = factors.size, cl, "Bernoulli", "logit", level = quote(site), prior.params = prior.params)
  glm.expand.Site <- nim_glm$process(LHS.Site, RHS.Site)
  lm.expand.Site <- lmPredictor$process(glm.expand.Site$LHS, glm.expand.Site$RHS)

  # Models for Phi and Gamma
  LHS.Phi <- substitute(z[1:L], list(L = nsite * nseason))
  RHS.Phi <- make.glm(phiform, factors = factors.size, cl, "Bernoulli", "logit", level = quote(phi), prior.params = prior.params )
  glm.expand.Phi <- nim_glm$process(LHS.Phi, RHS.Phi)
  lm.expand.Phi <- lmPredictor$process(glm.expand.Phi$LHS, glm.expand.Phi$RHS)  #See if We Need LHS here.

  LHS.Gamma <- substitute(z[1:L], list(L = nsite * nseason))
  RHS.Gamma <- make.glm(gammaform, factors = factors.size, cl,"Bernoulli" , "logit", level = quote(gam), prior.params = prior.params)
  glm.expand.Gamma <- nim_glm$process(LHS.Gamma, RHS.Gamma)
  lm.expand.Gamma <- lmPredictor$process(glm.expand.Gamma$LHS, glm.expand.Gamma$RHS)  #See if We Need LHS here.

  # Mean of Latent States, muZ
  latent.mean <- latent.mean.loop(nsite, nseason)

  # Latent State z
  latent.phi <- embedLinesInForLoop(substitute(S ~ dbern(THETA), list(S = LHS2BUGSterm(quote(z), index.name), THETA = LHS2BUGSterm(quote(psi), index.name))), index.name, start = 1,
                                    finish = nsite)  #latent state for first season
  latent.z <- latent.state.loop(nsite, nseason, index.name)  #make latent state for season s + 1 to nseason


  # Observation Model
  cl$droplast = FALSE #for observation model there is no markov structure
  LHS.P <- substitute(z[1:L], list(L = nsite * nseason * nsurvey))
  RHS.P <- make.glm(pform, factors = factors.size, cl, "Bernoulli", "logit", level = quote(p), prior.params = prior.params)
  glm.expand.P <- nim_glm$process(LHS.P, RHS.P)
  lm.expand.P <- lmPredictor$process(glm.expand.P$LHS, glm.expand.P$RHS)  #See if We Need LHS here.

  # Model for y
  Obs.Mod <- make.obs.loop(nsite, nseason, nsurvey)

  # Full BUGS Code Expansion
  full.expand <- embedLinesInCurlyBrackets(lines = list(lm.expand.Site$code, lm.expand.Phi$code, lm.expand.Gamma$code, latent.mean, latent.phi, latent.z, lm.expand.P$code, Obs.Mod))

  # Make Factor Variables Numeric
  indx <- sapply(tidy.data.sort, is.factor)
  tidy.data.sort[indx] <- lapply(tidy.data.sort[indx], function(x) as.numeric(as.factor(x)))


  #Combine Starting Values from all model formula
  combine.inits <- list(lm.expand.Site$inits, lm.expand.Phi$inits, lm.expand.Gamma$inits, lm.expand.P$inits)

  # Inititalize z, latent occupancy
  z <- c(apply(y, c(1, 3), max))

  #Make List of List for Initial Values
  inits.list <- list()
  for (i in 1:length(initmcmc)){
    start.values <- lapply(unlist(combine.inits), lmPred2init, initmcmc[i])
    inits.list[[i]] <- append(start.values, list(z = z))
  }

  # NIMBLE settings
  nimble::nimbleOptions(enableBUGSmodules = TRUE)
  nimble::nimbleOptions(verbose = FALSE)
  nimble::nimbleOptions(MCMCprogressBar = FALSE)

  nimMod.obj <- nimble::nimbleModel(code = full.expand, constants = as.list(tidy.data.sort))

  if (burnin >= niter) {
    stop("Burnin Needs to Be Less Than the Number of Iterations")
  }
  mod.output <- callNim(nimMod.obj, initial = inits.list, niter = niter, burn = burnin, chains = chains)


  if (is.list(mod.output)) {
    mod.out = do.call(rbind, mod.output)
  } else {
    mod.out = mod.output
  }

  results <- list(Summary = makeoutput(mod.output))

  # Q: Do we want to backtransform parameters from log and logit scale?


  if (returncode == TRUE) {
    results[["BUGScode"]] = full.expand
  }
  if (returnsamp == TRUE) {
    results[["Samples"]] = mod.out
  }
  return(results)
}


#look <- nimble.dynamic.occ( ~1, ~Season, ~ 1 + Season*elev, ~Season, y = detection, site = sitevars,site.season = seasonsitevars, site.season.survey = pvars,
             #                chains = 1, niter = 20, burnin = 10,  dropbase = FALSE)

