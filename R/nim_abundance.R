#'Estimation of Abundance from Counts Using Binomial Mixture Models.
#'
#'\code{nimble.abund} is used to fit Binomial Mixture Models using Bayesian
#'estimation.  The function returns model summary statistics, posterior draws,
#'and the BUGS code used to estimate the model.
#'
#'This is a function for fitting Binomial Mixture Models using a Bayesian
#'framework.  The model structure is: \deqn{N_{i} ~ Poisson(\lambda_{i})}
#'\deqn{y_{i,j} | N_{i} ~ Binomial(N_{i}, p_{i,j})} Where i is the number of
#'sites and j the number of temporal replicates.  N is a latent variable for
#'population size. Both \eqn{\lambda} and p can be modeled as a linear function
#'of covariates via a linked GLM.  The link for \eqn{\lambda} is the log link.
#'The link for p is the logit link.  Estimation is performed using NIMBLE.  In
#'addition to fitting the model, generated BUGS code for the full model, as well
#'as posterior draws, are returned.
#'
#'
#'@param siteformula  A linear model formula for the latent state:
#'  \eqn{log(\lambda_{i})}.  Formula must start with a tilde (~).  Add -1 to
#'  remove intercept. Random effects are specified with a pipe (see
#'  \href{https://www.rdocumentation.org/packages/lme4/versions/1.1-13/topics/lmer}{lmer}).
#'  If model contains an intercept, random effects will be centered at 0.
#'  Interactions are specified with * (e.g., x*y).
#'@param obsformula A linear model formula for the observation process:
#'  \eqn{logit(p_{i,j})}.
#'@param y The observed counts at site i during temporal replicate j.  Must be
#'  provided as either a matrix or dataframe. Rows are sites and columns are
#'  replicates.
#'@param sitevars A dataframe of site-specific covariates.  Must be a dataframe
#'  with named columns corresponding to variables.  Number of rows should equal
#'  the number of sites.
#'@param obsvars A named list with site and survey specific covariates.  Each
#'  list element is a matrix with rows corresponding to site and columns to
#'  replicate. Variables specified in the model formula must match the variable names in
#'  the list.
#'@param mixture There are three mixture options for the latent states: Poisson
#'  (default), Zero-Inflated Poisson (ZIP), and Negative Binomial (NB).  The ZIP
#'  mixture adds an extra hierarchial layer to the model, which is a Bernoulli
#'  random variable that divides sites between positive counts and sites with 0
#'  counts.  The NB link adds an extra term (rho) to model overdispersion.
#'@param priors One prior is specified for all model parameters in the
#'  siteformula and obsformula model statements.  The default prior is a
#'  Normal(mu = 0, sd = 1000).  Other options are t(mu = 0, tau = 1000, df = 5),
#'  Uniform(0, 1000), and Gamma(shape = .001, rate = .001).
#'@param dropbase  Drops the first level of a factor variable.  Default is TRUE.
#'@param droplast Drops the last level of a factor variable.  Default is FALSE.
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
#'@return Output is a named list with the following elements: Summary, BUGScode,
#'  and Samples.  In the Summary statement, variables from the siteformula will
#'  start with "s.", and variables in the obsformula will start with "o."  Note:
#'  parameters are on the transformed scale (log for siteformula variables;
#'  logit for obsformula variables). In addition to quantiles, the effective
#'  sample size \link[coda]{effectiveSize} and Gelman Rubin diagnoistic
#'  \link[coda]{gelman.diag} are provided from the coda  package. \cr \cr If the
#'  ZIP mixture is selected theta is returned corresponding to the Bernoulli
#'  random variable modeling if a site is suitable for a positive abundance
#'  count. \cr \cr If the NB mixture is selected log(alpha) is returned:
#'  the log of the dispersion parameter, alpha, in the NB distribution.
#'
#'
#'
#'@author Colin Lewis-Beck
#'
#' @examples
#' # Simulate Abundance Data
#' R <- 75
#' T <- 3
#' X <- runif(n = R, min = -1, max = 1)
#' A = c(rep('Red', 25), rep('Orange',20), rep('Blue', 30))
#' sim.covariates <- data.frame("X" = X, "A" = A)
#' lam <- exp(1 + 3*X)
#' N <- rpois(n = R, lambda = lam)
#' p <- plogis(-5*X)
#' y <- matrix(NA, nrow = R, ncol = T)
#' for (i in 1:T){y[,i] <- rbinom(n = R, size = N, prob =p)}
#'

#'#Fit Model
#'model <- nimble.abund(siteformula = ~ 1 + (1|A) + X, obsformula = ~
#'X, y = y, sitevars = sim.covariates, initmcmc = 1, chains = 1)
#'
#'
#'@export
nimble.abund <- function(siteformula = NULL, obsformula = NULL, y = NULL, sitevars = NULL, obsvars = NULL, mixture = c("Poisson", "ZIP", "NB"), priors = c("Normal", "t", "Uniform", "Gamma"),
                         dropbase = TRUE, droplast = FALSE, niter = 10000, burnin = 1000, initmcmc = c(.5, 1, 5, 10), chains = 4, returncode = TRUE,
                         returnsamp = TRUE, ...) {

  cl <- match.call.defaults()
  mf <- match.call(expand.dots = FALSE)

  obsform <- as.formula(mf$obsformula)
  stateform <- as.formula(mf$siteformula)

  if (length(obsform) == 0 | length(stateform) == 0) {
    stop("A State and Observation Formula Must be Provided")
  }

  if (length(initmcmc) != chains) {
    stop("Number of initial values must equal number of MCMC chains")
  }

  prior.params = list(...)
  argnames <- names(prior.params)

  mixture <- match.arg(mixture, c("Poisson", "ZIP", "NB"))
  cl$priors <- match.arg(priors, c("Normal", "t", "Uniform", "Gamma"))


  df <- maketidy(y, sitevars, obsvars)  #Munge 3 Data Objects into One Tidy DataFrame

  # Length of Tidy df
  len <- as.numeric(dim(df)[1])

  # No. Sites
  S <- as.numeric(dim(y)[1])

  # Convert Character Variables to Factors
  for (i in 1:dim(df)[2]) {
    df[, i] <- makeFac(df[, i], char.only = TRUE)
  }

  factors <- names(Filter(is.factor, df))
  factors.size <- lapply(factors, factor.bracket, df)


  # Make Site Level Formula
  LHS.Site <- substitute(N[1:L], list(L = S))
  RHS.Site <- make.glm(mf$siteformula, factors.size, cl = cl, mixture, "log", level = quote(site), prior.params = prior.params)

  # Make Obs Level Formula
  LHS.Obs <- substitute(y[1:S], list(S = len))
  RHS.Obs <- make.glm(mf$obsformula, factors.size, cl = cl, "Binomial", "logit", N = quote(N), Site = quote(Site), level = quote(obs),
                      prior.params = prior.params)

  # Use BUGS Modeules to Expand Code
  glm.expand.Site <- nim_glm$process(LHS.Site, RHS.Site)
  glm.expand.Obs <- nim_glm$process(LHS.Obs, RHS.Obs)

  lm.expand.Site <- lmPredictor$process(glm.expand.Site$LHS, glm.expand.Site$RHS)
  lm.expand.Obs <- lmPredictor$process(glm.expand.Obs$LHS, glm.expand.Obs$RHS)


  # Full BUGS Code Expansion
  full.expand <- embedLinesInCurlyBrackets(lines = list(glm.expand.Site$prob.mod, lm.expand.Site$code, glm.expand.Obs$prob.mod, lm.expand.Obs$code))
  # Get Data and Model Ready for nimbleModel and MCMC

  # Make Factor Variables Numeric for NIMBLE
  indx <- sapply(df, is.factor)
  df[indx] <- lapply(df[indx], function(x) as.numeric(as.factor(x)))

  # Inititalize N, the latent counts
  N <- apply(y, 1, max)
  N[is.na(N)] <- 0
  N = N + 1

  #Make List of List for Initial Values
  #For Now Users Cannot Change Initial Values for Theta or LogAlpha
  inits.list <- list()
  for (i in 1:length(initmcmc)){
    start.values <- lapply(append(lm.expand.Site$inits, lm.expand.Obs$inits), lmPred2init, initmcmc[i])
    if (mixture == "NB") {
      start.values <- append(start.values, list(logalpha = 1))
    }

    if (mixture == "ZIP") {
      start.values <- append(start.values, list(theta = 0.5))
    }

    inits.list[[i]] <- append(start.values, list(N = N))
  }

  # NIMBLE settings
  nimble::nimbleOptions(enableBUGSmodules = TRUE)
  nimble::nimbleOptions(verbose = FALSE)
  nimble::nimbleOptions(MCMCprogressBar = FALSE)

  nimMod.obj <- nimble::nimbleModel(code = full.expand, constants = as.list(df), data = list(y = df$Count))

    if (burnin >= niter) {
    stop("Burnin Needs to Be Less Than the Number of Iterations")
  }

  mod.output <- callNim(nimMod.obj, initial = inits.list, niter = niter, burn = burnin, chains = chains)

  if (is.list(mod.output)) {
    mod.out = do.call(rbind, mod.output)
  } else {
    mod.out = mod.output
  }


  results <- list(Summary = makeoutput(mod.output, names(lm.expand.Site$inits), names(lm.expand.Obs$inits)))

  #Q: Do we want to backtransform parameters from log and logit scale?

  if (returncode == TRUE) {
    results[["BUGScode"]] = full.expand
  }
  if (returnsamp == TRUE) {
    results[["Samples"]] = mod.out
  }
  return(results)
}
