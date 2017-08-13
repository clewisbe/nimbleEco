#'Estimation of Single-Season Occupancy Models
#'
#'\code{nimble.occ} is used to fit single species single-season occupancy models
#'using Bayesian estimation.  The function returns model summary statistics,
#'posterior draws, and the BUGS code used to estimate the model.
#'
#'This is a function for fitting Single Season-Occupany Models using a Bayeisan
#'framework.  The model structure is: \deqn{z_{i} ~
#'Bernoulli(\psi_{i})} \deqn{y_{i,j} | z_{i} ~ Bernoulli(z_{i} * p_{i,j})} Where
#'i is the number of sites and j the number of temporal replicates. z is a
#'binary latent variable for the true observed state (occupied or not occupied).
#'Both \eqn{\psi} and p can be modeled as a linear function of covariates via a
#'linked GLM.  The link for \eqn{\psi} and p is the logit link. Estimation is
#'performed using NIMBLE.  In addition to fitting the model, generated BUGS code
#'for the full model, as well as posterior draws, are returned.
#'
#'@param siteformula  A linear model formula for the latent state:
#'  \eqn{log(\lambda_{i})}.  Formula must start with a tilde (~).  Add -1 to
#'  remove intercept. Random effects are specified with a pipe (see
#'  \link[lme4]{glmer}).  If model contains an intercept random effects will be
#'  centered at 0.  Interactions can be specified with * (e.g., x*y).
#'@param obsformula A linear model formula for the observation process:
#'  \eqn{logit(p_{i,j})}.
#'@param y Detection and non detection data.  Must be provided as either a
#'  matrix or dataframe. Rows are sites and columns are replicates.  1 is
#'  detection and 0 is non detection.
#'@param sitevars A dataframe of site-specific covariates.  Must be a dataframe
#'  with named columns corresponding to variables.  Number of rows should equal
#'  the number of sites.
#'@param obsvars A named list with site and survey specific covariates.  Each
#'  list element is a matrix with rows corresponding to site and columns to
#'  replicate. Variables specified in the model must match the variable names in
#'  the list.
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
#'@return Output is a named list with the following elements: Summary,
#'  BUGScode, and Samples.  In the Summary statement, variables from the
#'  siteformula will start with "s.", and variables in the obsformula will start
#'  with "o."  Note: parameters are on the transformed scale (log for
#'  siteformula variables; logit for obsformula variables). In addition to
#'  quantiles, the effective sample size and Gelman Rubin diagnoistic are
#'  provided from the \link[coda]{coda} package.
#'
#'
#'@author Colin Lewis-Beck
#'
#' @examples
#' # Simualte Occupancy Data
#' R <- 200
#' T <- 3
#' X <- sort(runif(n = R, min = -1, max = 1))
#' Y <- sort(rgamma(n = R, shape = 1, rate = 1))
#' sim.covariates <- data.frame("X" = X, "Y" = Y)
#' psi <- plogis(-1 + 2*X)
#' p <- plogis(1 - 3*X + 4*Y + 4*X*Y)
#' z <- rbinom(n = R, size = 1, prob = psi)
#' y <- matrix(NA, nrow = R, ncol = T)
#' for (i in 1:T){y[,i] <- rbinom(n = R, size = 1, prob = z*p)}
#'
#'
#' #Fit Model
#' model <- nimble.occ(siteformula = ~ 1 + X, obsformula = ~ 1 + X + Y + X*Y, y = y, sitevars = sim.covariates, initmcmc = c(1, 1.5, 3, 5), chains = 4)
#'
#'
#'@export
nimble.occ <- function(siteformula = NULL, obsformula = NULL, y = NULL, sitevars = NULL, obsvars = NULL, priors = c("Normal", "t", "Uniform", "Gamma"),
                         dropbase = TRUE, droplast = FALSE, niter = 10000, burnin = 1000, initmcmc = c(.5, 1, 5, 10), chains = 4, returncode = TRUE, returnsamp = TRUE, ...) {

  cl <- match.call.defaults()
  mf <- match.call(expand.dots = FALSE)

  prior.params = list(...)
  argnames <- names(prior.params)

  cl$priors <- match.arg(priors, c("Normal", "t", "Uniform", "Gamma"))

  obsform <- as.formula(mf$obsformula)
  stateform <- as.formula(mf$siteformula)

  if (length(obsform) == 0 | length(stateform) == 0) {
    stop("A State and Observation Formula Must be Provided")
  }

  df <- maketidy(y, sitevars, obsvars)  #turn 3 data pieces into one tiny data frame


  # Length of Long DF
  L <- as.numeric(dim(df)[1])

  # No. Sites
  S <- as.numeric(dim(y)[1])

  # Convert Character Variables to Factors
  for (i in 1:dim(df)[2]) {
    df[, i] <- makeFac(df[, i], char.only = TRUE)
  }

  # Make Site Level Formula #
  LHS.Site <- substitute(z[1:L], list(L = S))
  factors <- names(Filter(is.factor, df))
  factors.size <- lapply(factors, factor.bracket, df)
  RHS.Site <- make.glm(mf$siteformula, factors.size, cl, "Bernoulli", "logit", level = quote(site), prior.params = prior.params)
  LHS.Obs <- substitute(y[1:S], list(S = L))
  RHS.Obs <- make.glm(mf$obsformula, factors.size, cl, "Bernoulli", "logit", N = quote(N), Site = quote(Site), level = quote(obs),
                      prior.params = prior.params)


  glm.expand.Site <- nim_glm$process(LHS.Site, RHS.Site)
  glm.expand.Obs <- nim_glm$process(LHS.Obs, RHS.Obs)

  lm.expand.Site <- lmPredictor$process(glm.expand.Site$LHS, glm.expand.Site$RHS)
  lm.expand.Obs <- lmPredictor$process(glm.expand.Obs$LHS, glm.expand.Obs$RHS)

  # Full BUGS Code Expansion
  full.expand <- embedLinesInCurlyBrackets(lines = list(glm.expand.Site$prob.mod, lm.expand.Site$code, glm.expand.Obs$prob.mod, lm.expand.Obs$code))

  # Get Data and Model Readyfor nimbleModel and MCMC

  # Make Factor Variables Numeric
  indx <- sapply(df, is.factor)
  df[indx] <- lapply(df[indx], function(x) as.numeric(as.factor(x)))

  # Inititalize z, the latent observed/non observed
  z <- apply(y, 1, max)
  z[is.na(z)] <- 0

  #Make List of List for Initial Values
  inits.list <- list()
  for (i in 1:length(initmcmc)){
    start.values <- lapply(append(lm.expand.Site$inits, lm.expand.Obs$inits), lmPred2init, initmcmc[i])
    inits.list[[i]] <- append(start.values, list(z = z))
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
