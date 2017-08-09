require(testthat)
require(methods)
require(nimbleEcology)

#Make an automatic named list
namedList <- function(...){
  names <- as.list(substitute(list(...)))[-1L]
  result <- list(...)
  names(result) <- names
  result
}

## Tests Abundance Models
#Models with Fixed Effects
mod1 <- list(site = quote(~ 1), obs = quote(~ 1))
mod2 <- list(site = quote(~ 1), obs = quote(~ xm))
mod3 <- list(site = quote(~ 1 + xm), obs = quote(~ 1))
mod4 <- list(site = quote(~ 1 + A), obs = quote(~ 1 + weekday))
mod5 <- list(site = quote(~ A), obs = quote(~ weekday))
mod6 <- list(site = quote(~ 1 + ym + (ym|A)), obs = quote(~ 1 + A))
mod7 <- list(site = quote(~ 1 + ym + xm*ym), obs = quote(~ 1 + A))

#Models with Random Effects
mod8 <- list(site  = quote(~ (1|A)), obs = quote(~ (1|weekday)))
mod9 <- list(site  = quote(~ 1 + xm + (xm|A)), obs = quote(~ 1 + (1|weekday)))
mod10 <- list(site = quote(~ 1 + xm), obs = quote(~ 1 + (xm|weekday)))
mod11 <- list(site = quote(~ 1 + (1|A)), obs = quote( ~ 1 + xm))
mod12 <- list(site = quote(~ 1 + ym + (ym|A)), obs = quote(~ 1 + A))


#Put all Test Cases in List
all.abund.mod <- namedList(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12)



test_mcmc <- function(model, data, inits = 1, iter = 10, burnin = 5 , mixture = "Poisson", seed = 1, gold = NULL, compare = TRUE){
  if (data == "abundance"){
  dat <- readRDS("./inst/testdata/abundance.sim.rds")
  }
  set.seed(seed)
  runtest = eval(substitute(nimble.abund(siteformula = S, obsformula = O , y = dat$y, sitevars = dat$sitevars, obsvars = dat$obsvars,,mixture = mixture, priors = "Normal",
                                       dropbase = FALSE, niter = iter, burnin = burnin, initmcmc = inits, chains = 1, returncode = TRUE, returnsamp = TRUE),
                          list(S = model$site, O = model$obs)))

  out <- Abun_samples <- as.matrix(runtest$Samples)
  if(compare == TRUE){
    out = expect_equal(out, gold, info = paste("Gold File and MCMC posterior samples are not equal") )
  }
  else{
  return(out)
  }
}


## Only Run to Generate New Gold Files ##

#gold_list_abun <- list()
#for (i in 1:12){
#  gold_list_abun[[i]] <- test_mcmc(all.abund.mod[[i]], "abundance", compare = FALSE)
#}
#saveRDS(gold_list_abun, file = "./inst/testdata/abundance.gold.rds")






