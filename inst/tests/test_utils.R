require(testthat)
require(nimble)
require(nimbleEcology)

#Make an automatic named list
namedList <- function(...){
  names <- as.list(substitute(list(...)))[-1L]
  result <- list(...)
  names(result) <- names
  result
}

#Models for Abundance and Occupancy Models
#Specifications with Fixed Effects
mod1 <- list(site = quote(~ 1), obs = quote(~ 1))
mod2 <- list(site = quote(~ 1), obs = quote(~ -1 + xm))
mod3 <- list(site = quote(~ 1 + xm), obs = quote(~ 1))
mod4 <- list(site = quote(~ 1 + A), obs = quote(~ 1 + weekday))
mod5 <- list(site = quote(~ -1 + A), obs = quote(~ -1 + weekday))
mod6 <- list(site = quote(~ 1 + ym + (ym|A)), obs = quote(~ 1 + A))
mod7 <- list(site = quote(~ 1 + ym + xm*ym), obs = quote(~ 1 + A))

#Specifications with Random Effects
mod8 <- list(site  = quote(~ -1 + (1|A)), obs = quote(~ -1 + (1|weekday)))
mod9 <- list(site  = quote(~ 1 + xm + (xm|A)), obs = quote(~ 1 + (1|weekday)))
mod10 <- list(site = quote(~ 1 + xm), obs = quote(~ 1 + (xm|weekday)))
mod11 <- list(site = quote(~ 1 + (1|A)), obs = quote( ~ 1 + xm))
mod12 <- list(site = quote(~ 1 + ym + (ym|A)), obs = quote(~ 1 + A))

#Models for Single Season Dynamic Occupancy (Need to Add Additional Checks)
mod.dyn1 <- list(psi = quote(~ 1), phi = quote(~ -1 + Season), gam = quote(~ -1 + Season), p = quote(~ -1 + Season))


#Put all Test Cases in List
all.abund.mod <- namedList(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12)
all.occup.mod <- namedList(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12)
all.dynam.occup.mod <-namedList(mod.dyn1)



test_mcmc <- function(model, data, inits = 1, iter = 10, burnin = 5 , mixture = "Poisson", seed = 1, gold = NULL, compare = TRUE, mod = NULL){
  if (data == "abundance"){
    dat <- readRDS(system.file("testdata", 'abundance.sim.rds', package="nimbleEcology"))
    M <- quote(nimble.abund)
  }
  if (data == "occupancy"){
    dat <- readRDS(system.file("testdata", 'occupancy.sim.rds', package="nimbleEcology"))
    M <- quote(nimble.occ)
  }
  if (data == "dynoccup"){
    dat <- readRDS(system.file("testdata", 'dynam.sim.rds', package="nimbleEcology"))
    M <- quote(nimble.dynamic.occ)
  }


  set.seed(seed)
  if(data == "abundance" | data == "occupancy"){
    runtest <- eval(substitute(MODEL(siteformula = S, obsformula = O , y = dat$y, sitevars = dat$sitevars, obsvars = dat$obsvars, mixture = mixture, priors = "Normal",
                                       dropbase = FALSE, niter = iter, burnin = burnin, initmcmc = inits, chains = 1, returncode = TRUE, returnsamp = TRUE),
                          list(S = model$site, O = model$obs, MODEL = M)))
  }

  if(data == "dynoccup"){
    runtest <- eval(substitute(MODEL(psiformula = S, phiformula = O , gammaformula = G, pformula = P, y = dat, priors = "Normal",
                                     dropbase = FALSE, niter = iter, burnin = burnin, initmcmc = inits, chains = 1, returncode = TRUE, returnsamp = TRUE),
                               list(S = model$psi, O = model$phi, G = model$gam, P = model$p, MODEL = M)))
  }

  out <- as.matrix(runtest$Samples)
  if(compare == TRUE){
    out = expect_equal(out, gold, info = paste("Gold File and MCMC posterior samples are not equal for model", mod))
  }
  else{
    return(out)
  }
}


## Only Uncomment to Generate New Gold Files ##

#gold_list_abun <- list()
#gold_list_occu <- list()
#gold_list_dynam <- list()
#for (i in 1:12){
#  gold_list_abun[[i]] <- test_mcmc(all.abund.mod[[i]], "abundance", compare = FALSE)
#  gold_list_occu[[i]] <- test_mcmc(all.occup.mod[[i]], "occupancy", compare = FALSE)
#   gold_list_dynam[[i]] <- test_mcmc(all.dynam.occup.mod[[i]], "dynoccup", compare = FALSE)
#  }
#saveRDS(gold_list_abun, file = "./inst/testdata/abundance.gold.rds")
#saveRDS(gold_list_occu, file = "./inst/testdata/occupancy.gold.rds")
#saveRDS(gold_list_dynam, file = "./inst/testdata/dynamic.single.gold.rds")






