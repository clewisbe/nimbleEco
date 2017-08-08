#load data
library(nimble)
library(coda)
#model.input <- readRDS("sim.occ.rds")
#model.input2 <- readRDS("model.input2.rds")


#Season Lag Model with Logit Link
nimCode <- nimbleCode({
  {
    for (i in 1:100) {
      logit(psi[i]) <- s.intercept
    }
    {
    s.intercept ~ dnorm(0, sd = 1)
    }
  }
  {
  for (i in 1:1000) {
    logit(phi[i]) <- Season.effect.phi[Season[i]]
  }
  {
  {
    for (i in 1:9) {
      Season.effect.phi[i] ~ dnorm(0, sd = 1)
    }
    Season.effect.phi[10] <- 0
  }
  }
  }
  {
  for (i in 1:1000) {
    logit(gamma[i]) <- Season.effect.gam[Season[i]]
  }
  {
  {
    for (i in 1:9) {
      Season.effect.gam[i] ~ dnorm(0, sd = 1)
    }
    Season.effect.gam[10] <- 0
  }
  }
  }
  for (s in 1:9) {
    muZ[(100 * (s - 1) + 1):(100 * s)] <- z[(100 * (s - 1) + 1):(100 * s)] * phi[(100 * (s - 1) + 1):(100 * s)] +
      (1 - z[(100 * (s - 1) + 1):(100 * s)]) * gamma[(100 * (s - 1) + 1):(100 * s)]
  }
  for (i in 1:100) {
    z[i] ~ dbern(psi[i])
  }
  for (i in 1:900) {
    z[i + 100] ~ dbern(muZ[i])
  }
  {
  for (i in 1:10000) {
    logit(p[i]) <- Season.effect.p[Season[i]]
  }
  {
  {
    for (i in 1:10) {
      Season.effect.p[i] ~ dnorm(0, sd = 1)
    }
  }
  }
  }
  for (i in 1:10000) {
    muy[i] <- z[SiteSeason[i]] * p[i]
    y[i] ~ dbern(muy[i])
  }
})

#Season Lag Model without Logit Link
nimCode2 <- nimbleCode({
  {
    for (i in 1:100) {
      psi[i] <- s.intercept
    }
    {
    s.intercept ~ dunif(0, 1)
    }
  }
  {
  for (i in 1:1000) {
    phi[i] <- Season.effect.phi[Season[i]]
  }
  {
  {
    for (i in 1:9) {
      Season.effect.phi[i] ~ dunif(0,1)
    }
    Season.effect.phi[10] <- 0
  }
  }
  }
  {
  for (i in 1:1000) {
    gamma[i] <- Season.effect.gam[Season[i]]
  }
  {
  {
    for (i in 1:9) {
      Season.effect.gam[i] ~ dunif(0, 1)
    }
      Season.effect.gam[10] <- 0
  }
  }
  }
  for (s in 1:9) {
    muZ[(100 * (s - 1) + 1):(100 * s)] <- (z[(100 * (s - 1) + 1):(100 * s)] * phi[(100 * (s - 1) + 1):(100 * s)]) +
      (1 - z[(100 * (s - 1) + 1):(100 * s)]) * gamma[(100 * (s - 1) + 1):(100 * s)]
  }
  for (i in 1:100) {
    z[i] ~ dbern(psi[i])
  }
  for (i in 1:900) {
    z[i + 100] ~ dbern(muZ[i])
  }
  {
  for (i in 1:10000) {
    p[i] <- Season.effect.p[Season[i]]
  }
  {
  {
    for (i in 1:10) {
      Season.effect.p[i] ~ dunif(0, 1)
    }
  }
  }
  }
  for (i in 1:10000) {
    muy[i] <- z[SiteSeason[i]] * p[i]
    y[i] ~ dbern(muy[i])
  }
})


#Lauren's Model with Double Indexing
ss.ms.occ <- nimbleCode({
  ## Specify priors
  psi1 ~ dunif(0, 1)

  for(year in 1:(nyear-1)){
    phi[year] ~ dunif(0, 1)
    gamma[year] ~ dunif(0, 1)
    p[year] ~ dunif(0, 1)
  }
  p[nyear] ~ dunif(0, 1)

  ## Ecological submodel: Define state conditional on parameters
  for (site in 1:nsite){
    z[site,1] ~ dbern(psi1)
    for (year in 2:nyear){
      muZ[site,year]<- z[site,year-1]*phi[year-1] + (1-z[site,year-1])*gamma[year-1]
      z[site,year] ~ dbern(muZ[site,year])
    }
  }

  ## Observation model
  for (site in 1:nsite){
    for (rep in 1:nrep){
      for (year in 1:nyear){
        muy[site,rep,year] <- z[site,year]*p[year]
        y[site,rep,year] ~ dbern(muy[site,rep,year])
      }
    }
  }

})

#Data for Colin Model
# z <- model.input$data$y
# tidy.data <- tidy.dynam.occ(z)
#
# tidy.data.sort = arrange(tidy.data, Survey)  #Allows Looping Over Seasons w/o Skipping Rows with Multiple Surveys Per Season
#
# tidy.data.sort$SiteSeason <- rep(1:(100 * 10), times = 10, each = 1)  #for nested indexing over latent states
# tidy.data.sort$Season <- as.numeric(tidy.data.sort$Season)
# tidy.data.sort$Site <- as.numeric(tidy.data.sort$Site)


# # #Make Starting Values for Parameters and Latents
#  start.values <- list(s.intercept = 0.1, Season.effect.phi = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2),
#                       Season.effect.gam = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), Season.effect.p = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1))
#  zstart <- c(apply(model.input2$inits$z, c(1, 3), max))
#  start.values <- append(start.values, list(z = zstart))
#
# # model.input2 <- list(data = tidy.data.sort, inits = start.values)
# #Call Nimble and Pass the BUGS code and data
# callNim <- function(modobj, niter, burn, chains){
#   nimMod.C <- compileNimble(modobj)
#
#   #Set Up MCMC
#   config.Mod <- configureMCMC(nimMod.C, print = TRUE)
#   mod.MCMC <- buildMCMC(config.Mod)
#   C.mod.MCMC <- compileNimble(mod.MCMC)
#   samplesList <- runMCMC(C.mod.MCMC, niter = niter, nburnin = burn, nchains = chains, returnCodaMCMC = TRUE)
#   return(samplesList)
# }
#
# #Take the MCMC samples and return summary statistics
# makeoutput <- function(codaobj, sitevars, obsvars){
#   quantile <- summary(codaobj)$quantiles
#   meansd <- summary(codaobj)$statistics[,1:2]
#   rhat <- tryCatch(gelman.diag(codaobj)$psrf, error = function(e) NA)
#   Effsamp <- effectiveSize(codaobj)
#   output  <- as.data.frame(cbind(meansd, quantile, rhat, Effsamp))
#   names(output)[names(output) == 'Point est.'] <- 'Rhat'
#   output$`Upper C.I.`<-NULL
#   return(list(output))
# }
#
#
#
# nimMod.obj_lauren <- nimbleModel(code = ss.ms.occ, inits = model.input$inits, constants = model.input$constants, data = as.list(model.input$data))
# nimMod.obj_colin <- nimbleModel(code = nimCode, inits = model.input2$inits, constants = model.input2$data)
# mod.output_lauren <- callNim(nimMod.obj_lauren, niter = 10000, burn = 1000, chains = 1)
# mod.output_colin <- callNim(nimMod.obj_colin, niter = 2500, burn = 100, chains = 1)
#
#
# #Should Get Matching Outputs for Both
# laurenout <- list(Summary = makeoutput(mod.output_lauren))
# colinout2 <- list(Summary = makeoutput(mod.output_colin))
