source(system.file(file.path('tests', 'test_utils.R'), package = 'nimbleEcology'))

context("Testing default MCMC Single Season Dynamic Occupancy Models")

RwarnLevel <- options('warn')$warn
options(warn = -1)
nimbleVerboseSetting <- nimbleOptions('verbose')

## If you do *not* want to write to results files
## comment out the sink() call below.  And consider setting verbose = FALSE
nimbleOptions(verbose = FALSE)
tempFileName <- 'mcmcTestDynOccuLog.Rout'

sink(tempFileName)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

## Tests Single Season Dynamic Occupancy Models
dyn.single.gold <- readRDS(system.file("testdata", 'dynamic.single.gold.rds', package="nimbleEcology"))

test_that("mcmc compare single dyn occu",{
  for (i in 1:1){
    test_mcmc(all.dynam.occup.mod[[i]], "dynoccup", seed = 1, gold = dyn.single.gold[[i]], compare = TRUE, mod = i)
  }
})

sink(NULL)

## Restore Defaults
options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)

