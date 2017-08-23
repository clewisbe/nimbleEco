#source(system.file(file.path('tests', 'test_utils.R'), package = 'nimbleEcology'))
source("test_utils.R")

context("Testing default MCMC Occupancy Models")

RwarnLevel <- options('warn')$warn
options(warn = -1)
nimbleVerboseSetting <- nimbleOptions('verbose')

## If you do *not* want to write to results files
## comment out the sink() call below.  And consider setting verbose = FALSE
nimbleOptions(verbose = FALSE)
tempFileName <- 'mcmcTestOccuLog.Rout'

sink(tempFileName)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

## Tests Occupancy Models
#occu.gold <- readRDS(system.file("testdata", 'occupancy.gold.rds', package="nimbleEcology"))
occu.gold <- readRDS("../testdata/occupancy.gold.rds")

test_that("mcmc compare occup",{
          for (i in 1:12){
          test_mcmc(all.abund.mod[[i]], "occupancy", seed = 1, gold = occu.gold[[i]], compare = TRUE, mod = i)
            }
        })

sink(NULL)

## Restore Defaults
options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)

