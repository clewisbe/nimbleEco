source(system.file(file.path('tests', 'test_utils.R'), package = 'nimbleEcology'))

context("Testing default MCMC Abundance Models")

RwarnLevel <- options('warn')$warn
options(warn = -1)
nimbleVerboseSetting <- nimbleOptions('verbose')

## If you do *not* want to write to results files
## comment out the sink() call below.  And consider setting verbose = FALSE
nimbleOptions(verbose = FALSE)
tempFileName <- 'mcmcTestAbunLog.Rout'

sink(tempFileName)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

## Tests Abundance Models
abund.gold <- readRDS(system.file("testdata", 'abundance.gold.rds',package="nimbleEcology"))

test_that("mcmc compare abund",{
          for (i in 8:8){
          test_mcmc(all.abund.mod[[i]], "abundance", seed = 1, gold = abund.gold[[i]], compare = TRUE, mod = i)
            }
        })

sink(NULL)

## Restore Defaults
options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)

