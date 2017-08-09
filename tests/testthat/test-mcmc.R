source(system.file(file.path('tests', 'test_utils.R'), package = 'nimbleEcology'))

context("Testing default MCMC Abundance Models")

RwarnLevel <- options('warn')$warn
options(warn = -1)
nimbleVerboseSetting <- nimbleOptions('verbose')

## If you do *not* want to write to results files
##    comment out the sink() call below.  And consider setting verbose = FALSE
## To record a new gold file, nimbleOptions('generateGoldFileForMCMCtesting') should contain the path to the directory where you want to put it
## e.g. nimbleOptions(generateGoldFileForMCMCtesting = getwd())
## Comparison to the gold file won't work until it is installed with the package.
nimbleOptions(verbose = FALSE)
tempFileName <- 'mcmcTestLog.Rout'

sink(tempFileName)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

## Tests Abundance Models
abund.gold <- readRDS(system.file("testdata", 'abundance.gold.rds',package="nimbleEcology"))

test_that("mcmc compare abund",{
          for (i in 1:1){
          test_mcmc(all.abund.mod[[i]], "abundance", seed = 1, gold = abund.gold[[i]], compare = TRUE)
            }
        })

sink(NULL)

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)

