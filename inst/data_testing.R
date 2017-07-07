library(devtools)
install_github("nimble-dev/nimble", ref = "devel", subdir = "packages/nimble")
library(nimble)
nimbleOptions(enableBUGSmodules = TRUE)

##Test Data
sillyData <- list( B = c(rep("colin",10)), z = rpois(10,2), q=rnorm(10), xm = rnorm(10), ym=rgamma(10,3), A = c(rep('Red',7), rep('Orange',1), rep('Blue', 2)))


#Works for Fixed Effects, Still a Bug for Random
#Returns Summary Stats and Nimble Code
test <- nimble.lm(mod = z ~ xm + A, dat = sillyData, priors="Normal", dropbase = FALSE, returncode = TRUE)


#test2 <- nim_lm$process(test$LHS, test$RHS)
#test3 <- lmPred$process(test2$LHS, test2$RHS)

