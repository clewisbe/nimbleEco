# Test Code Generation with Gold Files #
library(testthat)

## Code Testing Functions ##
Combine.Code <- function(...) {
  ll <- list(...)
  ll <- lapply(ll, function(X) as.list(X)[-1])
  ll <- do.call("c", ll)
  as.call(c(as.symbol("{"), ll))
}

#Make an automatic named list
namedList <- function(...){
  names <- as.list(substitute(list(...)))[-1L]
  result <- list(...)
  names(result) <- names
  result
}


#Status 1 Indicates the Test Failed
compareFilesUsingDiff <- function(trialFile, correctFile, main = "") {
  if(main == "") main <- paste0(trialFile, ' and ', correctFile, ' do not match\n')

  if(.Platform$OS.type == 'windows'){
    diffOutput <- system2('fc', c(trialFile, correctFile), stdout = TRUE)}

  else{
    diffOutput <- system2('diff', c(trialFile, correctFile), stdout = TRUE)}

  test_that(paste0(main, paste0(diffOutput, collapse = '\n')),
            expect_true(length(diffOutput) < 4)  #testing against length 0 always fails
  )
  invisible(NULL)
}



compareFilesByLine <- function(trialResults, correctResults, main = "") {
  test_that(paste0(main, ': same number of output lines'),
            expect_equal(length(trialResults), length(correctResults)))

  linesToTest <- min(length(trialResults), length(correctResults))
  mapply(function(lineno, trialLine, correctLine) {
    test_that(paste0(main, ": output line #", lineno),
              expect_identical(trialLine, correctLine))
  }, 1:linesToTest, trialResults, correctResults)
  invisible(NULL)
}


## Simulated Data for Abundance Models  ##

#Generate Dependent Varialbe for Abundance or Occupancy Models
# 20 Sites and 3 Visits
sim.y <- function(S = 20, V = 3, prob = .5, lam = 2, model){
  y <- array(dim = c(S, V))
  if (model == "abundance"){
  N <- rpois(n = S, lambda = lam)
  }
  else{
    N <- 1
  }
  for (j in 1:V){
    y[,j] <- rbinom(n = S, size = N, prob = .5)
  }
  return(as.data.frame(y))
}


y.a <- sim.y(model = "abundance")
y.o <- sim.y(model = "occupancy")

sitevars = as.data.frame(list(zm = rpois(20,2), xm = rnorm(20), ym = rgamma(20,3), A = c(rep('Red', 7), rep('Orange',3), rep('Blue', 10))))
x1 = rep('Mon', 20)
x2 = rep('Tues', 20)
x3 = rep('Wed', 20)
day = cbind(x1, x2, x3)

obsvars = list(growth = matrix(data = rexp(60, rate = 10), nrow = 20, ncol = 3), weekday = day)

abundance.sim <- list(y = y.a , sitevars = sitevars, obsvars =  obsvars)
occupancy.sim <- list(y = y.o, sitevars = sitevars, obsvars = obsvars)
saveRDS(abundance.sim, "abundance.sim.rds")
saveRDS(occupancy.sim, "occupancy.sim.rds")

## Generate Model Cases to Test ##
#First column is Site Model; 2nd is for Obs Model
mod1 <- c(quote(~ 1), quote(~ 1))
mod2 <- c(quote(~ 1), quote(~ xm))
mod3 <- c(quote(~ 1 + xm), quote(~ 1))
mod4 <- c(quote(~ 1 + A), quote(~ 1 + weekday))
mod5 <- c(quote(~ A), quote(~ weekday))
mod6 <- c(quote(~ 1 + ym + (ym|A)), quote(~ 1 + A))
mod7 <- c(quote(~ 1 + ym + xm*ym), quote(~ 1 + A))

#Models with Random Effects
mod8 <- c(quote(~ (1|A)), quote(~ (1|weekday)))
mod9 <- c(quote(~ 1 + xm + (xm|A)), quote(~ 1 + (1|weekday)))
mod10 <- c(quote(~ 1 + xm), quote(~1 + (m|weekday)))
mod11 <- c(quote(~ 1 + (1|A)), quote(~1 + xm))
mod12 <- c(quote(~ 1 + ym + (ym|A)), quote(~1 + A))


#Put all Test Cases in List
allmod <- namedList(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12)

testcode <- function(S, O, mixture){
  runtest = eval(substitute(nimble.abund(siteformula = S, obsformula = O , y = y, sitevars = sitevars, obsvars = obsvars, mixture = mixture, priors = "Normal",
                      dropbase = FALSE, niter = 10, burnin = 1, initmcmc = 1, chains = 2, returncode = TRUE, returnsamp = TRUE),
                      list(S = S, O = O)))
  return(runtest$BUGScode)
}


#Generate Code Over All Models #

testoutput <- list()
for (i in 1:12){
  model.id <- (paste0(quote(mod), i))
  model.formula <- allmod[[model.id]]
  testoutput[[i]] <- testcode(model.formula[[1]], model.formula[[2]], mixture ="ZIP")
}


goldfile <- Combine.Code(testoutput)
dput(goldfile, file="NMixPoiGold.txt")  #save .txt file for Binomial/Poisson N Mixture 7/11/17



#Compare Output Files
compareFilesUsingDiff(goldfile, goldfile2)
compareFilesByLine(goldfile, goldfile2)

## Simulate Data for Dynamic Occupancy Models ##
#Toy Data for Testing

#data <- readRDS("dyn.occu.test.dat.rds")
#detection <- data$data$y
#detection <- detection[,,-1] #Site, Survey, Season


#Simulated Covariates for Site and Season
myMat1 <-matrix(runif(9*100), ncol=9)
myMat2 <-matrix(rnorm(9*100), ncol=9)
seasonsitevars <- list(elev = myMat1, height = myMat2)

#Simulated Season 1 Site Vars
sitevars <- as.data.frame(matrix(rexp(100)))
colnames(sitevars) <- "length"

#ObsSiteLocationVars
pvars <- list(heightcm = detection, time = detection)




