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

compareFilesUsingDiff <- function(trialFile, correctFile, main = "") {
  if(main == "") main <- paste0(trialFile, ' and ', correctFile, ' do not match\n')
  diffOutput <- system2('diff', c(trialFile, correctFile), stdout = TRUE)
  test_that(paste0(main, paste0(diffOutput, collapse = '\n')),
            expect_true(length(diffOutput) == 0)
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


## Simulated Data  ##

# 20 Sites and 3 Visits
R <- 20
T <- 3
y <- array(dim=c(R,T))
N <- rpois(n=R, lambda =2)
for (j in 1:T){
  y[,j] <- rbinom(n=R, size=N, prob = .5)
}

y <- as.data.frame(y)

sitevars = list(z = rpois(20,2), xm = rnorm(20), ym = rgamma(20,3), A = c(rep('Red', 7), rep('Orange',3), rep('Blue', 10)))
sitevars = as.data.frame(sitevars)
x1 = rep('Mon', 20)
x2 = rep('Tues', 20)
x3 = rep('Wed', 20)
day = cbind(x1, x2, x3)

obsvars = list(m = matrix(data = rexp(60, rate = 10), nrow = 20, ncol = 3), weekday = day)

## Generate Model Cases to Test ##
mod1 <- c(quote(~1), quote(~1))
mod2 <- c(quote(~1), NULL)
mod3 <- c(quote(~1 + xm), quote(~1))
mod4 <- c(quote(~1 + (1|A)), quote(~1 + xm))
mod5 <- c(quote(~xm), quote(~1 + (1|A)))
mod6 <- c(quote(~1 + ym + (ym|A)), quote(~1 + A))


allmod <- namedList(mod1, mod2, mod3, mod4, mod5, mod6)


testcode <- function(S, O){
  runtest = eval(substitute(nimble.abund(siteformula = S, obsformula = O , y = y, sitevars = sitevars, obsvars = obsvars, mixture = "Poisson", priors = "Normal",
                      dropbase = FALSE, niter = 10, burnin = 1, initmcmc = 1, chains = 2, returncode = TRUE, returnsamp = TRUE),
                      list(S = S, O = O)))
  return(runtest$BUGScode)
}


#Generate Code Over All Models #

testoutput <- list()
for (i in 1:3){
  model.id <- (paste0(quote(mod), i))
  model.formula <- allmod[[model.id]]
  testoutput[[i]] <- testcode(model.formula[[1]], model.formula[[2]])
}


goldfile <- Combine.Code(testoutput)
goldfile2 <- goldfile


#Compare Output Files
compareFilesUsingDiff(goldfile, goldfile2)
compareFilesByLine(goldfile, goldfile2)


