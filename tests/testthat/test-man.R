library(testthat)

test_that('nimble.abund example works', {
    # Simulate Abundance Data
    R <- 75
    T <- 3
    X <- runif(n = R, min = -1, max = 1)
    A = c(rep('Red', 25), rep('Orange',20), rep('Blue', 30))
    sim.covariates <- data.frame("X" = X, "A" = A)
    lam <- exp(1 + 3*X)
    N <- rpois(n = R, lambda = lam)
    p <- plogis(-5*X)
    y <- matrix(NA, nrow = R, ncol = T)
    for (i in 1:T){y[,i] <- rbinom(n = R, size = N, prob =p)}

    #Fit Model
    model <- nimble.abund(
        siteformula = ~ 1 + (1|A) + X,
        obsformula = ~ X,
        y = y,
        sitevars = sim.covariates,
        initmcmc = 1,
        chains = 1)
})
