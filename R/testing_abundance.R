#Abundance
library(unmarked)
#library(reshape2)
#library(tidyr)
#library(dplyr)
#library(nimble)

#Mallard Data from Unmarked Packcage
#Mallard Data from Unmarked Packcage
y <- as.data.frame(mallard.y)
sitevars <- as.data.frame(mallard.site)
obsvars <- mallard.obs

#Example with Constant Parameters
R <- 200
T <- 3
y <- array(dim=c(R,T))
N <- rpois(n=R, lambda =2)
for (j in 1:T){
  y[,j] <- rbinom(n=R, size=N, prob = .5)
}



test <- nimble.abund(siteformula = ~ 1, ~ 1, y = y, sitevars = NULL, obsvars = NULL, mixture = "Poisson", priors = "Normal",
                     dropbase = FALSE, niter = 10000, burnin = 1, initmcmc = 1, chains = 4, returncode = TRUE, returnsamp = TRUE)

test2 <- nimble.occup(siteformula = ~ 1, ~ 1, y = y, sitevars = NULL, obsvars = NULL, mixture = "Bernoulli", priors = "Normal",
                     dropbase = FALSE, niter = 100, burnin = 1, initmcmc = 1, chains = 2, returncode = TRUE, returnsamp = TRUE)



