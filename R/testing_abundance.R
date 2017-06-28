#Abundance
library(unmarked)
library(reshape2)
library(tidyr)
library(dplyr)
library(nimble)

#Mallard Data from Unmarked Packcage
#Mallard Data from Unmarked Packcage
y <- as.data.frame(mallard.y)
sitecov <- as.data.frame(mallard.site)
obscov <- mallard.obs



test <- nimble.abund(~ 1 + elev + (1|Visit), ~ date, y = y, sitevars = sitecov, obsvars = obscov, mixture = "Poisson", priors = "Normal",
                     dropbase = FALSE, niter = 10000, burnin = 1000, initmcmc = 1, chains = 1, returncode = FALSE, returnsamp = FALSE)



#Getting Error C_setGraph not found##
