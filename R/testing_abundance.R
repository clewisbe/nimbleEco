#Abundance
library(unmarked)
library(reshape2)
library(tidyr)
library(dplyr)

#Mallard Data from Unmarked Packcage
y <- as.data.frame(mallard.y)
site <- as.data.frame(mallard.site)
obs <- mallard.obs



test <- nimble.abund(~ 1 + (1|ivel) ~ 1 + length + (length|forest), y = y, siteCov = site, ObsCov = obs, mixture = "Poisson", priors = "Normal",
                     dropbase = TRUE, niter = 10000, burnin = 1000, initmcmc = 1, chains = 1, returncode = FALSE, returnsamp = FALSE)






## IGNORE BELOW ##

visitid <- colnames(y)

#Make Site ID Variable
y$site <- as.numeric(rownames(y))
sitecov$site <- as.numeric(rownames(sitecov))

long.y <- gather(y, visit, count, 1:3)

site.merge <- left_join(long.y, sitecov, by = "site")

#Add Obs Level Covariates
evel = as.data.frame(obscov$ivel)
colnames(evel) <- visitid
evel$site <- as.numeric(rownames(evel))
long.evel <- gather(evel, visit, ivel, 1:3)

date = as.data.frame(obscov$date)
colnames(date) <- visitid
date$site <- as.numeric(rownames(date))
long.date <- gather(date, visit, date, 1:3)


out <- left_join(site.merge, long.evel, by = c("site", "visit"))
longdata <- left_join(out, long.date, by = c("site", "visit"))
head(longdata)
