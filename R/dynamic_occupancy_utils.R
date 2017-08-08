### NIMBLE Dynamic Occupancy Models ###

library(nimble)
# library(devtools) install_github('nimble-dev/nimble', ref = 'devel', subdir = 'packages/nimble')
library(tidyverse)
library(reshape2)
library(coda)
nimbleOptions(enableBUGSmodules = TRUE)

makeBUGSmodule <- function(fun) {
  ans <- structure(list(process = fun), class = "BUGSmodule")
  ans
}

# Takes 4 Data Objects and Returns Tidy Data Frame for Dynamic Occupancy Models Data should be entered as follows: y is array (site, survey, season), site is matrix or df,
# season.site is named list (each element a variable), and season.site.survey is a named list of arrays (site, survey, season)

tidy.dynam.occ <- function(y = NULL, site = NULL, season.site = NULL, season.site.survey = NULL) {

  if (is.null(y)) {
    stop("Must Provide 0/1 Observation Data, y")
  }

  # Get Number of Seasons
  nseason <- as.numeric(dim(y)[3])

  # Make Tidy Data for Y
  y.tidy <- melt(y, varnames = c("Site", "Survey", "Season"), value.name = "y")
  y.tidy$Site <- y.tidy$Site
  y.tidy$Season <- as.factor(y.tidy$Season)

  # Combine Site/Season Data
  site.season.groups <- y.tidy[y.tidy$Survey == 1, c(1, 3)]
  site.season.groups$Season <- as.factor(site.season.groups$Season)
  site.season.groups$Site <- site.season.groups$Site

  if (!is.null(site)) {
    covariates <- list()
    site$Site <- as.numeric(site.season.groups[site.season.groups$Season == 1, 1])
  }


  if (!is.null(season.site)) {
    covariates <- list()
    covariates[[1]] <- site.season.groups

    site.season.vars <- names(season.site)
    for (i in seq_along(season.site)) {
      mf = as.data.frame(season.site[[i]])
      if (nseason != dim(mf)[2]) {
        stop("Number of Seasons do Not Match Between Y and Covariates")
      }
      colnames(mf) <- seq(1, nseason, 1)
      mf$Site <- as.numeric(rownames(mf))
      covariates[[i + 1]] = melt(mf, id.vars = "Site", variable.name = "Season", value.name = site.season.vars[i])
    }

    season.site.tidy <- Reduce(function(dtf1, dtf2) left_join(dtf1, dtf2, by = c("Site", "Season")), covariates)
  }

  if (!is.null(season.site.survey)) {
    covariates <- list()
    covariates[[1]] <- as.data.frame(cbind(y.tidy$Site, y.tidy$Survey, y.tidy$Season))
    colnames(covariates[[1]]) <- c("Site", "Survey", "Season")
    obs.vars <- names(season.site.survey)

    for (i in seq_along(season.site.survey)) {
      if (is.array(season.site.survey[[i]]) == FALSE) {
        stop("Site Level Covariates Need to be Provided as An Array")
      }
      mf <- melt(season.site.survey[[i]], varnames = c("Site", "Survey", "Season"), value.name = obs.vars[i])
      if (dim(y.tidy)[1] != dim(mf)[1]) {
        stop("Number of Observations do Not Match Between Y and Covariates")
      }
      covariates[[i + 1]] = mf
    }
    season.site.survey.tidy <- Reduce(function(dtf1, dtf2) left_join(dtf1, dtf2, by = c("Site", "Survey", "Season")), covariates)
    season.site.survey.tidy$Season <- as.factor(season.site.survey.tidy$Season)
  } else {
    season.site.survey.tidy = NULL
  }


  if (!is.null(season.site.survey.tidy)) {
    final.tidy1 <- Reduce(function(dtf1, dtf2) left_join(dtf1, dtf2, by = c("Site", "Survey", "Season")), list(y.tidy, season.site.survey.tidy))
  } else {
    final.tidy1 <- y.tidy
  }

  if (!is.null(season.site)) {
    final.tidy2 <- Reduce(function(dtf1, dtf2) left_join(dtf1, dtf2, by = c("Site", "Season")), list(final.tidy1, season.site.tidy))
  } else {
    final.tidy2 <- final.tidy1
  }

  if (!is.null(site)) {
    final.tidy3 <- Reduce(function(dtf1, dtf2) left_join(dtf1, dtf2, by = c("Site")), list(final.tidy2, site))
  } else {
    final.tidy3 <- final.tidy2
  }

  # Ensure Site/Survey/Season are Factors
  cols <- c("Site", "Survey", "Season")
  final.tidy3[cols] <- lapply(final.tidy3[cols], factor)

  return(final.tidy3)
}



latent.mean.loop <- function(nsite, nseason) {
  fin = nseason - 1
  l1 <- substitute(muZ[(T * (s - 1) + 1):(T * s)] <- z[(T * (s - 1) + 1):(T * s)] * phi[(T * (s - 1) + 1):(T * s)] +
                      (1 - z[(T * (s - 1) + 1):(T * s)]) * gamma[(T * (s - 1) + 1):(T *s)], list(T = nsite))
  return(embedLinesInForLoop(l1, quote(s), start = 1, finish = fin))
}


latent.state.loop <- function(nsite, nseason, index) {
  end <- nsite * nseason - nsite
  state <- substitute(z[INDEX + NSITE] ~ DIST(muZ[INDEX]), list(INDEX = index, NSITE = nsite, DIST = quote(dbern)))
  return(embedLinesInForLoop(state, index, start = 1, finish = end))
}

make.obs.loop <- function(nsite, nseason, nsurvey) {
  l1 <- quote(muy[i] <- z[SiteSeason[i]] * p[i])  #Indexing for Latent State Since We Need to Loop Over Surveys
  l2 <- substitute(y[i] ~ DIST(muy[i]), list(DIST = quote(dbern)))
  return(embedLinesInForLoop(list(l1, l2), quote(i), start = 1, finish = nsite * nseason * nsurvey))
}
