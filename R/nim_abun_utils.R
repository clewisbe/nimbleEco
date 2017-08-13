### Utility Functions for Abundance N Mixture Models and Single Season Occupancy  ###

makeBUGSmodule <- function(fun) {
  ans <- structure(list(process = fun), class = "BUGSmodule")
  ans
}

# Helper function to gather data
gatherfctn <- function(mydata, key.col, val.col, gather.cols) {
  new.data <- tidyr::gather_(data = mydata, key_col = key.col, value_col = val.col, gather_cols = colnames(mydata)[gather.cols])
  return(new.data)
}


## Make Tidy Data From 3 Abundance Data Sets ## y is matrix, sitevars is matrix, obsvars is named list
maketidy <- function(y, sitevars = NULL, obsvars = NULL) {
  final.df <- list()

  if (is.null(y)) {
    stop("Count Data is Required")
  }
  y <- as.data.frame(y)
  L = dim(y)[2]
  VisitID <- paste("Visit", 1:L, sep = "")
  colnames(y) <- VisitID
  y$Site <- as.numeric(rownames(y))
  long.y <- tidyr::gather(y, Visit, Count, 1:L)


  # Site Variables
  if (!is.null(sitevars)) {
    sitevars <- as.data.frame(sitevars)
    sitevars$Site <- as.numeric(rownames(sitevars))

    # Check Dimensions
    if (dim(y)[1] != dim(sitevars)[1]) {
      stop("Number of Rows of Y and Site Variables Do Not Match")
    }

    # Join Response with Site Level Variables
    final.df[[1]] <- dplyr::left_join(long.y, sitevars, by = "Site")
  } else {
    final.df[[1]] <- long.y
  }

  # Combine Observation Level Data
  if (!is.null(obsvars)) {
    varnames.obs <- names(obsvars)

    for (i in seq_along(obsvars)) {
      mf = as.data.frame(obsvars[[i]])
      if (L != dim(mf)[2]) {
        stop("Number of Visits of Y and Obs Covariate Variables Do Not Match")
      }
      colnames(mf) <- VisitID
      mf$Site <- as.numeric(rownames(mf))
      final.df[[i + 1]] <- gatherfctn(mf, "Visit", varnames.obs[i], 1:L)
    }
  }
  return(Reduce(function(dtf1, dtf2) dplyr::left_join(dtf1, dtf2, by = c("Site", "Visit")), final.df))
}


# Match.Call Function that Picks Up Defaut Arguments
match.call.defaults <- function(...) {
  call <- evalq(match.call(expand.dots = FALSE), parent.frame(1))
  formals <- evalq(formals(), parent.frame(1))

  for (i in setdiff(names(formals), names(call))) call[i] <- list(formals[[i]])


  match.call(sys.function(sys.parent()), call)
}



## BUGS Code Generation Helper Functions ##

# Takes Y (as a name) and constructs Y[i]
LHS2BUGSterm <- function(LHScode, indexName) {
  substitute(X[I], list(X = LHScode, I = indexName))
}

# Make the Name X.effect/mu/tau
make.effect.name <- function(term, level) {
  if (level == quote(site)) {
    return(as.name(paste0(term, ".effect.s")))
  }
  if (level == quote(gam)) {
    return(as.name(paste0(term, ".effect.gam")))
  }
  if (level == quote(phi)) {
    return(as.name(paste0(term, ".effect.phi")))
  }
  if (level == quote(p)) {
    return(as.name(paste0(term, ".effect.p")))
  } else {
    return(as.name(paste0(term, ".effect.o")))
  }
}

# To avoid repeated names, seperate naming for Site and Obs level parameters
make.coef.name <- function(term, level) {
  if (as.character(level) == "site") {
    return(as.name(paste0("s.", term)))
  }
  if (as.character(level) == "gam") {
    return(as.name(paste0("gam.", term)))
  }
  if (as.character(level) == "phi") {
    return(as.name(paste0("phi.", term)))
  }
  if (as.character(level) == "p") {
    return(as.name(paste0("p.", term)))
  } else {
    as.name(paste0("o.", term))
  }
}

make.coef.obs <- function(term) {
  as.name(paste0("o.", term))
}

make.mean.name <- function(term) {
  as.name(paste0(term, ".mu"))
}

make.precision.name <- function(term) {
  as.name(paste0(term, ".tau"))
}

make.sigma.name <- function(term) {
  as.name(paste0(term, ".sigma"))
}

make.pred.name <- function(term) {
  as.name(paste0(term, ".pred"))
}

make.param.name <- function(paramval, xvalue) {
  substitute(B * X, list(B = paramval, X = xvalue))
}

# Make Normal Prior for Random Effect
make.random.prior <- function(term, mean, sig) {
  substitute(TERM ~ dnorm(MEAN, sd = SD), list(TERM = term, MEAN = mean, SD = sig))
}

# Function to Make Priors
make.prior <- function(term, prior, ...) {
  param = list(...)
  argnames <- names(param)
  if (!("location" %in% argnames)) {
    param$location = 0
  }
  if (!("scale" %in% argnames)) {
    param$scale = 1000
  }
  if (!("lb" %in% argnames)) {
    param$lb = 0
  }
  if (!("ub" %in% argnames)) {
    param$ub = 1000
  }
  if (!("df" %in% argnames)) {
    param$df = 5
  }
  if (!("ub" %in% argnames)) {
    param$ub = 1000
  }
  if (!("shape" %in% argnames)) {
    param$shape = 0.001
  }
  if (!("rate" %in% argnames)) {
    param$rate = 0.001
  }

  if (prior == "Normal") {
    substitute(TERM ~ dnorm(LOCATION, sd = SCALE), list(TERM = term, LOCATION = param$location, SCALE = param$scale))
  } else if (prior == "Uniform") {
    substitute(TERM ~ dunif(MIN, MAX), list(TERM = term, MIN = param$lb, MAX = param$ub))
  } else if (prior == "t") {
    substitute(TERM ~ dt(LOCATION, SCALE, DF), list(TERM = term, LOCATION = param$location, SCALE = param$scale, DF = param$df))
  } else if (prior == "Gamma") {
    substitute(TERM ~ dgamma(SHAPE, RATE), list(TERM = term, SHAPE = param$shape, RATE = param$rate))
  } else {
    substitute(TERM <- 0, list(TERM = term))
  }
}


# Lets Deparse Run Over Multiple Lines
safeDeparse <- function(expr) {
  ret <- paste(deparse(expr), collapse = "")
  # rm whitespace
  gsub("[[:space:]][[:space:]]+", " ", ret)
}

# Create A[1:L] where L is number of levels for factor
factor.bracket <- function(term, df) {
  if (is.factor(df[[term]])) {
    l = as.numeric(nlevels((df[[term]])))
    return(substitute(TERM[1:L], list(TERM = as.name(term), L = l)))
  }
}

# Pulls out formula from model frame
RHSForm <- function(form, as.form = FALSE) {
  rhsf <- form[[length(form)]]
  if (as.form)
    reformulate(deparse(rhsf)) else rhsf
}

# Remove Parentheses from term
removeParen <- function(term) {
  out = gsub("[()]", "", term)
  return(unlist(strsplit(out, split = "|", fixed = TRUE)))
}

# Adds Term to a BUGS RHS.
addTerm <- function(currentTerms, newTerm) {
  substitute(A + B, list(A = currentTerms, B = newTerm))
}

prodTerm <- function(currentTerms, newTerm) {
  substitute(A * B, list(A = currentTerms, B = newTerm))
}

# Make Factor Variable
makeFac <- function(x, char.only = FALSE) {
  if (!is.factor(x) && (!char.only || is.character(x)))
    factor(x) else x
}

# Put Block Indexing Around Term
bracket <- function(term, index) {
  substitute(TERM[INDEX], list(TERM = term, INDEX = index))
}

# Make Double Bracket for Factor A[i,j]
doublebracket <- function(term, index1, index2) {
  substitute(TERM[INDEX1, INDEX2], list(TERM = term, INDEX1 = index1, INDEX2 = index2))
}


# Removes Bracket Indexing for single [ ]
removeIndexing <- function(term) {
  return(gsub("\\[.*", "", term))
}

make.interaction <- function(interaction, factors, index.name) {
  if (interaction %in% factors) {
    out = substitute(TERM[INDEX], list(TERM = as.name(interaction), INDEX = index.name))
    out <- substitute((TERM), list(TERM = out))

  } else {
    out = substitute(TERM[INDEX], list(TERM = as.name(interaction), INDEX = index.name))
  }
  return(out)
}


no.groups <- function(term) {
  term <- as.character(term)
  splt <- as.numeric(strsplit(term, "[^0-9]+")[[1]])
  return(max(splt, na.rm = T) - min(splt, na.rm = T) + 1)
}


# This will take a list of code lines and embed them in {}
embedLinesInCurlyBrackets <- function(lines) {
  as.call(c(list(quote(`{`)), lines))
}


# Put Lines in a For Loop
embedLinesInForLoop <- function(lines, indexName, start = 1, finish) {
  linesInBrackets <- embedLinesInCurlyBrackets(lines)
  rangeCall <- substitute(A:B, list(A = start, B = finish))
  ans <- substitute(for (INDEXNAME in RANGE) STUFF, list(INDEXNAME = indexName, RANGE = rangeCall, STUFF = linesInBrackets))
  ans
}


# Find Fixed Effect Terms (from LME4 Package)
nobars <- function(term) {
  nb <- nobars_(term)  ## call recursive version
  if (is(term, "formula") && length(term) == 3 && is.symbol(nb)) {
    ## called with two-sided RE-only formula: construct response~1 formula
    nb <- reformulate("1", response = deparse(nb))
  }
  ## called with one-sided RE-only formula, or RHS alone
  if (is.null(nb)) {
    nb <- if (is(term, "formula"))
      ~1 else 1
  }
  nb
}

nobars_ <- function(term) {
  if (!anyBars(term))
    return(term)
  if (isBar(term))
    return(NULL)
  if (isAnyArgBar(term))
    return(NULL)
  if (length(term) == 2) {
    nb <- nobars_(term[[2]])
    if (is.null(nb))
      return(NULL)
    term[[2]] <- nb
    return(term)
  }
  nb2 <- nobars_(term[[2]])
  nb3 <- nobars_(term[[3]])
  if (is.null(nb2))
    return(nb3)
  if (is.null(nb3))
    return(nb2)
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}

isBar <- function(term) {
  if (is.call(term)) {
    if ((term[[1]] == as.name("|")) || (term[[1]] == as.name("||"))) {
      return(TRUE)
    }
  }
  FALSE
}

isAnyArgBar <- function(term) {
  if ((term[[1]] != as.name("~")) && (term[[1]] != as.name("("))) {
    for (i in seq_along(term)) {
      if (isBar(term[[i]]))
        return(TRUE)
    }
  }
  FALSE
}

anyBars <- function(term) {
  any(c("|", "||") %in% all.names(term))
}


# Find Random Effect Terms from Formula object
fb <- function(term) {
  if (is.name(term) || !is.language(term))
    return(NULL)
  if (identical(term[[1]], as.name("(")))
    return(fb(term[[2]]))
  stopifnot(is.call(term))
  if (term[[1]] == as.name("|"))
    return(term)
  if (length(term) == 2)
    return(fb(term[[2]]))
  c(fb(term[[2]]), fb(term[[3]]))
}

# Make nim_glm object to pass to BUGS Module for Processing
make.glm <- function(RHS, factors = "None", cl, mixture, link, N = NULL, Site = NULL, level = NULL, prior.params = NULL) {
  RHS.out <- call("nim_glm", RHSForm(RHS))
  RHS.out[[3]] <- factors
  RHS.out[[4]] <- as.name(cl$priors)
  RHS.out[[5]] <- prior.params
  RHS.out[[6]] <- cl$dropbase
  RHS.out[[7]] <- cl$droplast
  RHS.out[[8]] <- mixture
  RHS.out[[9]] <- link
  RHS.out[[10]] <- level
  names(RHS.out)[3:10] <- c("factors", "priors", "prior.shape", "dropbase", "droplast", "family", "link", "level")
  if (mixture %in% c("Binomial", "Bernoulli")) {
    RHS.out[[11]] <- list(N, Site)
    names(RHS.out)[11] <- c("args")
  }
  return(RHS.out)
}

## nim_glm BUGS Module ## Returns Full Code, LHS, RHS, and the Probability Model
nim_glm <- makeBUGSmodule(function(LHS, RHS) {
  RHSargs <- match.call(function(mod, factors, priors, prior.shape, family, link, dropbase, droplast, level, args) {
  }, RHS)

  # Set Up LHS
  index.name <- quote(i)

  if (RHSargs$family %in% c("Poisson", "NB", "ZIP", "Bernoulli")) {
    distn <- quote(dpois)
    link = quote(log)
    param = quote(lambda)
    pred.d <- substitute(DIST(PARAM), list(DIST = distn, PARAM = LHS2BUGSterm(param, index.name)))
  }

  if (RHSargs$family == "Binomial") {
    distn <- quote(dbinom)
    link = quote(logit)
    param = quote(p)
    pred.d <- substitute(DIST(PARAM, N), list(DIST = distn, PARAM = LHS2BUGSterm(param, index.name), N = bracket(RHSargs$args[[1]], bracket(RHSargs$args[[2]], quote(i)))))
  }

  meanfctn <- substitute(LHS ~ RHS, list(LHS = LHS2BUGSterm(LHS[[2]], index.name), RHS = pred.d))
  meanfctn.loop <- embedLinesInForLoop(meanfctn, index.name, start = 1, finish = LHS[[3]][[3]])

  if (RHSargs$family == "NB") {
    dispersion <- substitute(LAMBDA <- MU * DISP, list(LAMBDA = LHS2BUGSterm(param, index.name), DISP = LHS2BUGSterm(quote(rho), index.name), MU = LHS2BUGSterm(quote(mu), index.name)))
    rho.dist <- substitute(X ~ dgamma(alpha, alpha), list(X = LHS2BUGSterm(quote(rho), index.name)))
    meanfctn.no.prior <- embedLinesInForLoop(list(meanfctn, dispersion, rho.dist), index.name, start = 1, finish = LHS[[3]][[3]])

    rho.param <- substitute(ALPHA <- exp(LOGALPHA), list(ALPHA = quote(alpha), LOGALPHA = quote(logalpha)))
    rho.prior <- make.prior(quote(logalpha), "Normal")

    meanfctn.loop <- embedLinesInCurlyBrackets(list(rho.prior, rho.param, meanfctn.no.prior))
  }

  if (RHSargs$family == "ZIP") {
    zip <- substitute(S ~ dbern(THETA), list(S = LHS2BUGSterm(quote(s), index.name), THETA = quote(theta)))
    zip.prior <- make.prior(quote(theta), "Uniform", ub = 1)
    distn <- quote(dpois)
    link = quote(log)
    param = quote(lambda)
    pred.d <- substitute(DIST(PARAM * ZIP), list(DIST = distn, PARAM = LHS2BUGSterm(param, index.name), ZIP = LHS2BUGSterm(quote(s), index.name)))
    meanfctn <- substitute(LHS ~ RHS, list(LHS = LHS2BUGSterm(LHS[[2]], index.name), RHS = pred.d))
    meanfctn.no.prior <- embedLinesInForLoop(list(meanfctn, zip), index.name, start = 1, finish = LHS[[3]][[3]])
    meanfctn.loop <- embedLinesInCurlyBrackets(list(zip.prior, meanfctn.no.prior))
  }

  if (RHSargs$family == "Bernoulli") {
    if (RHSargs$level == quote(site)) {
      param <- quote(psi)
      presence <- substitute(dbern(PSI), list(PSI = LHS2BUGSterm(param, index.name)))
    }
    if (RHSargs$level == quote(phi)) {
      param <- quote(phi)
      presence <- substitute(dbern(PHI), list(PHI = LHS2BUGSterm(param, index.name)))
    }
    if (RHSargs$level == quote(p)) {
      param <- quote(p)
      presence <- substitute(dbern(P), list(P = LHS2BUGSterm(param, index.name)))
    }

    if (RHSargs$level == quote(gam)) {
      param <- quote(gamma)
      presence <- substitute(dbern(GAM), list(GAM = LHS2BUGSterm(param, index.name)))
    }
    if (RHSargs$level == quote(obs)) {
      param <- quote(p)
      presence <- substitute(dbern(Z * P), list(Z = bracket(quote(z), bracket(RHSargs$args[[2]], quote(i))), P = LHS2BUGSterm(param, index.name)))
    }
    meanfctn <- substitute(LHS ~ RHS, list(LHS = LHS2BUGSterm(LHS[[2]], index.name), RHS = presence))
    meanfctn.loop <- embedLinesInForLoop(meanfctn, index.name, start = 1, finish = LHS[[3]][[3]])
  }

  # Set up RHS
  lpred <- call("lmPred", RHSargs$mod)
  lpred[[3]] <- RHSargs$factors
  lpred[[4]] <- RHSargs$priors
  lpred[[5]] <- RHSargs$prior.shape
  lpred[[6]] <- as.name(RHSargs$link)
  lpred[[7]] <- RHSargs$dropbase
  lpred[[8]] <- RHSargs$droplast
  lpred[[9]] <- RHSargs$level
  names(lpred)[3:9] <- c("factors", "priors", "prior.shape", "link", "dropbase", "droplast", "level")

  if (RHSargs$family == "NB") {
    lpred.full <- substitute(LHS <- RHS, list(LHS = bracket(quote(mu), LHS[[3]]), RHS = lpred))
    LHS.Param <- bracket(quote(mu), LHS[[3]])
  } else {
    lpred.full <- substitute(LHS <- RHS, list(LHS = bracket(param, LHS[[3]]), RHS = lpred))
    LHS.Param <- bracket(param, LHS[[3]])
  }

  # Return Code
  newCode <- embedLinesInCurlyBrackets(lines = list(meanfctn.loop, lpred.full))
  return(list(code = newCode, LHS = LHS.Param, RHS = lpred, prob.mod = meanfctn.loop))
})


# Takes Term from lmPred and Converts it to a BUGS Term Note: : is not supported for interaction terms like for lm function, just *
lmPred2BUGSterm <- function(term, factors, index.name, level) {


  factors.name <- lapply(factors, removeIndexing)
  factors.name <- unlist(lapply(factors.name, "[[", 2))
  if (term == "1") {
    return(make.coef.name(quote(intercept), level))
  }
  if (length(strsplit(term, "*", fixed = TRUE)[[1]]) == 2) {
    interaction = unlist(strsplit(term, "*", fixed = TRUE))

    if (interaction[1] %in% factors.name & !(interaction[2] %in% factors.name)) {
      out = substitute(EFFECTNAME[TERM[INDEX]] * TERM2[INDEX], list(EFFECTNAME = make.effect.name(paste0(interaction[1], interaction[2]), level), TERM = as.name(interaction[1]),
                                                                    TERM2 = as.name(interaction[2]), INDEX = index.name))
      return(out)
    }

    if (interaction[2] %in% factors.name & !(interaction[1] %in% factors.name)) {
      out = substitute(EFFECTNAME[TERM[INDEX]] * TERM2[INDEX], list(EFFECTNAME = make.effect.name(paste0(interaction[1], interaction[2]), level), TERM = as.name(interaction[2]),
                                                                    TERM2 = as.name(interaction[1]), INDEX = index.name))
      return(out)
    }

    if (all(interaction %in% factors.name)) {
      out = substitute(EFFECTNAME[TERM[INDEX], TERM2[INDEX]], list(EFFECTNAME = make.effect.name(paste0(interaction[1], interaction[2]), level), TERM = as.name(interaction[2]),
                                                                   TERM2 = as.name(interaction[1]), INDEX = index.name))
      return(out)
    } else {
      out = substitute(EFFECTNAME * TERM[INDEX] * TERM2[INDEX], list(EFFECTNAME = make.coef.name(paste0(interaction[1], interaction[2]), level), TERM = as.name(interaction[1]),
                                                                     TERM2 = as.name(interaction[2]), INDEX = index.name))
      return(out)
    }
  }

  if (length(strsplit(term, "*", fixed = TRUE)[[1]]) > 2) {
    interaction = unlist(strsplit(term, "*", fixed = TRUE))

    tm.name <- make.coef.name(paste(interaction, collapse = ""), level)

    if (interaction[1] %in% factors) {
      out = substitute(EFFECTNAME[TERM[INDEX]], list(EFFECTNAME = tm.name, TERM = as.name(interaction[1]), INDEX = index.name))
    } else {
      out = substitute(EFFECTNAME * TERM[INDEX], list(EFFECTNAME = tm.name, TERM = as.name(interaction[1]), INDEX = index.name))
    }

    for (i in seq_along(interaction[-1])) {
      out <- prodTerm(out, make.interaction(interaction[[i + 1]], factors, index.name))
    }
    return(out)
  }

  if (term %in% factors.name) {
    ## If term is a factor use BUGS indexing create b2[block[i]]
    out = substitute(EFFECTNAME[TERM[INDEX]], list(EFFECTNAME = make.effect.name(term, level), TERM = as.name(term), INDEX = index.name))
    return(out)
  } else {
    ## otherwise create BUGS term like x[i]
    out = substitute(TERM[INDEX], list(TERM = as.name(term), INDEX = index.name))

    make.param.name(make.coef.name(term, level), out)
  }

}

# lmPred to BUGS Expansion for Random Effect Term
lmPred2BUGSRandom <- function(termparen, factors, index.name, level) {

  term <- removeParen(termparen)

  factors.name <- lapply(factors, removeIndexing)
  factors.name <- unlist(lapply(factors.name, "[[", 2))

  # Random Intercept
  if (term[[1]] == 1 & removeIndexing(term[[2]]) %in% factors.name) {

    out = substitute(EFFECTNAME[TERM[INDEX]], list(EFFECTNAME = make.coef.name(removeIndexing(term[[2]]), level), TERM = as.name(removeIndexing(term[[2]])), INDEX = index.name))
  }

  # Random Slope
  if (removeIndexing(term[[1]]) != 1 & removeIndexing(term[[2]]) %in% factors.name) {
    out = substitute(EFFECTNAME[TERM[INDEX]] * TERM2, list(EFFECTNAME = make.effect.name(paste0(term[1], term[2]), level), TERM = as.name(removeIndexing(term[[2]])), TERM2 = LHS2BUGSterm(as.name(removeIndexing(term[[1]])),
                                                                                                                                                                                           index.name), INDEX = index.name))

  }

  if (!(term[[2]] %in% factors.name)) {
    stop("Conditioned Term is Not a Factor")
  }
  return(out)
}

# Make Prior for Fixed Term
lmPredTerm2PRIOR <- function(term, RHSname, factors, index.name, prior, prior.shape, dropbase, droplast, int.check) {
  priors = list()

  factors.name <- lapply(factors, removeIndexing)
  factors.name <- unlist(lapply(factors.name, "[[", 2))


  if (length(strsplit(as.character(RHSname), "*", fixed = TRUE)[[1]]) > 0) {
    param.name <- as.name(removeIndexing((gsub(" ", "", strsplit(as.character(RHSname), "*", fixed = TRUE)[[1]][[1]]))))
  } else {
    param.name <- as.name(removeIndexing(RHSname))
  }

  if (term %in% factors.name) {
    ## if term is a factor, create priors for each level (note, dropbase=T will assign 0 prior to first level)

    fct.index = factors[which(as.character(term) == factors.name)]
    groups = no.groups(fct.index)  #Drop one group for model identification; Need to add Check if Intercept in Model Too

    if (dropbase == TRUE | int.check == 1) {
      priors[[1]] <- make.prior(LHS2BUGSterm(param.name, 1), prior = 0)
      st = 2
      idx = 2
    } else {
      st = 1
      idx = 1
    }

    p = do.call(make.prior, c(list(term = quote(LHS2BUGSterm(param.name, index.name)), prior = as.character(prior)), prior.shape))

    #Set Last Group to 0 if droplast = T (for dynamic occupancy models)
    if (droplast == TRUE) {
      groups = groups - 1
      idx.last = idx + 1
      priors[[idx.last]] <- make.prior(LHS2BUGSterm(param.name, groups + 1), prior = 0)
    }

    priors[[idx]] <- embedLinesInForLoop(p, index.name, start = st, finish = groups)

    initial <- list(groups)
    names(initial) <- quote(param.name)
    return(c(embedLinesInCurlyBrackets(priors), initial))
  }

  # Priors for Interaction Terms with Factor Variable
  if (length(strsplit(as.character(term), "*", fixed = TRUE)[[1]]) == 2) {
    interaction = unlist(strsplit(as.character(term), "*", fixed = TRUE))
    # Taking this Piece Our for Now.  Determining Sum to 0 constraints for 2 and 3 Level Anova gets messy quick.  if(interaction[1] %in% factors.name & interaction[2] %in%
    # factors.name){ #Do We want to Allow Users to Have Interaction Terms for Multiple Factors?  out <- list() group1 <- no.groups(factors[which(interaction[1] %in% factors.name)])
    # group2 <- no.groups(factors[which(interaction[2] %in% factors.name)]) if(dropbase==TRUE | int.check==1){ #Identifiability Constraints out[[1]] <-
    # make.prior(doublebracket(param.name, 1, 1), prior = 0) out[[2]] <- embedLinesInForLoop(make.prior(doublebracket(param.name, 1, index.name), prior = 0), index.name, start = 2,
    # finish = group1) out[[3]] <- embedLinesInForLoop(make.prior(doublebracket(param.name, index.name, 1), prior = 0), index.name, start = 2, finish = group2) out[[4]] <-
    # embedLinesInForLoop(embedLinesInForLoop(make.prior(doublebracket(param.name, index.name, quote(j)), prior = 'Normal'), index.name, start = 2, finish = group2), quote(j), start =
    # 2, finish = group1) } return(out) }


    if (interaction[1] %in% factors.name | interaction[2] %in% factors.name) {

      # Assume the Factor Group will Have a Lower Upper Index Number; Think About Better Way To Grab This
      fct.index = factors[which(interaction == factors.name)]
      groups = no.groups(fct.index)

      if (dropbase == TRUE) {
        priors[[1]] <- make.prior(LHS2BUGSterm(param.name, 1), prior = 0)
        st = 2
        idx = 2
      } else {
        st = 1
        idx = 1
      }

      p = do.call(make.prior, c(list(term = quote(LHS2BUGSterm(param.name, index.name)), prior = as.character(prior)), prior.shape))

      #Set Last Group to 0 if droplast = T (for dynamic occupancy models)
      if (droplast == TRUE) {
        groups = groups - 1
        idx.last = idx + 1
        priors[[idx.last]] <- make.prior(LHS2BUGSterm(param.name, groups + 1), prior = 0)
      }

      priors[[idx]] <- embedLinesInForLoop(p, index.name, start = st, finish = groups)
      initial <- list(groups)
      names(initial) <- quote(param.name)
      return(c(embedLinesInCurlyBrackets(priors), initial))
    } else {
      priors = do.call(make.prior, c(list(term = quote(param.name), prior = as.character(prior)), prior.shape))
      initial = list(1)
      names(initial) <- quote(param.name)
      return(c(priors, initial))
    }
  } else {
    priors = do.call(make.prior, c(list(term = quote(param.name), prior = as.character(prior)), prior.shape))
    initial = list(1)
    names(initial) <- quote(param.name)
    return(c(priors, initial))
  }
}

# Prior for RANDOM EFFECTS Term
lmPredRandom2PRIOR <- function(termparen, RHSname, factors, index.name, prior, prior.shape, dropbase, droplast, int.check) {

  priors = list()
  factors.name <- lapply(factors, removeIndexing)
  factors.name <- unlist(lapply(factors.name, "[[", 2))

  term <- removeParen(termparen)
  groups <- no.groups(factors[which(term[2] == factors.name)])
  param.name <- as.name(removeIndexing(RHSname))

  # Priors for Random Intercept
  if (term[1] == 1) {
    if (dropbase == TRUE | int.check == 1) {
      random <- make.random.prior(LHS2BUGSterm(param.name, index.name), mean = make.mean.name(param.name), sig = make.sigma.name(param.name))
      prior.mu <- make.prior(make.mean.name(param.name), 0)
    } else {
      random <- make.random.prior(LHS2BUGSterm(param.name, index.name), mean = make.mean.name(param.name), sig = make.sigma.name(param.name))
      prior.mu = do.call(make.prior, c(list(term = quote(make.mean.name(param.name)), prior = as.character(prior)), prior.shape))
    }

    #Set Last Group to 0 if droplast = T (for dynamic occupancy models)
    if (droplast == TRUE) {
      groups = groups - 1
      last.group.prior <- make.prior(LHS2BUGSterm(param.name, groups + 1), prior = 0)
    }

    random.loop <- embedLinesInForLoop(random, index.name, start = 1, finish = groups)

    # Priors
    prior.sd <- make.prior(make.sigma.name(param.name), "Uniform", lb = 0, ub = 1000)

    if(droplast == TRUE){
      priors <- embedLinesInCurlyBrackets(list(random.loop, prior.sd, prior.mu, last.group.prior))
    } else{
      priors <- embedLinesInCurlyBrackets(list(random.loop, prior.sd, prior.mu))
    }

    initial <- list(1, 1)
    initial <- setNames(initial, c(as.character(make.sigma.name(param.name)), as.character(make.mean.name(param.name))))
    # initial[[1]] <-cbind(as.character(make.sigma.name(param.name)), 1, as.character(make.mean.name(param.name)), 1 )
    return(c(priors, initial))
  } else {
    random <- make.random.prior(LHS2BUGSterm(param.name, index.name), mean = make.mean.name(param.name), sig = make.sigma.name(param.name))
    prior.mu <- make.prior(make.mean.name(param.name), 0)

    if (droplast == TRUE) {
      groups = groups - 1
      last.group.prior <- make.prior(LHS2BUGSterm(param.name, groups + 1), prior = 0)
    }

    last.group.prior <- invisible()
    random.loop <- embedLinesInForLoop(random, index.name, start = 1, finish = groups)

    # Priors
    prior.sd <- make.prior(make.sigma.name(param.name), "Uniform", lb = 0, ub = 1000)
    if(droplast == TRUE){
      priors <- embedLinesInCurlyBrackets(list(random.loop, prior.sd, prior.mu, last.group.prior))
    } else{
      priors <- embedLinesInCurlyBrackets(list(random.loop, prior.sd, prior.mu))
    }

    initial <- list(1, 1)
    initial <- setNames(initial, c(as.character(make.sigma.name(param.name)), as.character(make.mean.name(param.name))))
    return(c(priors, initial))
  }
}


## Make lmPredictor function ## Returns Full BUGS code and List of Parameters and Initital Values for MCMC
lmPredictor <- makeBUGSmodule(function(LHS, RHS) {

  RHSargs <- match.call(function(RHSmodel, factors, priors, prior.shape, family, link, level, dropbase, droplast) {
  }, RHS)

  link <- match.arg(as.character(RHSargs$link), c("identity", "log", "logit", "Poisson", "Bernoulli"))

  LHS.Index <- bracket(as.name(removeIndexing(LHS)[[2]]), quote(i))


  # LHS with Link
  if (link == "identity") {
    LHS.out <- LHS.Index
  }

  if (link == "log") {
    LHS.out <- substitute(LINK(LHS), list(LINK = quote(log), LHS = LHS.Index))
  }

  if (link == "logit") {
    LHS.out <- substitute(LINK(LHS), list(LINK = quote(logit), LHS = LHS.Index))
  }


  terms <- strsplit(safeDeparse(RHSargs$RHSmodel), split = "+", fixed = TRUE)[[1]]
  terms.list <- lapply(terms, FUN = function(x) {
    gsub(" ", "", x, fixed = TRUE)
  })


  RHS.fix.tm <- list()
  RHS.random.tm <- list()

  for (iTerm in seq_along(terms.list)) {
    if (grepl("|", terms.list[[iTerm]], fixed = TRUE) == FALSE) {
      RHS.fix.tm[[iTerm]] <- terms.list[[iTerm]]
    } else {
      RHS.random.tm[[iTerm]] <- terms.list[[iTerm]]
    }
  }

  # Drop Null Slots.
  RHS.fix.tm <- RHS.fix.tm[!sapply(RHS.fix.tm, is.null)]
  RHS.random.tm <- RHS.random.tm[!sapply(RHS.random.tm, is.null)]

  #Make Intercept a Default Option
  RM <- FALSE
  if(any(sapply(RHS.fix.tm, function(x)  -1 %in% x)) == TRUE){
    RHS.fix.tm <- RHS.fix.tm[-sapply(RHS.fix.tm, function(x)  -1 %in% x)]
    RM <- TRUE
  }

  if(any(sapply(RHS.fix.tm, function(x)  1 %in% x)) == FALSE & RM ==FALSE){
    RHS.fix.tm <- append(RHS.fix.tm, "1")
  }



  index.name <- quote(i)
  int.check <- sum(!(is.na(sapply(RHS.fix.tm, function(x) pmatch("1", x)))))


  # Generate Terms and Priors for Fixed Terms
  RHS.fix <- list()
  Priors.fix <- list()

  for (iTerm in seq_along(RHS.fix.tm)) {
    RHS.fix[[iTerm]] <- lmPred2BUGSterm(RHS.fix.tm[[iTerm]], RHSargs$factors, index.name, RHSargs$level)
    Priors.fix[[iTerm]] <- lmPredTerm2PRIOR(RHS.fix.tm[[iTerm]], RHS.fix[iTerm], factors = RHSargs$factors, index.name, prior = RHSargs$priors, prior.shape = RHSargs$prior.shape,
                                            dropbase = RHSargs$dropbase, droplast = RHSargs$droplast, int.check = int.check)
  }

  # Generate Terms and Priors for Random Terms
  RHS.random <- list()
  Priors.random <- list()
  for (iTerm in seq_along(RHS.random.tm)) {
    RHS.random[[iTerm]] <- lmPred2BUGSRandom(RHS.random.tm[[iTerm]], RHSargs$factors, index.name, RHSargs$level)
    Priors.random[[iTerm]] <- lmPredRandom2PRIOR(RHS.random.tm[[iTerm]], RHS.random[iTerm], factors = RHSargs$factors, index.name, prior = RHSargs$priors, prior.shape = RHSargs$prior.shape,
                                                 dropbase = RHSargs$dropbase, droplast = RHSargs$droplast, int.check = int.check)
  }


  Priors.Combo <- append(lapply(Priors.fix, "[[", 1), lapply(Priors.random, "[[", 1))
  Terms.Combo <- append(RHS.fix, RHS.random)

  # Initial.Values <- append(lapply(Priors.fix, '[[', 2), lapply(Priors.random, '[[', 2))
  Inits.Fix <- sapply(Priors.fix, function(x) x[2])
  Inits.Random <- append(sapply(Priors.random, function(x) x[2]), sapply(Priors.random, function(x) x[3]))
  Inits.All <- append(Inits.Fix, Inits.Random)

  RHS2 <- Terms.Combo[[1]]
  for (i in seq_along(Terms.Combo)[-1]) {
    RHS2 <- addTerm(RHS2, Terms.Combo[[i]])
  }

  # Make Model Form and Embed in Loop
  fullLine <- substitute(LHS <- RHS, list(LHS = LHS.out, RHS = RHS2))
  forLoop <- embedLinesInForLoop(fullLine, index.name, start = 1, finish = LHS[[3]][[3]])
  Priors.All <- embedLinesInCurlyBrackets(Priors.Combo)

  newCode <- embedLinesInCurlyBrackets(list(forLoop, Priors.All))

  return(list(code = newCode, inits = Inits.All))
})


## Functions to Get Code, Data, and Initial Values ready for NIMBLE ##

# Make Starting Values for a Variable Passed To nimbleModel
make.init <- function(term, initialvalue, size) {
  l = rep(initialvalue, size)
  init <- substitute(c(TERM, L), list(TERM = as.character(term), L = l))
  init
}


lmPred2init <- function(term, value) {
  out = rep(value, term)
  term <- out
  return(term)
}

## Functions that Interact with NIMBLE ##

# Call Nimble and Pass the BUGS code and data
callNim <- function(modobj, initial, niter, burn, chains) {
  nimMod.C <- nimble::compileNimble(modobj)

  # Set Up MCMC
  config.Mod <- nimble::configureMCMC(nimMod.C, print = FALSE)
  mod.MCMC <- nimble::buildMCMC(config.Mod)
  C.mod.MCMC <- nimble::compileNimble(mod.MCMC)
  samplesList <- nimble::runMCMC(C.mod.MCMC, inits = initial, niter = niter, nburnin = burn, nchains = chains, returnCodaMCMC = TRUE)
  return(samplesList)
}



# Take MCMC coda object and return summary statistics
makeoutput <- function(codaobj, sitevars, obsvars) {
  quantile <- summary(codaobj)$quantiles
  meansd <- summary(codaobj)$statistics[, 1:2]
  rhat <- tryCatch(coda::gelman.diag(codaobj, multivariate = FALSE)$psrf, error = function(e) NA)
  Effsamp <- coda::effectiveSize(codaobj)
  output <- as.data.frame(cbind(meansd, quantile, rhat, Effsamp))
  names(output)[names(output) == "Point est."] <- "Rhat"
  output$`Upper C.I.` <- NULL
  return(as.data.frame(output))
}

