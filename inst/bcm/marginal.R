# Models whose likelihood is a FINITE SUM over something discrete that
# was never observed: how many surveys were sent, when the mean changed,
# which country each person and each question came from.
#
# Every one of them is the structured protocol's whole-response loglik
# slot, and every one of them is the same shape as what the Stan ports
# do by hand: build the log probability of each value of the discrete
# thing, then log_sum_exp them. The difference is that Stan has to do it
# because it cannot sample a discrete parameter, while frmtmb has to do
# it because a maximum over a discrete parameter is not an estimate.
#
# See dev/structured-family-protocol.md for the slot contracts.
# vignette("bayesian-cognitive-modeling") shows these in place;
# tests/testthat/test-bcm-binomial.R,
# tests/testthat/test-bcm-data-analysis.R and
# tests/testthat/test-bcm-latent-mixtures.R fit this code. Load it with:
#
#   source(system.file("bcm", "marginal.R", package = "frmtmb"))

## ---- bcm-logsumexp ----
# Fold a growing set of log terms, on the tape. RTMB::logspace_add() is
# the pairwise operation; starting from NULL rather than from -Inf keeps
# the first term out of a subtraction with an infinity, which is the
# NaN this arithmetic is prone to.
bcm_lse_start <- function() NULL

bcm_lse_add <- function(acc, term) {
  if (is.null(acc)) term else RTMB::logspace_add(acc, term)
}

# Broadcast a dpar the core handed back as one value; see
# inst/rl/rw-delta.R for why rep(v, each = ) cannot be used here.
bcm_mbcast <- function(v, n) if (length(v) == 1L) v * rep(1, n) else v

## ---- bcm-survey ----
# Lee and Wagenmakers chapter 3.6. Five surveys were sent out and 16,
# 18, 22, 25 and 27 came back. Neither the number sent nor the return
# rate is known, and the number sent is a whole number with a uniform
# prior on 1..nmax.
#
# So the likelihood is
#
#   sum over n of  (1 / nmax) prod_j Binomial(k_j | n, theta)
#
# which is one number for the whole response, not a product over rows.
# The lower limit is max(k): fewer surveys than that cannot have come
# back, and those terms are zero rather than small, so they are dropped
# instead of being added as -Infinity.
bcm_survey <- function(nmax, link = "logit") {
  if (length(nmax) != 1L || nmax < 1 || nmax != round(nmax)) {
    stop("bcm_survey(nmax =) is the largest number of surveys that ",
         "could have been sent, one whole number", call. = FALSE)
  }
  fam <- frmtmb_family(
    "bcm_survey",
    dpars = "mu",
    links = list(mu = link),
    type = "discrete",
    lpdf = function(y, dpars, aterms, extra = NULL) {
      stop("A bcm_survey() return count has no row-wise log density: ",
           "every count is binomial in the SAME unknown number of ",
           "surveys, and that number is summed out of the whole ",
           "response at once. Use logLik() for the total, or ",
           "latent_probs() for the posterior over the number sent",
           call. = FALSE)
    },
    valid_y = function(y, aterms) {
      if (any(y < 0) || any(y != round(y))) {
        stop("bcm_survey(): the response is a count of returns",
             call. = FALSE)
      }
    },
    init_dpars = list(mu = function(y, aterms) 0.5),
    structure = bcm_survey_structure())
  fam[["survey"]] <- list(nmax = as.integer(nmax))
  fam
}

# The log probability of each possible number sent, as a plain vector at
# plain parameter values or as a list of AD scalars on the tape. One
# function, so the posterior latent_probs() reports and the total the
# objective maximizes cannot disagree.
bcm_survey_terms <- function(y, theta, nmax) {
  nmin <- max(y)
  lapply(nmin:nmax, function(nn) {
    -log(nmax) + sum(RTMB::dbinom(y, nn, theta, log = TRUE))
  })
}

bcm_survey_structure <- function() {
  frmtmb_structure(
    check_spec = function(resp, spec, av) {
      if (length(spec$responses) > 1L || isTRUE(spec$rescor)) {
        stop("bcm_survey() supports univariate models only",
             call. = FALSE)
      }
    },
    frame_block = function(resp, spec, av, mf, y, n) {
      nmax <- resp$family[["survey"]][["nmax"]]
      if (max(y) > nmax) {
        stop("bcm_survey(nmax = ", nmax, "): ", max(y), " surveys came ",
             "back, so at least that many were sent", call. = FALSE)
      }
      list(n = n, nmax = nmax, nmin = as.integer(max(y)))
    },
    loglik = function(y, dpars, aterms, weights, block, extra) {
      # every row shares the unknown number sent, so a row weight has no
      # meaning here and check_spec would have to refuse weights() to
      # give it one; `weights` is 1
      theta <- bcm_mbcast(dpars[["mu"]], block[["n"]])
      acc <- bcm_lse_start()
      for (term in bcm_survey_terms(y, theta, block[["nmax"]])) {
        acc <- bcm_lse_add(acc, term)
      }
      acc
    },
    unit = "the whole set of returns",
    latent_probs = function(fit, block) {
      rs <- single_response(fit, "a bcm_survey() fit")
      dp <- eval_dpars(fit)[[rs$resp_name]]
      y <- fit$frame[["y"]][[rs$resp_name]]
      theta <- rep(as.numeric(dp[["mu"]]), length.out = length(y))
      lp <- vapply(bcm_survey_terms(y, theta, block[["nmax"]]),
                   as.numeric, 0)
      p <- exp(lp - max(lp))
      m <- matrix(p / sum(p), nrow = 1L)
      colnames(m) <- as.character(block[["nmin"]]:block[["nmax"]])
      m
    },
    refusals = list(
      newdata_response = paste0(
        "predict(type = 'response') on a bcm_survey() fit is not ",
        "available for newdata: the expected number of returns depends ",
        "on a number of surveys that only the observed counts identify"),
      conditional_effects = paste0(
        "conditional_effects() is not available for a bcm_survey() fit: ",
        "the number of surveys sent is summed out of the whole response ",
        "at once, and the synthetic grid this function builds is not a ",
        "response"),
      osa = paste0(
        "residuals(type = 'osa') is not available for a bcm_survey() ",
        "fit: the tape holds one number for the whole response with no ",
        "registered observation vector")))
}

## ---- bcm-changepoint ----
# Lee and Wagenmakers chapter 5.4. A sequence of measurements has one
# mean before some moment and another after it, and the moment is not
# known.
#
# The Stan port declares that moment as a CONTINUOUS parameter and
# branches on `t[i] - tau < 0`. Its log density is then piecewise
# constant in tau, so its gradient with respect to tau is zero almost
# everywhere and NUTS has nothing to follow. The book's model puts a
# discrete uniform on the changepoint, and that is what this family
# does: sum over the n - 1 places the change could sit.
bcm_changepoint <- function(time) {
  time_expr <- substitute(time)
  fam <- frmtmb_family(
    "bcm_changepoint",
    dpars = c("mu1", "mu2", "sigma"),
    links = list(mu1 = "identity", mu2 = "identity", sigma = "log"),
    primary_dpars = "mu1",
    type = "continuous",
    lpdf = function(y, dpars, aterms, extra = NULL) {
      stop("A bcm_changepoint() measurement has no row-wise log ",
           "density: which mean it came from depends on where the ",
           "change is, and the change is summed out of the whole ",
           "sequence at once. Use logLik() for the total, or ",
           "latent_probs() for the posterior over the changepoint",
           call. = FALSE)
    },
    init_dpars = list(
      mu1 = function(y, aterms) mean(y[seq_len(max(1, length(y) %/% 4))]),
      mu2 = function(y, aterms) mean(y[seq(length(y) - length(y) %/% 4 + 1,
                                           length(y))]),
      sigma = function(y, aterms) stats::sd(y)),
    structure = bcm_changepoint_structure())
  fam[["cp"]] <- list(time_expr = time_expr)
  fam
}

# The log probability of each changepoint. Accumulated forward and
# backward so that the whole set costs one pass rather than one pass per
# candidate: A is the log density of everything up to tau under the
# first mean, B of everything after it under the second.
bcm_changepoint_terms <- function(y, mu1, mu2, sigma, ord) {
  yo <- y[ord]
  ll1 <- RTMB::dnorm(yo, mu1[ord], sigma[ord], log = TRUE)
  ll2 <- RTMB::dnorm(yo, mu2[ord], sigma[ord], log = TRUE)
  n <- length(yo)
  lpri <- -log(n - 1)
  A <- 0
  B <- sum(ll2)
  out <- vector("list", n - 1L)
  for (tau in seq_len(n - 1L)) {
    A <- A + ll1[tau]
    B <- B - ll2[tau]
    out[[tau]] <- lpri + A + B
  }
  out
}

bcm_changepoint_structure <- function() {
  frmtmb_structure(
    frame_vars = function(fam) list(fam[["cp"]][["time_expr"]]),
    check_spec = function(resp, spec, av) {
      if (length(spec$responses) > 1L || isTRUE(spec$rescor)) {
        stop("bcm_changepoint() supports univariate models only",
             call. = FALSE)
      }
    },
    frame_block = function(resp, spec, av, mf, y, n) {
      tv <- as.numeric(eval(resp$family[["cp"]][["time_expr"]], mf,
                            resp$formula_env))
      if (anyNA(tv)) {
        stop("bcm_changepoint(): the time variable has missing values, ",
             "so the order the change could fall in is undefined",
             call. = FALSE)
      }
      if (n < 3L) {
        stop("bcm_changepoint(): a change needs at least one ",
             "measurement on each side", call. = FALSE)
      }
      list(n = n, ord = order(tv), time = tv)
    },
    loglik = function(y, dpars, aterms, weights, block, extra) {
      n <- block[["n"]]
      acc <- bcm_lse_start()
      terms <- bcm_changepoint_terms(y, bcm_mbcast(dpars[["mu1"]], n),
                                     bcm_mbcast(dpars[["mu2"]], n),
                                     bcm_mbcast(dpars[["sigma"]], n),
                                     block[["ord"]])
      for (term in terms) acc <- bcm_lse_add(acc, term)
      acc
    },
    unit = "the whole sequence",
    latent_probs = function(fit, block) {
      rs <- single_response(fit, "a bcm_changepoint() fit")
      dp <- eval_dpars(fit)[[rs$resp_name]]
      y <- fit$frame[["y"]][[rs$resp_name]]
      n <- block[["n"]]
      g <- function(nm) rep(as.numeric(dp[[nm]]), length.out = n)
      lp <- vapply(bcm_changepoint_terms(y, g("mu1"), g("mu2"),
                                         g("sigma"), block[["ord"]]),
                   as.numeric, 0)
      p <- exp(lp - max(lp))
      m <- matrix(p / sum(p), nrow = 1L)
      colnames(m) <- as.character(seq_len(n - 1L))
      m
    },
    refusals = list(
      newdata_response = paste0(
        "predict(type = 'response') on a bcm_changepoint() fit is not ",
        "available for newdata: which mean a measurement belongs to ",
        "depends on where the change fell in the OBSERVED sequence"),
      conditional_effects = paste0(
        "conditional_effects() is not available for a ",
        "bcm_changepoint() fit: the changepoint is summed out of the ",
        "whole sequence at once"),
      osa = paste0(
        "residuals(type = 'osa') is not available for a ",
        "bcm_changepoint() fit: the tape holds one number for the whole ",
        "sequence with no registered observation vector")))
}

## ---- bcm-two-country ----
# Lee and Wagenmakers chapter 6.4, the two-country quiz. Each person
# comes from one of two countries and so does each question, neither is
# recorded, and a person answers a question from their own country
# correctly with rate alpha and one from the other country with rate
# beta.
#
# The Stan port of this model does not exist. Its file says so: "this
# model is not implemented", with a pointer to a mailing-list thread
# about double mixtures. The obstacle is real and it is the reason the
# sum below is written the way it is. Marginalizing a person's country
# and a question's country independently is wrong, because a cell reads
# BOTH, so the two sets of indicators do not factorize and there is no
# per-row log_mix to write.
#
# They do factorize CONDITIONALLY, which is the whole trick:
#
#   P(k) = sum over question patterns w of  (1/2)^nq
#            prod over people i of  sum over z_i of (1/2) P(k_i | z_i, w)
#
# Given the questions' countries, the people are independent. So the
# cost is 2^nq times the cost of one pass over the data, not
# 2^(np + nq). With eight questions that is 256 passes.
bcm_two_country <- function(person, question, max_patterns = 4096L) {
  person_expr <- substitute(person)
  question_expr <- substitute(question)
  fam <- frmtmb_family(
    "bcm_two_country",
    dpars = c("alpha", "beta"),
    links = list(alpha = "logit", beta = "logit"),
    primary_dpars = "alpha",
    type = "discrete",
    lpdf = function(y, dpars, aterms, extra = NULL) {
      stop("A bcm_two_country() answer has no row-wise log density: ",
           "whether it was an own-country question depends on the ",
           "person's country AND the question's, and neither was ",
           "recorded. Use logLik() for the total, or latent_probs() ",
           "for the posterior over each person's country", call. = FALSE)
    },
    valid_y = function(y, aterms) {
      if (!all(y %in% c(0, 1))) {
        stop("bcm_two_country(): the response is whether the answer was ",
             "correct, coded 0 or 1", call. = FALSE)
      }
    },
    init_dpars = list(alpha = function(y, aterms) 0.8,
                      beta = function(y, aterms) 0.2),
    structure = bcm_two_country_structure(max_patterns))
  fam[["tcq"]] <- list(person_expr = person_expr,
                       question_expr = question_expr)
  fam
}

# The log likelihood, and the per-person posterior country probability
# that falls out of the same sum. `want_probs` picks which one, because
# the two share every term and computing them separately is how they
# drift apart.
bcm_two_country_sum <- function(y, a, b, block, want_probs = FALSE) {
  W <- block[["patterns"]]
  cells <- block[["cells"]]
  np <- block[["np"]]
  nq <- block[["nq"]]
  # log P(y | own-country) and log P(y | other-country), per row
  lown <- y * log(a) + (1 - y) * log(1 - a)
  loth <- y * log(b) + (1 - y) * log(1 - b)
  lpri_w <- -nq * log(2)
  acc <- bcm_lse_start()
  per_person <- if (want_probs) matrix(0, np, 2L) else NULL
  wacc <- if (want_probs) numeric(nrow(W)) else NULL
  pacc <- if (want_probs) array(0, c(nrow(W), np, 2L)) else NULL
  for (p in seq_len(nrow(W))) {
    wp <- W[p, ]
    tot <- lpri_w
    for (i in seq_len(np)) {
      r <- cells[[i]]
      qi <- block[["q"]][r]
      # z = 1: this person's country is 1, so a question with w = 1 is
      # an own-country question
      # 1 where this question's country matches country 1. Arithmetic
      # rather than ifelse(): ifelse() assigns into a copy, which on an
      # advector is both slow and a way to lose the class.
      same1 <- as.numeric(wp[qi] == 1L)
      lz1 <- sum(same1 * lown[r] + (1 - same1) * loth[r])
      lz2 <- sum(same1 * loth[r] + (1 - same1) * lown[r])
      li <- RTMB::logspace_add(lz1 - log(2), lz2 - log(2))
      tot <- tot + li
      if (want_probs) {
        pacc[p, i, 1L] <- as.numeric(lz1 - log(2) - li)
        pacc[p, i, 2L] <- as.numeric(lz2 - log(2) - li)
      }
    }
    acc <- bcm_lse_add(acc, tot)
    if (want_probs) wacc[p] <- as.numeric(tot)
  }
  if (!want_probs) {
    return(acc)
  }
  wt <- exp(wacc - max(wacc))
  wt <- wt / sum(wt)
  for (i in seq_len(np)) {
    per_person[i, 1L] <- sum(wt * exp(pacc[, i, 1L]))
    per_person[i, 2L] <- sum(wt * exp(pacc[, i, 2L]))
  }
  # Which country is called 1 is not identified: every question pattern
  # has a complement under which every person's country flips, and the
  # two carry the same likelihood. So `person` is EXACTLY one half in
  # every cell, and that is the honest report rather than a failure.
  #
  # `agreement` is the quantity that survives the flip: the probability
  # that two people came from the same country, averaged over the
  # patterns. Given a pattern the people are independent, so within one
  # pattern it is just p_i1 p_j1 + p_i2 p_j2.
  agree <- matrix(0, np, np)
  for (p in seq_len(nrow(W))) {
    q1 <- exp(pacc[p, , 1L])
    q2 <- exp(pacc[p, , 2L])
    agree <- agree + wt[p] * (outer(q1, q1) + outer(q2, q2))
  }
  list(person = per_person, agreement = agree)
}

# The probability that each pair of people came from the same country.
# This is what the chapter's figure shows and what a label-switched
# posterior can still answer; see the note in bcm_two_country_sum().
bcm_two_country_agreement <- function(fit) {
  rs <- single_response(fit, "a bcm_two_country() fit")
  block <- frame_block_of(fit$frame, rs$resp_name)
  dp <- eval_dpars(fit)[[rs$resp_name]]
  y <- fit$frame[["y"]][[rs$resp_name]]
  n <- block[["n"]]
  a <- rep(as.numeric(dp[["alpha"]]), length.out = n)
  b <- rep(as.numeric(dp[["beta"]]), length.out = n)
  m <- bcm_two_country_sum(y, a, b, block, want_probs = TRUE)[["agreement"]]
  dimnames(m) <- list(block[["person_levels"]], block[["person_levels"]])
  m
}

bcm_two_country_structure <- function(max_patterns) {
  frmtmb_structure(
    frame_vars = function(fam) {
      list(fam[["tcq"]][["person_expr"]], fam[["tcq"]][["question_expr"]])
    },
    check_spec = function(resp, spec, av) {
      if (length(spec$responses) > 1L || isTRUE(spec$rescor)) {
        stop("bcm_two_country() supports univariate models only",
             call. = FALSE)
      }
    },
    frame_block = function(resp, spec, av, mf, y, n) {
      t_ <- resp$family[["tcq"]]
      pv <- factor(eval(t_[["person_expr"]], mf, resp$formula_env))
      qv <- factor(eval(t_[["question_expr"]], mf, resp$formula_env))
      np <- nlevels(pv)
      nq <- nlevels(qv)
      if (2^nq > max_patterns) {
        stop("bcm_two_country(): the sum runs over 2^", nq, " patterns ",
             "of question countries, which is more than the ",
             max_patterns, " this family will enumerate. The ",
             "marginalization is exact but exponential in the number of ",
             "QUESTIONS, so raise max_patterns deliberately or use a ",
             "sampler", call. = FALSE)
      }
      patterns <- as.matrix(expand.grid(
        rep(list(1:2), nq), KEEP.OUT.ATTRS = FALSE))
      dimnames(patterns) <- NULL
      list(n = n, np = np, nq = nq,
           p = as.integer(pv), q = as.integer(qv),
           cells = split(seq_len(n), pv),
           patterns = patterns,
           person_levels = levels(pv), question_levels = levels(qv))
    },
    loglik = function(y, dpars, aterms, weights, block, extra) {
      n <- block[["n"]]
      a <- bcm_mbcast(dpars[["alpha"]], n)
      b <- bcm_mbcast(dpars[["beta"]], n)
      bcm_two_country_sum(y, a, b, block)
    },
    unit = "the whole quiz",
    latent_probs = function(fit, block) {
      rs <- single_response(fit, "a bcm_two_country() fit")
      dp <- eval_dpars(fit)[[rs$resp_name]]
      y <- fit$frame[["y"]][[rs$resp_name]]
      n <- block[["n"]]
      a <- rep(as.numeric(dp[["alpha"]]), length.out = n)
      b <- rep(as.numeric(dp[["beta"]]), length.out = n)
      m <- bcm_two_country_sum(y, a, b, block,
                               want_probs = TRUE)[["person"]]
      dimnames(m) <- list(block[["person_levels"]],
                          c("country1", "country2"))
      m
    },
    refusals = list(
      newdata_response = paste0(
        "predict(type = 'response') on a bcm_two_country() fit is not ",
        "available for newdata: whether a question is an own-country ",
        "question is decided by the OBSERVED answers of every person"),
      conditional_effects = paste0(
        "conditional_effects() is not available for a ",
        "bcm_two_country() fit: the countries are summed out of the ",
        "whole quiz at once"),
      osa = paste0(
        "residuals(type = 'osa') is not available for a ",
        "bcm_two_country() fit: the tape holds one number for the whole ",
        "quiz with no registered observation vector")))
}

## ---- bcm-planes ----
# Lee and Wagenmakers chapter 5.6, capture-recapture. x planes were
# marked, a later sample of n contained k marked ones, and the size of
# the fleet has a uniform prior on 1..tmax.
#
# This is NOT a family, because the model has NO FREE PARAMETER. Its
# Stan port says so in its own way: it runs with
# algorithm = "Fixed_param", because there is nothing to sample. frm()
# estimates parameters and a model with none is not a fit, so the port
# refuses the model by name and computes the posterior directly, which
# is the whole content of the chapter.
bcm_planes_posterior <- function(x, n, k, tmax) {
  tmin <- x + n - k
  if (tmin > tmax) {
    stop("bcm_planes_posterior(): at least ", tmin, " planes are needed ",
         "to have marked ", x, " and then seen ", n, " with ", k,
         " marked, which is more than tmax = ", tmax, call. = FALSE)
  }
  tt <- tmin:tmax
  # k of the n sampled were among the x marked, out of a fleet of t
  lp <- -log(tmax) + stats::dhyper(k, x, tt - x, n, log = TRUE)
  p <- exp(lp - max(lp))
  list(t = tt, lp = lp, logmarg = max(lp) + log(sum(p)),
       prob = p / sum(p))
}
