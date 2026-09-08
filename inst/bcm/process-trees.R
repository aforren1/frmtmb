# Multinomial process trees: a small number of latent rates, a fixed
# tree that turns them into the probabilities of a few observable
# categories, and a multinomial count in each category.
#
# vignette("bayesian-cognitive-modeling") shows these in place;
# tests/testthat/test-bcm-mpt.R and
# tests/testthat/test-bcm-data-analysis.R fit this code. Load it with:
#
#   source(system.file("bcm", "process-trees.R", package = "frmtmb"))
#
# Both families are ROWWISE: one row is one respondent's category
# counts, and its probability reads nothing else. They need a custom
# family because the map from rates to category probabilities is the
# model, and no link function has that shape; they do not need
# frmtmb_structure().
#
# The response is an n-by-K count MATRIX, which is the shape the core's
# own multinomial() family takes, so `cbind(k1, k2, k3, k4) ~ ...` is
# the spelling.

## ---- bcm-tree-lpdf ----
# The multinomial log density from LOG category probabilities.
#
# Log probabilities rather than probabilities, because a tree can put a
# category's probability at 1e-300 for a respondent who never produced
# it, and `0 * log(0)` on the tape is NaN rather than the zero the
# arithmetic wants. Every branch below is assembled with logspace_add()
# for the same reason.
bcm_tree_lpdf <- function(y, logtheta) {
  n <- rowSums(y)
  ll <- lgamma(n + 1) - rowSums(lgamma(y + 1))
  for (j in seq_along(logtheta)) {
    ll <- ll + y[, j] * logtheta[[j]]
  }
  ll
}

bcm_tree_valid_y <- function(name, K) {
  function(y, aterms) {
    if (!is.matrix(y) || ncol(y) != K) {
      stop(name, "(): the response must be an n by ", K,
           " matrix of category counts, spelled cbind(k1, ..., k", K,
           ")", call. = FALSE)
    }
    if (any(y < 0) || any(y != round(y))) {
      stop(name, "(): category counts must be non-negative integers",
           call. = FALSE)
    }
  }
}

## ---- bcm-mpt ----
# The pair-clustering tree of Riefer and Batchelder, as Lee and
# Wagenmakers chapter 14 states it. A word pair is either clustered
# (probability c) or not; a clustered pair is recalled together with
# probability r; an unclustered word is recalled on its own with
# probability u. The four observable categories are
#
#   1  both words recalled, adjacently   c r
#   2  both recalled, not adjacently     (1 - c) u^2
#   3  one word recalled                 (1 - c) 2 u (1 - u)
#   4  neither recalled                  c (1 - r) + (1 - c) (1 - u)^2
#
# Category 4 is the one that makes this a TREE rather than a
# reparameterized multinomial: two different paths end there, so its
# probability is a sum and the map is not invertible category by
# category.
#
# The links are probits, because that is the scale the book's
# hierarchical versions put the individual differences on. With
# `1 + (1 | p | id)` on each of the three, the correlation of the
# random-effect block IS the latent-trait correlation the case study
# reports.
bcm_mpt_pairs <- function(link = "probit") {
  frmtmb_family(
    "bcm_mpt_pairs",
    dpars = c("c", "r", "u"),
    links = list(c = link, r = link, u = link),
    primary_dpars = "c",
    type = "discrete",
    lpdf = function(y, dpars, aterms) {
      lc <- log(dpars[["c"]])
      l1c <- log(1 - dpars[["c"]])
      lr <- log(dpars[["r"]])
      l1r <- log(1 - dpars[["r"]])
      lu <- log(dpars[["u"]])
      l1u <- log(1 - dpars[["u"]])
      bcm_tree_lpdf(y, list(
        lc + lr,
        l1c + 2 * lu,
        l1c + log(2) + lu + l1u,
        RTMB::logspace_add(lc + l1r, l1c + 2 * l1u)))
    },
    valid_y = bcm_tree_valid_y("bcm_mpt_pairs", 4L),
    init_dpars = list(c = function(y, aterms) 0.5,
                      r = function(y, aterms) 0.5,
                      u = function(y, aterms) 0.5),
    sim = function(dpars, aterms, n) {
      cc <- rep(as.numeric(dpars[["c"]]), length.out = n)
      rr <- rep(as.numeric(dpars[["r"]]), length.out = n)
      uu <- rep(as.numeric(dpars[["u"]]), length.out = n)
      size <- aterms[["trials"]]
      if (is.null(size)) {
        stop("simulate(): a bcm_mpt_pairs() fit needs trials() to know ",
             "how many word pairs each respondent saw", call. = FALSE)
      }
      size <- rep(as.numeric(size), length.out = n)
      out <- matrix(0, n, 4L)
      for (i in seq_len(n)) {
        p <- c(cc[i] * rr[i],
               (1 - cc[i]) * uu[i]^2,
               (1 - cc[i]) * 2 * uu[i] * (1 - uu[i]),
               cc[i] * (1 - rr[i]) + (1 - cc[i]) * (1 - uu[i])^2)
        out[i, ] <- stats::rmultinom(1L, size[i], p)
      }
      out
    })
}

# The three rates at the estimates, which is what the case study plots.
# `dpar_linpred()` would give the linear predictors; this applies the
# link so the numbers are on the book's scale.
bcm_mpt_rates <- function(fit) {
  rs <- single_response(fit, "a bcm_mpt_pairs() fit")
  dp <- eval_dpars(fit)[[rs$resp_name]]
  vapply(c("c", "r", "u"), function(nm) as.numeric(dp[[nm]])[1L], 0)
}

## ---- bcm-kappa ----
# The kappa coefficient of agreement, Lee and Wagenmakers chapter 5.3.
# Two methods each say one or zero. The four counts are
#
#   1  both say one            alpha beta
#   2  objective one, other zero    alpha (1 - beta)
#   3  objective zero, other one    (1 - alpha) (1 - gamma)
#   4  both say zero           (1 - alpha) gamma
#
# where alpha is the rate at which the objective method says one, beta
# the surrogate's hit rate and gamma its correct-rejection rate. Unlike
# the pair-clustering tree this map IS invertible, so a saturated fit
# reproduces the four cell proportions exactly and the interest is in
# the derived quantities rather than in the fit.
bcm_kappa <- function(link = "logit") {
  frmtmb_family(
    "bcm_kappa",
    dpars = c("alpha", "beta", "gamma"),
    links = list(alpha = link, beta = link, gamma = link),
    primary_dpars = "alpha",
    type = "discrete",
    lpdf = function(y, dpars, aterms) {
      la <- log(dpars[["alpha"]])
      l1a <- log(1 - dpars[["alpha"]])
      lb <- log(dpars[["beta"]])
      l1b <- log(1 - dpars[["beta"]])
      lg <- log(dpars[["gamma"]])
      l1g <- log(1 - dpars[["gamma"]])
      bcm_tree_lpdf(y, list(la + lb, la + l1b, l1a + l1g, l1a + lg))
    },
    valid_y = bcm_tree_valid_y("bcm_kappa", 4L),
    init_dpars = list(alpha = function(y, aterms) 0.5,
                      beta = function(y, aterms) 0.5,
                      gamma = function(y, aterms) 0.5),
    sim = function(dpars, aterms, n) {
      a <- rep(as.numeric(dpars[["alpha"]]), length.out = n)
      b <- rep(as.numeric(dpars[["beta"]]), length.out = n)
      g <- rep(as.numeric(dpars[["gamma"]]), length.out = n)
      size <- aterms[["trials"]]
      if (is.null(size)) {
        stop("simulate(): a bcm_kappa() fit needs trials() to know how ",
             "many cases were rated", call. = FALSE)
      }
      size <- rep(as.numeric(size), length.out = n)
      out <- matrix(0, n, 4L)
      for (i in seq_len(n)) {
        p <- c(a[i] * b[i], a[i] * (1 - b[i]),
               (1 - a[i]) * (1 - g[i]), (1 - a[i]) * g[i])
        out[i, ] <- stats::rmultinom(1L, size[i], p)
      }
      out
    })
}

# Chance-corrected agreement at the estimates: the book's xi, psi and
# kappa, computed from the three rates rather than from the counts, so
# that a model with covariates on the rates reports one kappa per row.
bcm_kappa_summary <- function(fit) {
  rn <- single_response(fit, "a bcm_kappa() fit")$resp_name
  dp <- eval_dpars(fit)[[rn]]
  a <- as.numeric(dp[["alpha"]])
  b <- as.numeric(dp[["beta"]])
  g <- as.numeric(dp[["gamma"]])
  pi1 <- a * b
  pi2 <- a * (1 - b)
  pi3 <- (1 - a) * (1 - g)
  pi4 <- (1 - a) * g
  xi <- a * b + (1 - a) * g
  psi <- (pi1 + pi2) * (pi1 + pi3) + (pi2 + pi4) * (pi3 + pi4)
  list(alpha = a, beta = b, gamma = g, xi = xi, psi = psi,
       kappa = (xi - psi) / (1 - psi))
}
