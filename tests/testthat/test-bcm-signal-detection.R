# Lee and Wagenmakers, Bayesian Cognitive Modeling, chapter 11: signal
# detection theory.
#
# Stan programs adapted from the BSD-3 ports at
# github.com/stan-dev/example-models, directory
# Bayesian_Cognitive_Modeling/CaseStudies/SignalDetection.
#
# The frmtmb spelling is a BINOMIAL PROBIT GLM. Equal-variance Gaussian
# SDT says the hit rate is Phi(d/2 - c) and the false-alarm rate is
# Phi(-d/2 - c), so writing one row per (case, trial type) with
#
#   half = +1/2 on signal trials and -1/2 on noise trials
#   bias = -1 everywhere
#
# makes the linear predictor d * half + c * bias, and the two regression
# coefficients ARE discriminability and criterion, on the book's signs.
# The alternative coding, an intercept plus a signal indicator, is the
# more familiar GLM one but puts the random effects on the wrong pair:
# it makes the criterion and d' correlated by construction, and the
# book's hierarchical model gives them independent group distributions.

bcm_sdt_long <- function(h, f, s, n) {
  k <- length(h)
  data.frame(
    y = c(h, f),
    N = c(s, n),
    case = factor(rep(seq_len(k), 2)),
    half = rep(c(0.5, -0.5), each = k),
    bias = -1)
}

# SDT_1: one discriminability and one criterion per case.
bcm_sdt1_code <- function() {
  paste(
    "data {",
    "  int<lower=1> k;",
    "  array[k] int<lower=0> h; array[k] int<lower=0> f;",
    "  array[k] int<lower=0> s; array[k] int<lower=0> n;",
    "}",
    "parameters {",
    "  vector[k] d;",
    "  vector[k] c;",
    "}",
    "model {",
    "  vector[k] thetah;",
    "  vector[k] thetaf;",
    "  for (i in 1 : k) {",
    "    thetah[i] = Phi(d[i] / 2 - c[i]);",
    "    thetaf[i] = Phi(-d[i] / 2 - c[i]);",
    "  }",
    "  target += normal_lpdf(d | 0, inv_sqrt(0.5));",
    "  target += normal_lpdf(c | 0, inv_sqrt(2));",
    "  target += binomial_lpmf(h | s, thetah);",
    "  target += binomial_lpmf(f | n, thetaf);",
    "}",
    sep = "\n")
}

# SDT_2: discriminability and criterion drawn from independent group
# distributions. The vague regularizers of the original, normal(0,
# 1 / sqrt(.001)) on the group means and gamma(.001, .001) on the group
# precisions, are dropped: they exist so the sampler has a proper
# target, and dropping them on both sides makes this a likelihood
# identity.
bcm_sdt2_code <- function() {
  paste(
    "data {",
    "  int<lower=1> k;",
    "  array[k] int<lower=0> h; array[k] int<lower=0> f;",
    "  int<lower=0> s; int<lower=0> n;",
    "}",
    "parameters {",
    "  vector[k] d;",
    "  vector[k] c;",
    "  real muc; real mud;",
    "  real<lower=0> sigmac; real<lower=0> sigmad;",
    "}",
    "model {",
    "  vector[k] thetah;",
    "  vector[k] thetaf;",
    "  for (i in 1 : k) {",
    "    thetah[i] = Phi(d[i] / 2 - c[i]);",
    "    thetaf[i] = Phi(-d[i] / 2 - c[i]);",
    "  }",
    "  target += normal_lpdf(c | muc, sigmac);",
    "  target += normal_lpdf(d | mud, sigmad);",
    "  target += binomial_lpmf(h | s, thetah);",
    "  target += binomial_lpmf(f | n, thetaf);",
    "}",
    sep = "\n")
}

# SDT_3: the same model written with parameter expansion. `xi` is a
# redundant multiplicative parameter that a Gibbs sampler uses to move
# faster; the likelihood is SDT_2's, so under maximum likelihood xi and
# the raw scale are a ridge and only their product is identified. The
# map below therefore fixes xi at 1, which is the point of the ridge
# where the two programs are the same function.
bcm_sdt3_code <- function() {
  paste(
    "data {",
    "  int<lower=1> k;",
    "  array[k] int<lower=0> h; array[k] int<lower=0> f;",
    "  int<lower=0> s; int<lower=0> n;",
    "}",
    "parameters {",
    "  real muc; real mud;",
    "  real<lower=0> sigmacnew; real<lower=0> sigmadnew;",
    "  real<lower=0, upper=1> xic; real<lower=0, upper=1> xid;",
    "  vector[k] deltac; vector[k] deltad;",
    "}",
    "model {",
    "  vector[k] c = muc + xic * deltac;",
    "  vector[k] d = mud + xid * deltad;",
    "  vector[k] thetah;",
    "  vector[k] thetaf;",
    "  for (i in 1 : k) {",
    "    thetah[i] = Phi(d[i] / 2 - c[i]);",
    "    thetaf[i] = Phi(-d[i] / 2 - c[i]);",
    "  }",
    "  target += normal_lpdf(deltac | 0, sigmacnew);",
    "  target += normal_lpdf(deltad | 0, sigmadnew);",
    "  target += binomial_lpmf(h | s, thetah);",
    "  target += binomial_lpmf(f | n, thetaf);",
    "}",
    sep = "\n")
}

# ---------------------------------------------------------------------
# SDT_1: three independent cases.
#
# The book's demo data set is not used here. Its third case is
# h = 10 / 10 and f = 0 / 10, which separates completely, so the maximum
# likelihood estimate of d is infinite; the Bayesian fit is finite only
# because of the prior on d. The Lehrner et al. (1995) data of the same
# script have no separated case, and the prior below is the one the
# original program uses, carried on both sides.
# ---------------------------------------------------------------------

bcm_sdt1_counts <- function() {
  m <- matrix(c(148, 29, 32, 151,
                150, 40, 30, 140,
                150, 51, 40, 139), nrow = 3, ncol = 4, byrow = TRUE)
  list(h = m[, 1], f = m[, 2], s = m[, 1] + m[, 3], n = m[, 2] + m[, 4])
}

bcm_sdt1_prior <- function(k = 3) {
  Reduce(`+`, c(
    lapply(seq_len(k), function(i) {
      set_prior(sprintf("normal(0, %.15g)", sqrt(2)), class = "b",
                coef = paste0("case", i, ":half"))
    }),
    lapply(seq_len(k), function(i) {
      set_prior(sprintf("normal(0, %.15g)", 1 / sqrt(2)), class = "b",
                coef = paste0("case", i, ":bias"))
    })))
}

test_that("SDT_1 is a probit GLM whose coefficients are d and c", {
  skip_unless_bcm("binomial-extras.R")
  cn <- bcm_sdt1_counts()
  d <- bcm_sdt_long(cn$h, cn$f, cn$s, cn$n)
  fit <- frm(y | trials(N) ~ 0 + case:half + case:bias,
             family = bcm_binomial_probit(), data = d,
             prior = bcm_sdt1_prior())
  b <- fixef(fit)$mu
  dd <- unname(b[paste0("case", 1:3, ":half")])
  cc <- unname(b[paste0("case", 1:3, ":bias")])
  # the fitted hit and false-alarm rates are the book's transforms
  expect_equal(pnorm(dd / 2 - cc), unname(fitted(fit)[1:3] / d$N[1:3]),
               tolerance = 0.02)
  expect_equal(pnorm(-dd / 2 - cc), unname(fitted(fit)[4:6] / d$N[4:6]),
               tolerance = 0.02)
  # all three cases discriminate
  expect_true(all(dd > 1))
})

test_that("SDT_1 matches its Stan program", {
  skip_unless_stan_identity()
  skip_unless_bcm("binomial-extras.R")
  cn <- bcm_sdt1_counts()
  d <- bcm_sdt_long(cn$h, cn$f, cn$s, cn$n)
  fit <- frm(y | trials(N) ~ 0 + case:half + case:bias,
             family = bcm_binomial_probit(), data = d,
             prior = bcm_sdt1_prior())
  stan_lp_check(
    bcm_sdt1_code(),
    data = list(k = 3L, h = cn$h, f = cn$f, s = cn$s, n = cn$n),
    fit = fit,
    pars = function(f) {
      b <- fixef(f)$mu
      list(d = unname(b[paste0("case", 1:3, ":half")]),
           c = unname(b[paste0("case", 1:3, ":bias")]))
    },
    # both programs carry the book's normal priors on d and c
    const = 0)
})

# ---------------------------------------------------------------------
# SDT_2 and SDT_3: 40 subjects, 4 signal and 4 noise trials each.
#
# The induction condition of Heit and Rotello (2005), as the case study
# distributes it. Most subjects sit at a boundary (4 hits out of 4, 0
# false alarms out of 4), so an unpooled fit has no finite estimate; the
# group distribution is what makes every subject's pair finite, and it
# is a random effect here.
# ---------------------------------------------------------------------

bcm_sdt2_counts <- function() {
  m <- matrix(c(
    3, 1, 1, 3, 4, 0, 0, 4, 4, 4, 0, 0, 4, 1, 0, 3, 4, 3, 0, 1,
    4, 4, 0, 0, 4, 1, 0, 3, 4, 3, 0, 1, 4, 0, 0, 4, 4, 2, 0, 2,
    3, 2, 1, 2, 4, 2, 0, 2, 4, 0, 0, 4, 3, 4, 1, 0, 3, 4, 1, 0,
    4, 1, 0, 3, 3, 2, 1, 2, 4, 1, 0, 3, 3, 1, 1, 3, 4, 1, 0, 3,
    3, 2, 1, 2, 4, 0, 0, 4, 3, 4, 1, 0, 4, 4, 0, 0, 4, 2, 0, 2,
    4, 4, 0, 0, 4, 4, 0, 0, 4, 4, 0, 0, 4, 4, 0, 0, 4, 0, 0, 4,
    2, 3, 2, 1, 4, 0, 0, 4, 4, 4, 0, 0, 2, 3, 2, 1, 4, 0, 0, 4,
    4, 3, 0, 1, 4, 4, 0, 0, 4, 1, 0, 3, 4, 1, 0, 3, 4, 4, 0, 0),
    nrow = 40, ncol = 4, byrow = TRUE)
  list(h = m[, 1], f = m[, 2], s = 4L, n = 4L, k = nrow(m))
}

bcm_sdt2_data <- function() {
  cn <- bcm_sdt2_counts()
  d <- bcm_sdt_long(cn$h, cn$f, rep(cn$s, cn$k), rep(cn$n, cn$k))
  names(d)[names(d) == "case"] <- "id"
  d
}

test_that("SDT_2 is a probit GLMM with independent d and c effects", {
  skip_unless_bcm("binomial-extras.R")
  d <- bcm_sdt2_data()
  fit <- frm(y | trials(N) ~ 0 + half + bias + (0 + half + bias || id),
             family = bcm_binomial_probit(), data = d)
  b <- fixef(fit)$mu
  expect_gt(unname(b["half"]), 0)
  # every subject's pair is finite, which is what the group
  # distribution buys: the raw per-subject estimates are not
  expect_true(all(is.finite(as.matrix(ranef(fit)[[1L]]))))
  expect_equal(nrow(ranef(fit)[[1L]]), 40L)
})

test_that("SDT_2 matches its Stan program", {
  skip_unless_stan_identity()
  skip_unless_bcm("binomial-extras.R")
  cn <- bcm_sdt2_counts()
  d <- bcm_sdt2_data()
  fit <- frm(y | trials(N) ~ 0 + half + bias + (0 + half + bias || id),
             family = bcm_binomial_probit(), data = d)
  stan_lp_check(
    bcm_sdt2_code(),
    data = list(k = cn$k, h = cn$h, f = cn$f, s = cn$s, n = cn$n),
    fit = fit,
    pars = function(f) {
      b <- fixef(f)$mu
      list(d = unname(b["half"]) + frm_u_term(f, "0 + half | id"),
           c = unname(b["bias"]) + frm_u_term(f, "0 + bias | id"),
           mud = unname(b["half"]), muc = unname(b["bias"]),
           sigmad = frm_sd_term(f, "0 + half | id"),
           sigmac = frm_sd_term(f, "0 + bias | id"))
    },
    # the outer estimate maximizes the Laplace-approximated marginal, so
    # only the inner block is at a stationary point of the joint
    inner = c("d", "c"),
    const = 0)
})

test_that("SDT_3's parameter expansion is a ridge, not a model", {
  skip_unless_stan_identity()
  skip_unless_bcm("binomial-extras.R")
  cn <- bcm_sdt2_counts()
  d <- bcm_sdt2_data()
  fit <- frm(y | trials(N) ~ 0 + half + bias + (0 + half + bias || id),
             family = bcm_binomial_probit(), data = d)
  stan_lp_check(
    bcm_sdt3_code(),
    data = list(k = cn$k, h = cn$h, f = cn$f, s = cn$s, n = cn$n),
    fit = fit,
    pars = function(f) {
      b <- fixef(f)$mu
      list(mud = unname(b["half"]), muc = unname(b["bias"]),
           sigmadnew = frm_sd_term(f, "0 + half | id"),
           sigmacnew = frm_sd_term(f, "0 + bias | id"),
           # the point of the ridge where SDT_3 IS SDT_2
           xid = 1, xic = 1,
           deltad = frm_u_term(f, "0 + half | id"),
           deltac = frm_u_term(f, "0 + bias | id"))
    },
    inner = c("deltad", "deltac"),
    const = 0)
})
