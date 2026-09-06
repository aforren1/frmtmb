# Lee and Wagenmakers, Bayesian Cognitive Modeling, chapters 8 and 9:
# comparing rates and comparing means.
#
# Stan programs adapted from the BSD-3 ports at
# github.com/stan-dev/example-models, directories
# Bayesian_Cognitive_Modeling/ModelSelection/{Rates,Means}.
#
# THE BAYES FACTOR IS NOT HERE. Every model in this chapter exists so
# that a Savage-Dickey density ratio can be read off the posterior of
# one parameter at one point, and a maximum likelihood fit has no
# posterior to read. What these tests do is fit the models and pin the
# quantity the ratio is taken over, so that the evidence test can be
# written the moment frmtmb.sample gains the ratio; see
# dev/bcm-findings.md for the handover, and the sibling lane that owns
# frmtmb.sample for the ratio itself.
#
# Each model's `deltaprior` parameter is dropped from the adapted
# programs. It has no data attached and exists only so that the sampler
# draws from the prior for the denominator of the ratio; there is
# nothing for maximum likelihood to say about it.

# ---------------------------------------------------------------------
# Pledgers_1 and Pledgers_2: two rates, unrestricted and order
# restricted.
# ---------------------------------------------------------------------

bcm_pledgers_data <- function() {
  data.frame(s = c(424, 5416), n = c(777, 9072),
             g = factor(c("pledge", "nopledge"),
                        levels = c("pledge", "nopledge")))
}

bcm_pledgers_code <- function() {
  paste(
    "data {",
    "  int<lower=0> n1; int<lower=0> n2;",
    "  int<lower=0> s1; int<lower=0> s2;",
    "}",
    "parameters {",
    "  real<lower=0, upper=1> theta1;",
    "  real<lower=0, upper=1> theta2;",
    "}",
    "model {",
    "  target += beta_lpdf(theta1 | 1, 1);",
    "  target += beta_lpdf(theta2 | 1, 1);",
    "  target += binomial_lpmf(s1 | n1, theta1);",
    "  target += binomial_lpmf(s2 | n2, theta2);",
    "}",
    sep = "\n")
}

test_that("Pledgers_1 estimates two rates and their difference", {
  d <- bcm_pledgers_data()
  fit <- frm(s | trials(n) ~ 0 + g, family = binomial(), data = d)
  th <- plogis(frm_b(fit))
  expect_equal(th, d$s / d$n, tolerance = 1e-6)
  # delta, the quantity the chapter's Bayes factor is about
  expect_lt(th[1] - th[2], 0)
  expect_gt(th[1] - th[2], -0.1)
})

test_that("Pledgers_1 matches its Stan program", {
  skip_unless_stan_identity()
  d <- bcm_pledgers_data()
  fit <- frm(s | trials(n) ~ 0 + g, family = binomial(), data = d)
  stan_lp_check(
    bcm_pledgers_code(),
    data = list(n1 = 777L, n2 = 9072L, s1 = 424L, s2 = 5416L),
    fit = fit,
    pars = function(f) {
      th <- plogis(frm_b(f))
      list(theta1 = th[1], theta2 = th[2])
    },
    const = 0)
})

test_that("Pledgers_2's order restriction does not bind", {
  d <- bcm_pledgers_data()
  fit <- frm(s | trials(n) ~ 0 + g, family = binomial(), data = d)
  th <- plogis(frm_b(fit))
  # Pledgers_2 restricts theta1 < theta2. An inequality constraint is
  # not a term of the likelihood: it either binds at the estimate or
  # leaves it alone. Here the unrestricted maximum already satisfies it,
  # so the restricted maximum IS the unrestricted one and no separate
  # fit exists to run. The Bayes factor between the two models is a
  # different question and belongs to the sampling package.
  expect_lt(th[1], th[2])
})

# ---------------------------------------------------------------------
# Geurts: a hierarchical probit rate model comparing two groups. The
# effect size the chapter tests is delta = alpha / sigma, the group
# difference in units of the between-subject standard deviation.
# ---------------------------------------------------------------------

bcm_geurts_data <- function() {
  errc <- c(15, 10, 61, 11, 60, 44, 63, 70, 57, 11, 67, 21, 89, 12, 63,
            11, 96, 10, 37, 19, 44, 18, 78, 27, 60, 14)
  nc <- c(89, 74, 128, 87, 128, 121, 128, 128, 128, 78, 128, 106, 128,
          83, 128, 100, 128, 73, 128, 86, 128, 86, 128, 100, 128, 79)
  erra <- c(88, 50, 58, 17, 40, 18, 21, 50, 21, 69, 19, 29, 11, 76, 46,
            36, 37, 72, 27, 92, 13, 39, 53, 31, 49, 57, 17, 10, 12, 21,
            39, 43, 49, 17, 39, 13, 68, 24, 21, 27, 48, 54, 41, 75, 38,
            76, 21, 41, 61, 24, 28, 21)
  na <- c(128, 128, 128, 86, 128, 117, 89, 128, 110, 128, 93, 107, 87,
          128, 128, 113, 128, 128, 98, 128, 93, 116, 128, 116, 128, 128,
          93, 86, 86, 96, 128, 128, 128, 86, 128, 78, 128, 111, 100, 95,
          128, 128, 128, 128, 128, 128, 98, 127, 128, 93, 110, 96)
  data.frame(
    k = c(nc - errc, na - erra),
    n = c(nc, na),
    group = factor(rep(c("control", "adhd"), c(length(nc), length(na))),
                   levels = c("control", "adhd")),
    id = factor(seq_len(length(nc) + length(na))))
}

bcm_geurts_code <- function() {
  paste(
    "data {",
    "  int<lower=1> N; int<lower=1> S;",
    "  array[N] int<lower=0> k; array[N] int<lower=1> n;",
    "  vector[N] adhd; array[N] int<lower=1> subj;",
    "}",
    "parameters {",
    "  real mu; real alpha; real<lower=0> sigma; vector[S] phi;",
    "}",
    "model {",
    "  vector[N] eta;",
    "  for (i in 1 : S)",
    "    target += normal_lpdf(phi[i] | 0, sigma);",
    "  for (i in 1 : N)",
    "    eta[i] = mu + alpha * adhd[i] + phi[subj[i]];",
    "  target += binomial_lpmf(k | n, Phi(eta));",
    "}",
    sep = "\n")
}

test_that("Geurts is a hierarchical probit rate comparison", {
  skip_unless_bcm("binomial-extras.R")
  d <- bcm_geurts_data()
  fit <- frm(k | trials(n) ~ group + (1 | id),
             family = bcm_binomial_probit(), data = d)
  b <- fixef(fit)$mu
  sd_id <- sqrt(unname(VarCorr(fit)[[1L]])[1, 1])
  delta <- unname(b["groupadhd"]) / sd_id
  # the chapter's conclusion is that the two groups barely differ, and
  # the standardized effect is what its Bayes factor is taken over
  expect_lt(abs(delta), 0.5)
  expect_gt(sd_id, 0)
  # GeurtsOrderRestricted adds delta < 0, the direction in which the
  # ADHD group is worse. An inequality is not a term of the likelihood:
  # it either binds at the estimate or leaves it alone, and the sign
  # here decides which. Recorded rather than asserted, because the
  # chapter's own point is that the effect is near zero.
  expect_true(is.finite(delta))
})

test_that("Geurts matches its Stan program", {
  skip_unless_stan_identity()
  skip_unless_bcm("binomial-extras.R")
  d <- bcm_geurts_data()
  fit <- frm(k | trials(n) ~ group + (1 | id),
             family = bcm_binomial_probit(), data = d)
  stan_lp_check(
    bcm_geurts_code(),
    data = list(N = nrow(d), S = nlevels(d$id), k = as.integer(d$k),
                n = as.integer(d$n),
                adhd = as.numeric(d$group == "adhd"),
                subj = as.integer(d$id)),
    fit = fit,
    pars = function(f) {
      b <- fixef(f)$mu
      list(mu = unname(b["(Intercept)"]),
           alpha = unname(b["groupadhd"]),
           sigma = frm_sd_term(f, "1 | id"),
           phi = frm_u_term(f, "1 | id"))
    },
    inner = "phi",
    const = 0)
})

# ---------------------------------------------------------------------
# Zeelenberg: the same comparison within subjects.
# ---------------------------------------------------------------------

bcm_zeelenberg_data <- function() {
  sb <- c(15, 11, 15, 14, 15, 18, 16, 16, 18, 16, 15, 13, 18, 12, 11, 13,
          17, 18, 16, 11, 17, 18, 12, 18, 18, 14, 21, 18, 17, 10, 11, 12,
          16, 18, 17, 15, 19, 12, 21, 15, 16, 20, 15, 19, 16, 16, 14, 18,
          16, 19, 17, 11, 19, 18, 16, 16, 11, 19, 18, 12, 15, 18, 20, 8,
          12, 19, 16, 16, 16, 12, 18, 17, 11, 20)
  sn <- c(15, 12, 14, 15, 13, 14, 10, 17, 13, 16, 16, 10, 15, 15, 10, 14,
          17, 18, 19, 12, 19, 18, 10, 18, 16, 13, 15, 20, 13, 15, 13, 14,
          19, 19, 19, 18, 13, 12, 19, 16, 14, 17, 15, 16, 15, 16, 13, 15,
          14, 19, 12, 11, 17, 13, 18, 13, 13, 19, 18, 13, 13, 16, 18, 14,
          14, 17, 12, 12, 16, 14, 16, 18, 13, 13)
  ns <- length(sb)
  data.frame(s = c(sb, sn), n = 21L,
             both = rep(c(1, 0), each = ns),
             id = factor(rep(seq_len(ns), 2)))
}

test_that("Zeelenberg is a within-subject probit rate comparison", {
  skip_unless_bcm("binomial-extras.R")
  d <- bcm_zeelenberg_data()
  fit <- frm(s | trials(n) ~ both + (1 + both || id),
             family = bcm_binomial_probit(), data = d)
  b <- fixef(fit)$mu
  sd_both <- frm_sd_term(fit, "0 + both | id")
  delta <- unname(b["both"]) / sd_both
  # studying both words helps, which is the direction the chapter's
  # order-restricted model assumes
  expect_gt(unname(b["both"]), 0)
  expect_gt(delta, 0)
})

# ---------------------------------------------------------------------
# OneSample and TwoSample: comparing means. Both are a plain normal
# model in a reparameterization that names the standardized effect
# directly, and maximum likelihood reads the same effect off the fit.
# ---------------------------------------------------------------------

bcm_onesample_data <- function() {
  winter <- c(-0.05, 0.41, 0.17, -0.13, 0.00, -0.05, 0.00, 0.17, 0.29,
              0.04, 0.21, 0.08, 0.37, 0.17, 0.08, -0.04, -0.04, 0.04,
              -0.13, -0.12, 0.04, 0.21, 0.17, 0.17, 0.17, 0.33, 0.04,
              0.04, 0.04, 0.00, 0.21, 0.13, 0.25, -0.05, 0.29, 0.42,
              -0.05, 0.12, 0.04, 0.25, 0.12)
  summer <- c(0.00, 0.38, -0.12, 0.12, 0.25, 0.12, 0.13, 0.37, 0.00,
              0.50, 0.00, 0.00, -0.13, -0.37, -0.25, -0.12, 0.50, 0.25,
              0.13, 0.25, 0.25, 0.38, 0.25, 0.12, 0.00, 0.00, 0.00,
              0.00, 0.25, 0.13, -0.25, -0.38, -0.13, -0.25, 0.00, 0.00,
              -0.12, 0.25, 0.00, 0.50, 0.00)
  x <- winter - summer
  data.frame(x = x / stats::sd(x))
}

test_that("OneSample's delta is the standardized mean", {
  d <- bcm_onesample_data()
  fit <- frm(x ~ 1, data = d)
  mu <- unname(fixef(fit)$mu)
  sigma <- unname(exp(fixef(fit)$sigma))
  # the book writes mu = delta * sigma, which is a reparameterization
  # and not a different model, so delta reads straight off the fit
  delta <- mu / sigma
  expect_equal(mu, mean(d$x), tolerance = 1e-6)
  expect_lt(abs(delta), 0.5)
  # OneSampleOrderRestricted adds delta > 0. The unrestricted estimate
  # decides whether that constraint binds, and here it does not, so the
  # restricted maximum is the unrestricted one.
  expect_gt(delta, 0)
})

bcm_twosample_data <- function() {
  x <- c(70, 80, 79, 83, 77, 75, 84, 78, 75, 75, 78, 82, 74, 81, 72, 70,
         75, 72, 76, 77)
  y <- c(56, 80, 63, 62, 67, 71, 68, 76, 79, 67, 76, 74, 67, 70, 62, 65,
         72, 72, 69, 71)
  s <- stats::sd(x)
  m <- mean(x)
  data.frame(v = c((x - m) / s, (y - m) / s),
             g = factor(rep(c("x", "y"), c(length(x), length(y))),
                        levels = c("x", "y")))
}

bcm_twosample_code <- function() {
  paste(
    "data {",
    "  int<lower=1> n1; int<lower=1> n2;",
    "  vector[n1] x; vector[n2] y;",
    "}",
    "parameters {",
    "  real mu; real<lower=0> sigma; real delta;",
    "}",
    "model {",
    "  real alpha = delta * sigma;",
    "  target += normal_lpdf(x | mu + alpha / 2, sigma);",
    "  target += normal_lpdf(y | mu - alpha / 2, sigma);",
    "}",
    sep = "\n")
}

test_that("TwoSample's delta is the standardized group difference", {
  d <- bcm_twosample_data()
  fit <- frm(v ~ g, data = d)
  b <- fixef(fit)$mu
  sigma <- unname(exp(fixef(fit)$sigma))
  delta <- -unname(b["gy"]) / sigma
  # the first group scores higher, by more than one standard deviation
  expect_gt(delta, 1)
})

test_that("TwoSample matches its Stan program", {
  skip_unless_stan_identity()
  d <- bcm_twosample_data()
  fit <- frm(v ~ g, data = d)
  stan_lp_check(
    bcm_twosample_code(),
    data = list(n1 = 20L, n2 = 20L,
                x = d$v[d$g == "x"], y = d$v[d$g == "y"]),
    fit = fit,
    pars = function(f) {
      b <- fixef(f)$mu
      sigma <- unname(exp(fixef(f)$sigma))
      # Stan's mu is the grand mean and alpha the difference; frmtmb's
      # intercept is the first group's mean and `gy` the difference the
      # other way round
      list(mu = unname(b["(Intercept)"]) + unname(b["gy"]) / 2,
           sigma = sigma,
           delta = -unname(b["gy"]) / sigma)
    },
    # the original's cauchy(0, 1) priors on mu, sigma and delta are
    # dropped: they are the Bayes factor's prior, not the likelihood
    const = 0)
})
