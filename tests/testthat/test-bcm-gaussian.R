# Lee and Wagenmakers, Bayesian Cognitive Modeling, chapter 4:
# inferring means and standard deviations.
#
# Stan programs adapted from the BSD-3 ports at
# github.com/stan-dev/example-models, directory
# Bayesian_Cognitive_Modeling/ParameterEstimation/Gaussian. Sampling
# statements become `target +=` (see helper-stan.R); the `<lower, upper>`
# constraints of the originals are kept and written out as the
# `uniform_lpdf` terms they imply, so that the stated constant is
# visible in the program rather than hidden in a declaration.

bcm_gaussian_code <- function() {
  paste(
    "data {",
    "  int<lower=1> n;",
    "  vector[n] x;",
    "}",
    "parameters {",
    "  real mu;",
    "  real<lower=0, upper=10> sigma;",
    "}",
    "model {",
    "  target += normal_lpdf(mu | 0, sqrt(1000));",
    "  target += uniform_lpdf(sigma | 0, 10);",
    "  target += normal_lpdf(x | mu, sigma);",
    "}",
    sep = "\n")
}

# The Seven Scientists, with the book's gamma prior on the precisions
# replaced by a normal prior on the log standard deviations. The note in
# the test below says why: the gamma prior is what makes the BAYESIAN
# model proper, and maximum likelihood needs a penalty of its own.
bcm_seven_code <- function() {
  paste(
    "data {",
    "  int<lower=1> n;",
    "  vector[n] x;",
    "  real<lower=0> tau;",
    "}",
    "parameters {",
    "  real mu;",
    "  vector[n] logsigma;",
    "}",
    "model {",
    "  target += normal_lpdf(logsigma | 0, tau);",
    "  target += normal_lpdf(x | mu, exp(logsigma));",
    "}",
    sep = "\n")
}

bcm_iq_code <- function() {
  paste(
    "data {",
    "  int<lower=1> n; int<lower=1> m;",
    "  matrix[n, m] x;",
    "}",
    "parameters {",
    "  vector<lower=0, upper=300>[n] mu;",
    "  real<lower=0, upper=100> sigma;",
    "}",
    "model {",
    "  target += uniform_lpdf(mu | 0, 300);",
    "  target += uniform_lpdf(sigma | 0, 100);",
    "  for (i in 1:n)",
    "    target += normal_lpdf(x[i] | mu[i], sigma);",
    "}",
    sep = "\n")
}

bcm_gaussian_data <- function() {
  data.frame(x = c(1.1, 1.9, 2.3, 1.8))
}

# ---------------------------------------------------------------------
# Gaussian: one mean, one standard deviation.
# ---------------------------------------------------------------------

test_that("Gaussian is an intercept-only normal", {
  d <- bcm_gaussian_data()
  fit <- frm(x ~ 1, data = d)
  # frmtmb maximizes the likelihood, so sigma is the ML (divide by n)
  # estimate rather than the sample standard deviation
  expect_equal(unname(fixef(fit)$mu), mean(d$x), tolerance = 1e-6)
  expect_equal(unname(exp(fixef(fit)$sigma)),
               sqrt(mean((d$x - mean(d$x))^2)), tolerance = 1e-6)
})

test_that("Gaussian matches its Stan program", {
  skip_unless_stan_identity()
  d <- bcm_gaussian_data()
  # The book's prior on mu is carried through set_prior(), so it is the
  # same penalty on both sides and contributes nothing to the constant.
  fit <- frm(x ~ 1, data = d,
             prior = set_prior(sprintf("normal(0, %.15g)", sqrt(1000)),
                               class = "Intercept"))
  stan_lp_check(
    bcm_gaussian_code(),
    data = list(n = 4L, x = d$x),
    fit = fit,
    pars = function(f) list(mu = unname(fixef(f)$mu),
                            sigma = unname(exp(fixef(f)$sigma))),
    # the only prior frmtmb does not carry is the uniform on sigma,
    # whose log density is -log(10) everywhere inside the interval
    const = -log(10))
})

# ---------------------------------------------------------------------
# Seven Scientists: one mean, one standard deviation per measurement.
#
# The likelihood is UNBOUNDED. Profiling out the seven standard
# deviations leaves -sum(log|x_i - mu|) - 7/2, which goes to +Infinity
# as mu approaches any observation, so there is no interior maximum and
# no local one either: between two observations the profile has a
# minimum, not a maximum. The book's gamma(.001, .001) prior on the
# precisions is what makes the Bayesian model proper; maximum likelihood
# needs a penalty of its own, so the fit below is a MAP fit and the Stan
# program carries the same penalty.
# ---------------------------------------------------------------------

bcm_seven_data <- function() {
  data.frame(x = c(-27.020, 3.570, 8.191, 9.898, 9.603, 9.945, 10.056),
             obs = factor(1:7))
}

test_that("Seven Scientists has an unbounded likelihood", {
  d <- bcm_seven_data()
  # The profile likelihood at sigma_i = |x_i - mu|, evaluated along a
  # sequence of mu approaching the first observation. Nothing is fitted:
  # this is the arithmetic that says a fit cannot exist.
  prof <- function(mu) {
    r <- abs(d$x - mu)
    sum(stats::dnorm(d$x, mu, r, log = TRUE))
  }
  eps <- 10^-(2:12)
  ll <- vapply(d$x[1] + eps, prof, 0)
  expect_true(all(diff(ll) > 0))
  # the divergence is logarithmic: one more decade of eps buys log(10)
  # more log likelihood, forever, so there is no supremum to converge to
  expect_equal(mean(diff(ll)), log(10), tolerance = 1e-6)
  expect_gt(ll[length(ll)] - ll[1], 20)
})

test_that("Seven Scientists fits as a MAP model", {
  d <- bcm_seven_data()
  fit <- frm(bf(x ~ 1, sigma ~ 0 + obs), data = d,
             prior = set_prior("normal(0, 1)", class = "b", dpar = "sigma"))
  # the book's point: the first scientist is discounted, so the estimate
  # sits with the other six rather than at the mean of all seven
  expect_gt(unname(fixef(fit)$mu), 5)
  expect_lt(unname(fixef(fit)$mu), 11)
  expect_lt(unname(fixef(fit)$mu), mean(d$x) + 20)
  # and that scientist's own standard deviation is the largest
  s <- exp(unname(fixef(fit)$sigma))
  expect_equal(which.max(s), 1L)
})

test_that("Seven Scientists matches its Stan program", {
  skip_unless_stan_identity()
  d <- bcm_seven_data()
  fit <- frm(bf(x ~ 1, sigma ~ 0 + obs), data = d,
             prior = set_prior("normal(0, 1)", class = "b", dpar = "sigma"))
  stan_lp_check(
    bcm_seven_code(),
    data = list(n = 7L, x = d$x, tau = 1),
    fit = fit,
    pars = function(f) list(mu = unname(fixef(f)$mu),
                            logsigma = unname(fixef(f)$sigma)),
    # both programs carry the same normal(0, 1) on the log standard
    # deviations, so there is nothing left over
    const = 0)
})

# ---------------------------------------------------------------------
# IQ: three people, three repeated measurements each, one common
# standard deviation.
# ---------------------------------------------------------------------

bcm_iq_matrix <- function() {
  matrix(c(90, 95, 100, 105, 110, 115, 150, 155, 160),
         nrow = 3, ncol = 3, byrow = TRUE)
}

bcm_iq_data <- function() {
  x <- bcm_iq_matrix()
  data.frame(y = as.vector(t(x)),
             person = factor(rep(seq_len(nrow(x)), each = ncol(x))))
}

test_that("IQ is one mean per person with a common sigma", {
  d <- bcm_iq_data()
  fit <- frm(y ~ 0 + person, data = d)
  expect_equal(unname(fixef(fit)$mu), c(95, 110, 155), tolerance = 1e-6)
  expect_equal(unname(exp(fixef(fit)$sigma)),
               sqrt(mean((d$y - rep(c(95, 110, 155), each = 3))^2)),
               tolerance = 1e-6)
})

test_that("IQ matches its Stan program", {
  skip_unless_stan_identity()
  d <- bcm_iq_data()
  fit <- frm(y ~ 0 + person, data = d)
  stan_lp_check(
    bcm_iq_code(),
    data = list(n = 3L, m = 3L, x = bcm_iq_matrix()),
    fit = fit,
    pars = function(f) list(mu = frm_b(f),
                            sigma = unname(exp(fixef(f)$sigma))),
    # Stan's `<lower, upper>` declarations are uniform priors that
    # frmtmb does not carry. Three means on (0, 300) and one standard
    # deviation on (0, 100) put a flat -3 log(300) - log(100) into every
    # value of log_prob, and nothing into its gradient.
    const = -3 * log(300) - log(100))
})
