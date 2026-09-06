# Lee and Wagenmakers, Bayesian Cognitive Modeling, chapter 3:
# inferring a rate.
#
# The Stan programs below are adapted from the BSD-3 licensed ports at
# github.com/stan-dev/example-models, directory
# Bayesian_Cognitive_Modeling/ParameterEstimation/Binomial. Two
# adaptations are made throughout, and only these two:
#
#   1. every `y ~ dist(...)` becomes `target += dist_lpmf(y | ...)`, so
#      Stan drops no normalizing constant (see helper-stan.R);
#   2. generated quantities are dropped, because log_prob() never
#      reaches them and a predictive is checked from R instead.
#
# The models are the book's published equations.

# Rate_1.stan and Rate_3.stan / Rate_5.stan differ only in how many
# counts share the rate, and the frmtmb spelling differs only in how
# many rows the data frame has. One program covers all three.
bcm_rate_common_code <- function() {
  paste(
    "data {",
    "  int<lower=1> m;",
    "  array[m] int<lower=1> n;",
    "  array[m] int<lower=0> k;",
    "}",
    "parameters {",
    "  real<lower=0, upper=1> theta;",
    "}",
    "model {",
    "  target += beta_lpdf(theta | 1, 1);",
    "  target += binomial_lpmf(k | n, theta);",
    "}",
    sep = "\n")
}

# Rate_2.stan: one rate per group, no pooling.
bcm_rate_two_code <- function() {
  paste(
    "data {",
    "  int<lower=1> n1; int<lower=1> n2;",
    "  int<lower=0> k1; int<lower=0> k2;",
    "}",
    "parameters {",
    "  real<lower=0, upper=1> theta1;",
    "  real<lower=0, upper=1> theta2;",
    "}",
    "model {",
    "  target += beta_lpdf(theta1 | 1, 1);",
    "  target += beta_lpdf(theta2 | 1, 1);",
    "  target += binomial_lpmf(k1 | n1, theta1);",
    "  target += binomial_lpmf(k2 | n2, theta2);",
    "}",
    sep = "\n")
}

# ---------------------------------------------------------------------
# Rate_1: k = 5 successes out of n = 10.
# ---------------------------------------------------------------------

test_that("Rate_1 is an intercept-only binomial", {
  d <- data.frame(k = 5, n = 10)
  fit <- frm(k | trials(n) ~ 1, family = binomial(), data = d)
  # the book's theta, read back off the logit scale
  expect_equal(plogis(unname(fixef(fit)$mu)), 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(logLik(fit)), dbinom(5, 10, 0.5, log = TRUE),
               tolerance = 1e-8)
})

test_that("Rate_1 matches its Stan program", {
  skip_unless_stan_identity()
  d <- data.frame(k = 5, n = 10)
  fit <- frm(k | trials(n) ~ 1, family = binomial(), data = d)
  res <- stan_lp_check(
    bcm_rate_common_code(),
    data = list(m = 1L, n = as.array(10L), k = as.array(5L)),
    fit = fit,
    pars = function(f) list(theta = plogis(unname(fixef(f)$mu))),
    # beta(1, 1) is the uniform density on the unit interval, so its log
    # density is 0 and frmtmb's flat objective is the same function.
    const = 0)
  expect_equal(res$measured_const, 0, tolerance = 1e-6)
})

# ---------------------------------------------------------------------
# Rate_2: two rates, k1 = 5 / n1 = 10 and k2 = 7 / n2 = 10. The book's
# quantity is delta = theta1 - theta2.
# ---------------------------------------------------------------------

test_that("Rate_2 is a binomial with one rate per group", {
  d <- data.frame(k = c(5, 7), n = c(10, 10),
                  g = factor(c("g1", "g2")))
  fit <- frm(k | trials(n) ~ 0 + g, family = binomial(), data = d)
  th <- plogis(unname(fixef(fit)$mu))
  expect_equal(th, c(0.5, 0.7), tolerance = 1e-6)
  expect_equal(th[1] - th[2], -0.2, tolerance = 1e-6)
})

test_that("Rate_2 matches its Stan program", {
  skip_unless_stan_identity()
  d <- data.frame(k = c(5, 7), n = c(10, 10),
                  g = factor(c("g1", "g2")))
  fit <- frm(k | trials(n) ~ 0 + g, family = binomial(), data = d)
  stan_lp_check(
    bcm_rate_two_code(),
    data = list(n1 = 10L, n2 = 10L, k1 = 5L, k2 = 7L),
    fit = fit,
    pars = function(f) {
      th <- plogis(frm_b(f))
      list(theta1 = th[1], theta2 = th[2])
    },
    const = 0)
})

# ---------------------------------------------------------------------
# Rate_3: one rate for both counts.
# ---------------------------------------------------------------------

test_that("Rate_3 pools two counts into one rate", {
  d <- data.frame(k = c(5, 7), n = c(10, 10))
  fit <- frm(k | trials(n) ~ 1, family = binomial(), data = d)
  expect_equal(plogis(unname(fixef(fit)$mu)), 0.6, tolerance = 1e-6)
})

test_that("Rate_3 matches its Stan program", {
  skip_unless_stan_identity()
  d <- data.frame(k = c(5, 7), n = c(10, 10))
  fit <- frm(k | trials(n) ~ 1, family = binomial(), data = d)
  stan_lp_check(
    bcm_rate_common_code(),
    data = list(m = 2L, n = c(10L, 10L), k = c(5L, 7L)),
    fit = fit,
    pars = function(f) list(theta = plogis(unname(fixef(f)$mu))),
    const = 0)
})

# ---------------------------------------------------------------------
# Rate_4: k = 1 out of n = 15, with prior and posterior predictives.
#
# Rate_4.stan declares a second parameter, `thetaprior`, that no data
# reaches: it exists so that the sampler draws from the prior and the
# script can histogram a prior predictive. Maximum likelihood has
# nothing to say about it, so the adapted program drops it and the
# predictive is drawn from R.
# ---------------------------------------------------------------------

test_that("Rate_4 is Rate_1 with a different count", {
  d <- data.frame(k = 1, n = 15)
  fit <- frm(k | trials(n) ~ 1, family = binomial(), data = d)
  expect_equal(plogis(unname(fixef(fit)$mu)), 1 / 15, tolerance = 1e-6)
})

test_that("Rate_4's posterior predictive is a binomial at the estimate", {
  d <- data.frame(k = 1, n = 15)
  fit <- frm(k | trials(n) ~ 1, family = binomial(), data = d)
  set.seed(4)
  draws <- unlist(simulate(fit, nsim = 4000))
  expect_true(all(draws >= 0 & draws <= 15))
  # E[k] = n * theta_hat = 1, which is the count that was observed
  expect_equal(mean(draws), 1, tolerance = 0.15)
})

test_that("Rate_4 matches its Stan program", {
  skip_unless_stan_identity()
  d <- data.frame(k = 1, n = 15)
  fit <- frm(k | trials(n) ~ 1, family = binomial(), data = d)
  stan_lp_check(
    bcm_rate_common_code(),
    data = list(m = 1L, n = as.array(15L), k = as.array(1L)),
    fit = fit,
    pars = function(f) list(theta = plogis(unname(fixef(f)$mu))),
    const = 0)
})

# ---------------------------------------------------------------------
# Rate_5: k1 = 0 of 10 and k2 = 10 of 10 share one rate. The book uses
# it to show a posterior predictive that misses the data badly, which
# maximum likelihood reproduces exactly: theta_hat is 0.5 and neither
# count is anywhere near 5.
# ---------------------------------------------------------------------

test_that("Rate_5 pools two contradictory counts", {
  d <- data.frame(k = c(0, 10), n = c(10, 10))
  fit <- frm(k | trials(n) ~ 1, family = binomial(), data = d)
  expect_equal(plogis(unname(fixef(fit)$mu)), 0.5, tolerance = 1e-6)
  set.seed(5)
  draws <- do.call(rbind, simulate(fit, nsim = 2000))
  # the predictive puts almost no mass on the observed pair
  hit <- mean(draws[, 1] == 0 & draws[, 2] == 10)
  expect_lt(hit, 0.01)
})

test_that("Rate_5 matches its Stan program", {
  skip_unless_stan_identity()
  d <- data.frame(k = c(0, 10), n = c(10, 10))
  fit <- frm(k | trials(n) ~ 1, family = binomial(), data = d)
  stan_lp_check(
    bcm_rate_common_code(),
    data = list(m = 2L, n = c(10L, 10L), k = c(0L, 10L)),
    fit = fit,
    pars = function(f) list(theta = plogis(unname(fixef(f)$mu))),
    const = 0)
})

# ---------------------------------------------------------------------
# Survey: five surveys came back with 16, 18, 22, 25 and 27 responses,
# and neither the number sent nor the return rate is known.
#
# The number sent is a whole number, and a maximum over a whole number
# is not an estimate, so it is SUMMED OUT of the whole response at once.
# That is the structured protocol's loglik slot, and the family is in
# inst/bcm/marginal.R. The Stan program does the same sum, for the
# different reason that Stan cannot sample a discrete parameter.
# ---------------------------------------------------------------------

bcm_survey_data <- function() {
  data.frame(k = c(16L, 18L, 22L, 25L, 27L))
}

bcm_survey_code <- function() {
  paste(
    "data {",
    "  int<lower=0> nmax; int<lower=0> m;",
    "  array[m] int<lower=0> k;",
    "}",
    "parameters {",
    "  real<lower=0, upper=1> theta;",
    "}",
    "model {",
    "  int nmin = max(k);",
    "  vector[nmax - nmin + 1] lp;",
    "  for (i in 1 : (nmax - nmin + 1))",
    "    lp[i] = -log(nmax) + binomial_lpmf(k | nmin + i - 1, theta);",
    "  target += log_sum_exp(lp);",
    "}",
    sep = "\n")
}

test_that("Survey sums out the number of surveys sent", {
  skip_unless_bcm("marginal.R")
  d <- bcm_survey_data()
  fit <- frm(k ~ 1, family = bcm_survey(nmax = 500L), data = d)
  theta <- plogis(unname(fixef(fit)$mu))
  expect_gt(theta, 0)
  expect_lt(theta, 1)
  # the posterior over the number sent is a proper distribution over
  # max(k)..nmax, and its mode sits above the largest observed return
  pr <- latent_probs(fit)
  expect_equal(sum(pr), 1, tolerance = 1e-10)
  expect_equal(ncol(pr), 500L - 27L + 1L)
  nhat <- as.integer(colnames(pr)[which.max(pr)])
  expect_gte(nhat, 27L)
  # the book's own reading: rate and number trade off, so the number is
  # very poorly determined and the distribution is wide
  expect_lt(max(pr), 0.2)
})

test_that("Survey has no row-wise density to offer", {
  skip_unless_bcm("marginal.R")
  d <- bcm_survey_data()
  fit <- frm(k ~ 1, family = bcm_survey(nmax = 500L), data = d)
  expect_error(residuals(fit, type = "osa"), "osa")
})

test_that("Survey matches its Stan program", {
  skip_unless_stan_identity()
  skip_unless_bcm("marginal.R")
  d <- bcm_survey_data()
  fit <- frm(k ~ 1, family = bcm_survey(nmax = 500L), data = d)
  stan_lp_check(
    bcm_survey_code(),
    data = list(nmax = 500L, m = nrow(d), k = as.integer(d$k)),
    fit = fit,
    pars = function(f) list(theta = plogis(unname(fixef(f)$mu))),
    const = 0)
})
