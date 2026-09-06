# Lee and Wagenmakers, Bayesian Cognitive Modeling, chapter 5: data
# analysis (correlation, agreement, change detection, capture-recapture,
# censored counts).
#
# Stan programs adapted from the BSD-3 ports at
# github.com/stan-dev/example-models, directory
# Bayesian_Cognitive_Modeling/ParameterEstimation/DataAnalysis.
#
# Two adaptations recur and are stated once here:
#
#   1. sampling statements become `target += ..._lpdf(...)`, so Stan
#      drops no normalizing constant (see helper-stan.R);
#   2. the vague regularizers the originals put on their location and
#      precision parameters, normal(0, 1 / sqrt(.001)) and
#      gamma(.001, .001), are dropped. They exist so that a sampler has
#      a proper target; maximum likelihood does not need them, and
#      dropping them on both sides is what makes the comparison a
#      LIKELIHOOD identity with a zero constant. Where a prior is kept
#      it is carried through set_prior() on the frmtmb side too; see
#      test-bcm-gaussian.R for that route.

# bcm_corr1_code() is in helper-bcm.R: ESP's Ability and
# OptionalStopping are the same program on different columns.

# Correlation_2: the true values are latent and the observations carry a
# known measurement standard deviation.
bcm_corr2_code <- function() {
  paste(
    "data {",
    "  int<lower=0> n;",
    "  array[n] vector[2] x;",
    "  vector[2] sigmaerror;",
    "}",
    "parameters {",
    "  vector[2] mu;",
    "  vector<lower=0>[2] sigma;",
    "  real<lower=-1, upper=1> r;",
    "  array[n] vector[2] y;",
    "}",
    "model {",
    "  matrix[2, 2] T;",
    "  T[1, 1] = square(sigma[1]);",
    "  T[1, 2] = r * sigma[1] * sigma[2];",
    "  T[2, 1] = T[1, 2];",
    "  T[2, 2] = square(sigma[2]);",
    "  target += multi_normal_lpdf(y | mu, T);",
    "  for (i in 1:n)",
    "    target += normal_lpdf(x[i] | y[i], sigmaerror);",
    "}",
    sep = "\n")
}

bcm_corr_matrix <- function() {
  matrix(c(0.8, 102, 1.0, 98, 0.5, 100, 0.9, 105, 0.7, 103, 0.4, 110,
           1.2, 99, 1.4, 87, 0.6, 113, 1.1, 89, 1.3, 93),
         nrow = 11, ncol = 2, byrow = TRUE)
}

bcm_corr_data <- function(sigmaerror = NULL) {
  x <- bcm_corr_matrix()
  d <- data.frame(x1 = x[, 1], x2 = x[, 2],
                  id = factor(seq_len(nrow(x))))
  if (!is.null(sigmaerror)) {
    d$s1 <- sigmaerror[1]
    d$s2 <- sigmaerror[2]
  }
  d
}

# ---------------------------------------------------------------------
# Correlation_1: the Pearson correlation as a bivariate normal.
# ---------------------------------------------------------------------

test_that("Correlation_1 is a bivariate normal with rescor", {
  d <- bcm_corr_data()
  fit <- frm(mvbf(bf(x1 ~ 1), bf(x2 ~ 1)) + set_rescor(TRUE), data = d)
  # the maximum likelihood correlation IS the sample correlation, which
  # is the frequentist point estimate the book plots as a dashed line
  expect_equal(rescor_matrix(fit)[1, 2],
               stats::cor(d$x1, d$x2), tolerance = 1e-5)
  expect_equal(unname(fixef(fit)$x1_mu), mean(d$x1), tolerance = 1e-6)
  expect_equal(unname(fixef(fit)$x2_mu), mean(d$x2), tolerance = 1e-6)
})

test_that("Correlation_1 matches its Stan program", {
  skip_unless_stan_identity()
  d <- bcm_corr_data()
  fit <- frm(mvbf(bf(x1 ~ 1), bf(x2 ~ 1)) + set_rescor(TRUE), data = d)
  stan_lp_check(
    bcm_corr1_code(),
    data = list(n = nrow(d), x = bcm_corr_matrix()),
    fit = fit,
    pars = function(f) {
      list(mu = c(unname(fixef(f)$x1_mu), unname(fixef(f)$x2_mu)),
           sigma = c(unname(exp(fixef(f)$x1_sigma)),
                     unname(exp(fixef(f)$x2_sigma))),
           r = rescor_matrix(f)[1, 2])
    },
    const = 0)
})

# ---------------------------------------------------------------------
# Correlation_2: the same correlation with known measurement error.
#
# The book's y[i] is the latent true pair and x[i] is what was measured.
# frmtmb spells the latent pair as a CORRELATED RANDOM INTERCEPT shared
# by the two responses, (1 | p | id), and the known measurement standard
# deviation as se(sd, sigma = FALSE), which replaces the residual scale
# instead of adding to it. The correlation the book reports is then the
# correlation of the random-effect block.
#
# `mi(sd)` is the other spelling of a measurement model in this grammar,
# and it does not reach this model: mi() refuses rescor = TRUE, and with
# a correlated random intercept instead there would be two variances per
# response, the latent one and the response's own sigma, which one
# observation per subject cannot separate. That is a property of the
# design, not a gap in the grammar.
# ---------------------------------------------------------------------

bcm_corr2_formula <- function() {
  mvbf(bf(x1 | se(s1, sigma = FALSE) ~ 1 + (1 | p | id)),
       bf(x2 | se(s2, sigma = FALSE) ~ 1 + (1 | p | id)))
}

test_that("Correlation_2 attenuates less than Correlation_1", {
  se <- c(0.03, 1)
  d <- bcm_corr_data(se)
  fit <- frm(bcm_corr2_formula(), data = d)
  r <- stats::cov2cor(unname(VarCorr(fit)[[1L]]))[1, 2]
  raw <- stats::cor(d$x1, d$x2)
  # removing measurement noise from the observed variances can only make
  # the estimated latent correlation more extreme
  expect_lt(r, raw)
  expect_gt(r, -1)
})

test_that("Correlation_2 matches its Stan program", {
  skip_unless_stan_identity()
  se <- c(0.03, 1)
  d <- bcm_corr_data(se)
  fit <- frm(bcm_corr2_formula(), data = d)
  stan_lp_check(
    bcm_corr2_code(),
    data = list(n = nrow(d), x = bcm_corr_matrix(), sigmaerror = se),
    fit = fit,
    pars = function(f) {
      vc <- unname(VarCorr(f)[[1L]])
      mu <- c(unname(fixef(f)$x1_mu), unname(fixef(f)$x2_mu))
      u <- frm_u(f)
      list(mu = mu, sigma = sqrt(diag(vc)),
           r = stats::cov2cor(vc)[1, 2],
           # the Stan program's latent pair is the mean plus this
           # subject's deviation
           y = cbind(mu[1] + u[, 1], mu[2] + u[, 2]))
    },
    inner = "y",
    const = 0)
})

# ---------------------------------------------------------------------
# Kappa: chance-corrected agreement between two rating methods.
#
# Four counts from three underlying rates, which is a small processing
# tree and is the family in inst/bcm/process-trees.R. The map from rates
# to cells is invertible here, so the fit is saturated and reproduces
# the four proportions exactly; what the chapter is about is the derived
# quantity, and bcm_kappa_summary() computes it from the estimates
# rather than from the counts.
# ---------------------------------------------------------------------

bcm_kappa_data <- function() {
  data.frame(y1 = 14L, y2 = 4L, y3 = 5L, y4 = 210L)
}

bcm_kappa_code <- function() {
  paste(
    "data {",
    "  array[4] int<lower=0> y;",
    "}",
    "parameters {",
    "  real<lower=0, upper=1> alpha;",
    "  real<lower=0, upper=1> beta;",
    "  real<lower=0, upper=1> gamma;",
    "}",
    "model {",
    "  vector[4] pi;",
    "  pi[1] = alpha * beta;",
    "  pi[2] = alpha * (1 - beta);",
    "  pi[3] = (1 - alpha) * (1 - gamma);",
    "  pi[4] = (1 - alpha) * gamma;",
    "  target += beta_lpdf(alpha | 1, 1);",
    "  target += beta_lpdf(beta | 1, 1);",
    "  target += beta_lpdf(gamma | 1, 1);",
    "  target += multinomial_lpmf(y | pi);",
    "}",
    sep = "\n")
}

test_that("Kappa reproduces Cohen's point estimate", {
  skip_unless_bcm("process-trees.R")
  d <- bcm_kappa_data()
  fit <- frm(cbind(y1, y2, y3, y4) ~ 1, family = bcm_kappa(), data = d)
  s <- bcm_kappa_summary(fit)
  y <- as.numeric(d[1, ])
  n <- sum(y)
  p0 <- (y[1] + y[4]) / n
  pe <- ((y[1] + y[2]) * (y[1] + y[3]) +
           (y[2] + y[4]) * (y[3] + y[4])) / n^2
  cohen <- (p0 - pe) / (1 - pe)
  # a saturated tree reproduces the sample proportions, so the maximum
  # likelihood kappa IS Cohen's kappa
  expect_equal(s$kappa[1], cohen, tolerance = 1e-5)
  expect_gt(s$kappa[1], 0.7)
})

test_that("Kappa matches its Stan program", {
  skip_unless_stan_identity()
  skip_unless_bcm("process-trees.R")
  d <- bcm_kappa_data()
  fit <- frm(cbind(y1, y2, y3, y4) ~ 1, family = bcm_kappa(), data = d)
  stan_lp_check(
    bcm_kappa_code(),
    data = list(y = as.integer(as.numeric(d[1, ]))),
    fit = fit,
    pars = function(f) {
      s <- bcm_kappa_summary(f)
      list(alpha = s$alpha[1], beta = s$beta[1], gamma = s$gamma[1])
    },
    const = 0)
})

# ---------------------------------------------------------------------
# ChangeDetection: one sequence, two means, and a changepoint nobody
# saw.
#
# The book puts a discrete uniform on the changepoint, and the
# likelihood is the sum over the n - 1 places it could sit. The family
# is in inst/bcm/marginal.R, which also records why the reference port's
# CONTINUOUS tau is not the same model: its log density is piecewise
# constant in tau, so its gradient there is zero almost everywhere.
#
# The case study's own thousand-point series is not embedded; a series
# of the same shape is generated instead, and the test asks the model to
# find the change that was put there.
# ---------------------------------------------------------------------

bcm_changepoint_data <- function(n = 120L, tau = 70L, mu = c(1, 3),
                                 sigma = 1, seed = 11L) {
  set.seed(seed)
  data.frame(t = seq_len(n),
             c = stats::rnorm(n, ifelse(seq_len(n) <= tau, mu[1], mu[2]),
                              sigma))
}

bcm_changepoint_code <- function() {
  paste(
    "data {",
    "  int<lower=3> n; vector[n] c;",
    "}",
    "parameters {",
    "  real mu1; real mu2; real<lower=0> sigma;",
    "}",
    "model {",
    "  vector[n - 1] lp;",
    "  vector[n] ll1;",
    "  vector[n] ll2;",
    "  for (i in 1 : n) {",
    "    ll1[i] = normal_lpdf(c[i] | mu1, sigma);",
    "    ll2[i] = normal_lpdf(c[i] | mu2, sigma);",
    "  }",
    "  for (tau in 1 : (n - 1))",
    "    lp[tau] = -log(n - 1) + sum(ll1[1 : tau])",
    "              + sum(ll2[(tau + 1) : n]);",
    "  target += log_sum_exp(lp);",
    "}",
    sep = "\n")
}

test_that("ChangeDetection finds the change it was given", {
  skip_unless_bcm("marginal.R")
  d <- bcm_changepoint_data()
  fit <- frm(c ~ 1, family = bcm_changepoint(t), data = d)
  dp <- eval_dpars(fit)[["c"]]
  mu <- c(as.numeric(dp[["mu1"]])[1L], as.numeric(dp[["mu2"]])[1L])
  expect_equal(mu, c(1, 3), tolerance = 0.3)
  pr <- latent_probs(fit)
  expect_equal(sum(pr), 1, tolerance = 1e-10)
  # the posterior over the changepoint concentrates where it was put
  expect_lt(abs(as.integer(colnames(pr)[which.max(pr)]) - 70L), 5L)
})

test_that("ChangeDetection matches its Stan program", {
  skip_unless_stan_identity()
  skip_unless_bcm("marginal.R")
  d <- bcm_changepoint_data()
  fit <- frm(c ~ 1, family = bcm_changepoint(t), data = d)
  stan_lp_check(
    bcm_changepoint_code(),
    data = list(n = nrow(d), c = d$c),
    fit = fit,
    pars = function(f) {
      dp <- eval_dpars(f)[["c"]]
      list(mu1 = as.numeric(dp[["mu1"]])[1L],
           mu2 = as.numeric(dp[["mu2"]])[1L],
           sigma = as.numeric(dp[["sigma"]])[1L])
    },
    const = 0)
})

# ---------------------------------------------------------------------
# Planes: REFUSED as a fit, by name.
#
# The model has no free parameter. Ten planes were marked, a later
# sample of five contained four marked ones, and the fleet size has a
# discrete uniform prior; the whole content of the chapter is the
# posterior over that size. Its own Stan port says the same thing in its
# own way, by running with algorithm = "Fixed_param".
#
# frm() estimates parameters, and a model with none is not a fit, so
# there is no formula to write. What CAN be checked is the arithmetic,
# and it is checked against the reference program directly: Stan
# evaluates the same log_sum_exp with an empty parameter vector.
# ---------------------------------------------------------------------

bcm_planes_code <- function() {
  paste(
    "data {",
    "  int<lower=1> x; int<lower=1> n; int<lower=0> k; int<lower=1> tmax;",
    "}",
    "transformed data {",
    "  int tmin = x + n - k;",
    "}",
    "parameters {",
    "}",
    "model {",
    "  vector[tmax - tmin + 1] lp;",
    "  for (i in 1 : (tmax - tmin + 1))",
    "    lp[i] = -log(tmax)",
    "            + hypergeometric_lpmf(k | n, x, tmin + i - 1 - x);",
    "  target += log_sum_exp(lp);",
    "}",
    sep = "\n")
}

test_that("Planes has no parameter to estimate", {
  skip_unless_bcm("marginal.R")
  p <- bcm_planes_posterior(x = 10, n = 5, k = 4, tmax = 50)
  expect_equal(sum(p$prob), 1, tolerance = 1e-12)
  # the smallest fleet consistent with the data is x + n - k
  expect_equal(min(p$t), 11)
  expect_equal(max(p$t), 50)
  # the mode is the Lincoln-Petersen estimate, floor(x n / k): the fleet
  # size that makes four recaptures out of five most likely
  expect_equal(p$t[which.max(p$prob)], floor(10 * 5 / 4))
})

test_that("Planes matches its Stan program", {
  skip_unless_stan_identity()
  skip_unless_bcm("marginal.R")
  p <- bcm_planes_posterior(x = 10, n = 5, k = 4, tmax = 50)
  stan_lp_check(
    bcm_planes_code(),
    data = list(x = 10L, n = 5L, k = 4L, tmax = 50L),
    ours = p$logmarg,
    pars = list(),
    const = 0,
    # there is no parameter, so there is no gradient to make a claim
    # about; the value IS the whole model
    grad = "none")
})

# ---------------------------------------------------------------------
# ChaSaSoon: one observed score of 30 out of 50, and 949 earlier
# attempts known only to have scored between 15 and 25.
#
# A BAND of counts rather than a marginalization. cens("interval") is
# the grammar's word for it and does not reach here: cens() refuses
# every family whose type is "discrete", CDF or not. So the band is the
# family's own business, and inst/bcm/binomial-extras.R spells it as
# "the response is the smallest count consistent with the row and
# vint(hi) the largest", which is the exact density when the two agree
# and the band probability otherwise. weights() carries the repeats.
# ---------------------------------------------------------------------

bcm_chasasoon_data <- function() {
  data.frame(z = c(30L, 15L), hi = c(30L, 25L),
             w = c(1, 949), n = 50L)
}

bcm_chasasoon_code <- function() {
  paste(
    "data {",
    "  int<lower=0> nfails; int<lower=0> n; int<lower=0> z;",
    "}",
    "parameters {",
    "  real<lower=0.25, upper=1> theta;",
    "}",
    "model {",
    "  target += binomial_lpmf(z | n, theta);",
    "  target += nfails * log(binomial_cdf(25 | n, theta)",
    "                         - binomial_cdf(14 | n, theta));",
    "}",
    sep = "\n")
}

bcm_chasasoon_formula <- function() {
  bf(z | trials(n) + vint(hi) + weights(w) ~ 1)
}

test_that("ChaSaSoon is a band of counts with a repeat count", {
  skip_unless_bcm("binomial-extras.R")
  d <- bcm_chasasoon_data()
  fit <- frm(bcm_chasasoon_formula(), family = bcm_binomial_band(),
             data = d)
  theta <- plogis(unname(fixef(fit)$mu))
  # the book's answer is about 0.34: well below the 30/50 of the one
  # success, and above the 0.25 a random guesser would score
  expect_gt(theta, 0.25)
  expect_lt(theta, 0.45)
})

test_that("a band of one count is the ordinary binomial density", {
  skip_unless_bcm("binomial-extras.R")
  d <- data.frame(z = c(5L, 7L), hi = c(5L, 7L), w = 1, n = 10L)
  band <- frm(bcm_chasasoon_formula(), family = bcm_binomial_band(),
              data = d)
  plain <- frm(z | trials(n) ~ 1, family = binomial(), data = d)
  expect_equal(as.numeric(logLik(band)), as.numeric(logLik(plain)),
               tolerance = 1e-8)
})

test_that("cens() refuses a discrete family whatever CDF it supplies", {
  skip_unless_bcm("binomial-extras.R")
  d <- bcm_chasasoon_data()
  d$cc <- c("none", "interval")
  d$ub <- c(NA_integer_, 25L)
  # the core binomial is refused for having no CDF ...
  expect_error(
    frm(z | trials(n) + cens(cc, ub) + weights(w) ~ 1,
        family = binomial(), data = d),
    "family with a CDF")
  # ... and a family that supplies one is refused for being discrete,
  # which is why the band above is written into the family instead
  expect_error(
    frm(z | trials(n) + vint(hi) + cens(cc, ub) + weights(w) ~ 1,
        family = bcm_binomial_band(), data = d),
    "discrete families")
})

test_that("ChaSaSoon matches its Stan program", {
  skip_unless_stan_identity()
  skip_unless_bcm("binomial-extras.R")
  d <- bcm_chasasoon_data()
  fit <- frm(bcm_chasasoon_formula(), family = bcm_binomial_band(),
             data = d)
  stan_lp_check(
    bcm_chasasoon_code(),
    data = list(nfails = 949L, n = 50L, z = 30L),
    fit = fit,
    pars = function(f) list(theta = plogis(unname(fixef(f)$mu))),
    # the original's <lower=.25, upper=1> is a uniform prior Stan does
    # not add to the target, and the estimate is inside it
    const = 0)
})
