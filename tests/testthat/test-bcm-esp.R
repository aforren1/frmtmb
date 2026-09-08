# Lee and Wagenmakers, Bayesian Cognitive Modeling, chapter 16:
# extrasensory perception.
#
# Stan programs adapted from the BSD-3 ports at
# github.com/stan-dev/example-models, directory
# Bayesian_Cognitive_Modeling/CaseStudies/ESP.
#
# Ability and OptionalStopping are Correlation_1's program applied to a
# different pair of columns, so their frmtmb spelling is the one
# test-bcm-data-analysis.R validates against that program: a bivariate
# normal with set_rescor(TRUE). Only OptionalStopping's data is embedded
# here; Ability's is a hundred pairs of proportions that would add
# nothing the correlation model has not already been checked on.
#
# Extraversion is the chapter's own model and is here in full.
#
# What the chapter is FOR is a Bayes factor at r = 0 or delta = 0, and
# that is not here. See test-bcm-model-selection.R for why, and
# dev/bcm-findings.md for the handover to frmtmb.sample.

# ---------------------------------------------------------------------
# OptionalStopping: the correlation between an experiment's sample size
# and its effect size, across nine Bem experiments.
# ---------------------------------------------------------------------

bcm_optstop_data <- function() {
  data.frame(N = c(100, 150, 97, 99, 100, 150, 200, 100, 50),
             E = c(0.25, 0.20, 0.25, 0.20, 0.22, 0.15, 0.09, 0.19, 0.42))
}

test_that("OptionalStopping is a bivariate normal on nine experiments", {
  d <- bcm_optstop_data()
  fit <- frm(mvbf(bf(N ~ 1), bf(E ~ 1)) + set_rescor(TRUE), data = d)
  r <- rescor_matrix(fit)[1, 2]
  expect_equal(r, stats::cor(d$N, d$E), tolerance = 1e-5)
  # the chapter's point: a strong NEGATIVE correlation, which is what
  # optional stopping would produce
  expect_lt(r, -0.7)
})

test_that("OptionalStopping matches Correlation_1's Stan program", {
  skip_unless_stan_identity()
  d <- bcm_optstop_data()
  fit <- frm(mvbf(bf(N ~ 1), bf(E ~ 1)) + set_rescor(TRUE), data = d)
  stan_lp_check(
    bcm_corr1_code(),
    data = list(n = nrow(d), x = cbind(d$N, d$E)),
    fit = fit,
    pars = function(f) {
      list(mu = c(unname(fixef(f)$N_mu), unname(fixef(f)$E_mu)),
           sigma = c(unname(exp(fixef(f)$N_sigma)),
                     unname(exp(fixef(f)$E_sigma))),
           r = rescor_matrix(f)[1, 2])
    },
    const = 0)
})

# ---------------------------------------------------------------------
# Extraversion: a latent pair per person, one half driving a binomial
# rate and the other an extraversion score with known measurement error.
#
# The chapter writes the second response as normal around
# 100 * Phi(theta_2). Dividing it by 100 makes the mean Phi(theta_2)
# exactly, so BOTH responses are Phi of a linear predictor and the
# latent pair is an ordinary correlated random intercept shared by two
# responses, (1 | p | id). The two families are the probit binomial and
# the probit gaussian of inst/bcm/binomial-extras.R, and the known
# measurement standard deviation rides on se(sd), scaled by the same
# 100.
#
# The known measurement standard deviation rides on se(sd), which is the
# term that means it: a family is given the term by declaring that it
# reads it. It rode on vreal(sd) until the core stopped gating se() on
# the family NAME. See inst/bcm/binomial-extras.R and
# dev/custom-findings.md.
# ---------------------------------------------------------------------

bcm_extraversion_data <- function() {
  k <- c(36, 32, 36, 36, 28, 40, 40, 24, 36, 36, 28, 40, 28, 36, 20, 24,
         24, 16, 20, 32, 40, 32, 36, 24, 28, 44, 40, 36, 40, 32, 32, 40,
         28, 20, 24, 32, 24, 24, 20, 28, 24, 28, 28, 32, 20, 44, 16, 36,
         32, 28, 24, 32, 40, 28, 32, 32, 28, 24, 28, 40, 28, 20, 20, 20,
         24, 24, 36, 28, 20, 20, 40, 32, 20, 36, 28, 28, 24, 20, 28, 32,
         48, 24, 32, 32, 40, 40, 40, 36, 36, 32, 20, 28, 40, 32, 20, 20,
         16, 16, 28, 40)
  x <- c(50, 80, 79, 56, 50, 80, 53, 84, 74, 67, 50, 45, 62, 65, 71, 71,
         68, 63, 67, 58, 72, 73, 63, 54, 63, 70, 81, 71, 66, 74, 70, 84,
         66, 73, 78, 64, 54, 74, 62, 71, 70, 79, 66, 64, 62, 63, 60, 56,
         72, 72, 79, 67, 46, 67, 77, 55, 63, 44, 84, 65, 41, 62, 64, 51,
         46, 53, 26, 67, 73, 39, 62, 59, 75, 65, 60, 69, 63, 69, 55, 63,
         86, 70, 67, 54, 80, 71, 71, 55, 57, 41, 56, 78, 58, 76, 54, 50,
         61, 60, 32, 67)
  data.frame(k = k, nt = 60L, xs = x / 100, sx = 3 / 100,
             id = factor(seq_along(k)))
}

bcm_extraversion_formula <- function() {
  mvbf(bf(k | trials(nt) ~ 1 + (1 | p | id)) + binomial(link = "probit"),
       bf(xs | se(sx) ~ 1 + (1 | p | id)) + gaussian(link = "probit"))
}

bcm_extraversion_code <- function() {
  paste(
    "data {",
    "  int<lower=1> S; int<lower=1> nt;",
    "  array[S] int<lower=0> k; vector[S] xs; real<lower=0> sx;",
    "}",
    "parameters {",
    "  vector[2] mu;",
    "  vector<lower=0>[2] sigma;",
    "  real<lower=-1, upper=1> r;",
    "  array[S] vector[2] theta;",
    "}",
    "model {",
    "  matrix[2, 2] T;",
    "  T[1, 1] = square(sigma[1]);",
    "  T[1, 2] = r * sigma[1] * sigma[2];",
    "  T[2, 1] = T[1, 2];",
    "  T[2, 2] = square(sigma[2]);",
    "  target += multi_normal_lpdf(theta | mu, T);",
    "  for (i in 1 : S) {",
    "    target += binomial_lpmf(k[i] | nt, Phi(theta[i, 1]));",
    "    target += normal_lpdf(xs[i] | Phi(theta[i, 2]), sx);",
    "  }",
    "}",
    sep = "\n")
}

test_that("Extraversion is a correlated latent pair across responses", {
  skip_unless_bcm("binomial-extras.R")
  d <- bcm_extraversion_data()
  fit <- frm(bcm_extraversion_formula(), data = d)
  vc <- unname(VarCorr(fit)[[1L]])
  r <- stats::cov2cor(vc)[1, 2]
  expect_equal(dim(vc), c(2L, 2L))
  # the chapter's conclusion is that psychic performance and
  # extraversion are barely related
  expect_lt(abs(r), 0.5)
  # both latent scales are estimated and finite
  expect_true(all(sqrt(diag(vc)) > 0))
})

test_that("Extraversion matches its Stan program", {
  skip_unless_stan_identity()
  skip_unless_bcm("binomial-extras.R")
  d <- bcm_extraversion_data()
  fit <- frm(bcm_extraversion_formula(), data = d)
  stan_lp_check(
    bcm_extraversion_code(),
    data = list(S = nrow(d), nt = 60L, k = as.integer(d$k),
                xs = d$xs, sx = 0.03),
    fit = fit,
    pars = function(f) {
      vc <- unname(VarCorr(f)[[1L]])
      mu <- c(unname(fixef(f)$k_mu), unname(fixef(f)$xs_mu))
      u <- frm_u(f)
      list(mu = mu, sigma = sqrt(diag(vc)),
           r = stats::cov2cor(vc)[1, 2],
           theta = cbind(mu[1] + u[, 1], mu[2] + u[, 2]))
    },
    inner = "theta",
    const = 0)
})
