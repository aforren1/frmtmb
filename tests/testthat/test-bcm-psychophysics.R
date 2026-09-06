# Lee and Wagenmakers, Bayesian Cognitive Modeling, chapter 12:
# psychophysical functions.
#
# Stan programs adapted from the BSD-3 ports at
# github.com/stan-dev/example-models, directory
# Bayesian_Cognitive_Modeling/CaseStudies/PsychophysicalFunctions.
# Sampling statements become `target += ..._lpdf(...)`; the vague
# normal(0, 1 / sqrt(.001)) regularizers on the group means are dropped
# (see the note at the top of test-bcm-data-analysis.R); and the
# padded-matrix data layout of the originals, with its -99 sentinel for
# a subject who saw fewer stimulus levels, becomes one long vector,
# because that is the shape a formula reads.

bcm_psy1_code <- function() {
  paste(
    "data {",
    "  int<lower=1> N; int<lower=1> S;",
    "  array[N] int<lower=0> r; array[N] int<lower=1> n;",
    "  vector[N] xc; array[N] int<lower=1> subj;",
    "}",
    "parameters {",
    "  real mua; real mub;",
    "  real<lower=0> sigmaa; real<lower=0> sigmab;",
    "  vector[S] alpha; vector[S] beta;",
    "}",
    "model {",
    "  vector[N] theta;",
    "  target += normal_lpdf(alpha | mua, sigmaa);",
    "  target += normal_lpdf(beta | mub, sigmab);",
    "  for (i in 1 : N)",
    "    theta[i] = inv_logit(alpha[subj[i]] + beta[subj[i]] * xc[i]);",
    "  target += binomial_lpmf(r | n, theta);",
    "}",
    sep = "\n")
}

# The contaminant model with the per-cell contaminant rate integrated
# out. See inst/bcm/binomial-extras.R for the integral and for why it
# has to be done rather than sampled: left free, the contaminant fits
# every cell exactly and maximum likelihood puts all the mass on it.
#
# Two adaptations beyond that. The group distribution on the
# contamination rate sits on its LOGIT rather than its probit, for the
# numerical reason inst/bcm/binomial-extras.R gives; and `sigmap`, which
# the original bounds above at 3, is unbounded, because that bound is a
# prior and the fit does not approach it.
bcm_psy2_code <- function() {
  paste(
    "data {",
    "  int<lower=1> N; int<lower=1> S;",
    "  array[N] int<lower=0> r; array[N] int<lower=1> n;",
    "  vector[N] xc; array[N] int<lower=1> subj;",
    "}",
    "parameters {",
    "  real mua; real mub; real mup;",
    "  real<lower=0> sigmaa; real<lower=0> sigmab; real<lower=0> sigmap;",
    "  vector[S] alpha; vector[S] beta; vector[S] logitphi;",
    "}",
    "model {",
    "  target += normal_lpdf(alpha | mua, sigmaa);",
    "  target += normal_lpdf(beta | mub, sigmab);",
    "  target += normal_lpdf(logitphi | mup, sigmap);",
    "  for (i in 1 : N) {",
    "    real theta = inv_logit(alpha[subj[i]] + beta[subj[i]] * xc[i]);",
    "    real phi = inv_logit(logitphi[subj[i]]);",
    "    target += log_mix(phi, -log(n[i] + 1),",
    "                      binomial_lpmf(r[i] | n[i], theta));",
    "  }",
    "}",
    sep = "\n")
}

bcm_psy_stan_data <- function(d) {
  list(N = nrow(d), S = nlevels(d$subj),
       r = as.integer(d$r), n = as.integer(d$n),
       xc = as.numeric(d$xc), subj = as.integer(d$subj))
}

# ---------------------------------------------------------------------
# PsychophysicalFunction1: a logistic psychometric function per subject,
# with the intercept and the slope each drawn from a group distribution.
# ---------------------------------------------------------------------

test_that("PsychophysicalFunction1 is a binomial logit GLMM", {
  d <- bcm_psy_data()
  fit <- frm(r | trials(n) ~ xc + (xc || subj),
             family = binomial(), data = d)
  b <- fixef(fit)$mu
  # longer intervals are called long more often, so the slope is positive
  expect_gt(unname(b["xc"]), 0)
  expect_equal(nrow(ranef(fit)[[1L]]), 8L)
  # the point of estimation at all: the just-noticeable difference is a
  # function of the slope, and every subject has a finite one
  jnd <- qlogis(0.84) / (unname(b["xc"]) + frm_u_term(fit, "0 + xc | subj"))
  expect_true(all(is.finite(jnd)))
  expect_true(all(jnd > 0))
})

test_that("PsychophysicalFunction1 matches its Stan program", {
  skip_unless_stan_identity()
  d <- bcm_psy_data()
  fit <- frm(r | trials(n) ~ xc + (xc || subj),
             family = binomial(), data = d)
  stan_lp_check(
    bcm_psy1_code(),
    data = bcm_psy_stan_data(d),
    fit = fit,
    pars = function(f) {
      b <- fixef(f)$mu
      list(mua = unname(b["(Intercept)"]), mub = unname(b["xc"]),
           sigmaa = frm_sd_term(f, "1 | subj"),
           sigmab = frm_sd_term(f, "0 + xc | subj"),
           alpha = unname(b["(Intercept)"]) + frm_u_term(f, "1 | subj"),
           beta = unname(b["xc"]) + frm_u_term(f, "0 + xc | subj"))
    },
    inner = c("alpha", "beta"),
    const = 0)
})

# ---------------------------------------------------------------------
# PsychophysicalFunction2: the same function with a contaminant process.
# ---------------------------------------------------------------------

bcm_psy2_formula <- function() {
  bf(r | trials(n) ~ xc + (xc || subj), phi ~ 1 + (1 | subj))
}

test_that("PsychophysicalFunction2 adds a contaminant component", {
  skip_unless_bcm("binomial-extras.R")
  d <- bcm_psy_data()
  fit <- frm(bcm_psy2_formula(), family = bcm_contaminant(), data = d)
  ph <- plogis(unname(fixef(fit)$phi["(Intercept)"]))
  # a contamination rate, not a second psychometric function: it stays
  # well away from one, which is what the closed-form integral in
  # inst/bcm/binomial-extras.R buys
  expect_gt(ph, 0)
  expect_lt(ph, 0.5)
  # allowing contaminants makes the psychometric function steeper,
  # which is the book's figure 12.7: the just-noticeable difference
  # shrinks
  plain <- frm(r | trials(n) ~ xc + (xc || subj),
               family = binomial(), data = d)
  expect_gt(unname(fixef(fit)$mu["xc"]), unname(fixef(plain)$mu["xc"]))
})

test_that("PsychophysicalFunction2 matches its Stan program", {
  skip_unless_stan_identity()
  skip_unless_bcm("binomial-extras.R")
  d <- bcm_psy_data()
  fit <- frm(bcm_psy2_formula(), family = bcm_contaminant(), data = d)
  stan_lp_check(
    bcm_psy2_code(),
    data = bcm_psy_stan_data(d),
    fit = fit,
    pars = function(f) {
      b <- fixef(f)$mu
      bp <- fixef(f)$phi
      # three independent blocks on the same grouping factor, two from
      # the mu formula and one from the phi formula, so every one is
      # addressed by the coefficient it carries
      list(mua = unname(b["(Intercept)"]), mub = unname(b["xc"]),
           mup = unname(bp["(Intercept)"]),
           sigmaa = frm_sd_term(f, "1 | subj"),
           sigmab = frm_sd_term(f, "0 + xc | subj"),
           sigmap = frm_sd_term(f, "phi: 1 | subj"),
           alpha = unname(b["(Intercept)"]) + frm_u_term(f, "1 | subj"),
           beta = unname(b["xc"]) + frm_u_term(f, "0 + xc | subj"),
           logitphi = unname(bp["(Intercept)"]) +
             frm_u_term(f, "phi: 1 | subj"))
    },
    inner = c("alpha", "beta", "logitphi"),
    const = 0)
})
