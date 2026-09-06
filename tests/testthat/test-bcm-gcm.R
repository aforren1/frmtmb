# Lee and Wagenmakers, Bayesian Cognitive Modeling, chapter 17: the
# generalized context model.
#
# Stan programs adapted from the BSD-3 ports at
# github.com/stan-dev/example-models, directory
# Bayesian_Cognitive_Modeling/CaseStudies/GCM, on Kruschke's (1993)
# data: eight stimuli in two categories, forty subjects, eight
# presentations of each stimulus each.
#
# The likelihood does not factorize over rows. A stimulus is called
# category A with a probability that normalizes its similarity to every
# stimulus in the set, so one row's contribution reads the whole set,
# and that is the structured protocol's loglik slot. The family is in
# inst/bcm/similarity.R.
#
# One adaptation: the original bounds the generalization gradient c on
# (0, 5). That bound is a prior rather than part of the likelihood, so
# the family puts c on a log link and lets the fit go where it goes;
# the test asserts that it stays inside the interval anyway.

bcm_gcm_d1 <- function() {
  matrix(c(
    0, 0, 1.00005, 1.00005, 1.95205, 1.95205, 3.1131, 3.1131,
    0, 0, 1.00005, 1.00005, 1.95205, 1.95205, 3.1131, 3.1131,
    1.00005, 1.00005, 0, 0, 0.952, 0.952, 2.11305, 2.11305,
    1.00005, 1.00005, 0, 0, 0.952, 0.952, 2.11305, 2.11305,
    1.95205, 1.95205, 0.952, 0.952, 0, 0, 1.16105, 1.16105,
    1.95205, 1.95205, 0.952, 0.952, 0, 0, 1.16105, 1.16105,
    3.1131, 3.1131, 2.11305, 2.11305, 1.16105, 1.16105, 0, 0,
    3.1131, 3.1131, 2.11305, 2.11305, 1.16105, 1.16105, 0, 0),
    nrow = 8, ncol = 8, byrow = TRUE)
}

bcm_gcm_d2 <- function() {
  matrix(c(
    0, 1.17505, 0.829, 2.23, 0.829, 2.23, 0, 1.17505,
    1.17505, 0, 2.00405, 1.05495, 2.00405, 1.05495, 1.17505, 0,
    0.829, 2.00405, 0, 3.059, 0, 3.059, 0.829, 2.00405,
    2.23, 1.05495, 3.059, 0, 3.059, 0, 2.23, 1.05495,
    0.829, 2.00405, 0, 3.059, 0, 3.059, 0.829, 2.00405,
    2.23, 1.05495, 3.059, 0, 3.059, 0, 2.23, 1.05495,
    0, 1.17505, 0.829, 2.23, 0.829, 2.23, 0, 1.17505,
    1.17505, 0, 2.00405, 1.05495, 2.00405, 1.05495, 1.17505, 0),
    nrow = 8, ncol = 8, byrow = TRUE)
}

bcm_gcm_a <- function() c(1L, 1L, 1L, 2L, 1L, 2L, 2L, 2L)

bcm_gcm_family <- function() {
  bcm_gcm(stimulus = stim, d1 = bcm_gcm_d1(), d2 = bcm_gcm_d2(),
          a = bcm_gcm_a())
}

# One Stan program for both models, with the number of subjects switched
# by data: nsubj = 1 and every row in one subject is GCM_1.
bcm_gcm_code <- function() {
  paste(
    "data {",
    "  int<lower=1> N; int<lower=1> nstim; int<lower=1> nsubj;",
    "  array[N] int<lower=0> y; array[N] int<lower=1> t;",
    "  array[N] int<lower=1> stim; array[N] int<lower=1> subj;",
    "  array[nstim] int<lower=1, upper=2> a;",
    "  matrix[nstim, nstim] d1; matrix[nstim, nstim] d2;",
    "  real<lower=0, upper=1> b;",
    "}",
    "parameters {",
    "  vector<lower=0>[nsubj] c;",
    "  vector<lower=0, upper=1>[nsubj] w;",
    "}",
    "model {",
    "  for (i in 1 : N) {",
    "    real S1 = 0;",
    "    real S2 = 0;",
    "    for (j in 1 : nstim) {",
    "      real s = exp(-c[subj[i]]",
    "                   * (w[subj[i]] * d1[stim[i], j]",
    "                      + (1 - w[subj[i]]) * d2[stim[i], j]));",
    "      if (a[j] == 1) S1 += s; else S2 += s;",
    "    }",
    "    target += binomial_lpmf(y[i] | t[i],",
    "                            b * S1 / (b * S1 + (1 - b) * S2));",
    "  }",
    "}",
    sep = "\n")
}

bcm_gcm_stan_data <- function(d, nsubj, subj) {
  list(N = nrow(d), nstim = 8L, nsubj = nsubj,
       y = as.integer(d$y), t = as.integer(d$t),
       stim = as.integer(d$stim), subj = as.integer(subj),
       a = bcm_gcm_a(), d1 = bcm_gcm_d1(), d2 = bcm_gcm_d2(), b = 0.5)
}

# ---------------------------------------------------------------------
# GCM_1: the forty subjects pooled into one set of counts.
# ---------------------------------------------------------------------

bcm_gcm1_data <- function() {
  data.frame(y = c(245L, 218L, 255L, 126L, 182L, 71L, 102L, 65L),
             t = 320L, stim = 1:8)
}

test_that("GCM_1 estimates a gradient and an attention weight", {
  skip_unless_bcm(c("binomial-extras.R", "similarity.R"))
  d <- bcm_gcm1_data()
  fit <- frm(y | trials(t) ~ 1, family = bcm_gcm_family(), data = d)
  dp <- eval_dpars(fit)[["y"]]
  cc <- as.numeric(dp[["c"]])[1L]
  ww <- as.numeric(dp[["w"]])[1L]
  # the book's headline numbers: a generalization gradient of about 1.6
  # and attention split unevenly between the two dimensions
  expect_gt(cc, 0)
  expect_lt(cc, 5)
  expect_gt(ww, 0)
  expect_lt(ww, 1)
  # the fitted category-A counts track the observed ones
  expect_gt(stats::cor(fitted(fit), d$y), 0.95)
})

test_that("GCM_1 matches its Stan program", {
  skip_unless_stan_identity()
  skip_unless_bcm(c("binomial-extras.R", "similarity.R"))
  d <- bcm_gcm1_data()
  fit <- frm(y | trials(t) ~ 1, family = bcm_gcm_family(), data = d)
  stan_lp_check(
    bcm_gcm_code(),
    data = bcm_gcm_stan_data(d, 1L, rep(1L, nrow(d))),
    fit = fit,
    pars = function(f) {
      dp <- eval_dpars(f)[["y"]]
      list(c = as.array(as.numeric(dp[["c"]])[1L]),
           w = as.array(as.numeric(dp[["w"]])[1L]))
    },
    # the beta(1, 1) on w is the uniform density on the unit interval
    const = 0)
})

# ---------------------------------------------------------------------
# GCM_2: one gradient and one attention weight per subject.
#
# Ten of the forty subjects, which is what exercises the per-subject
# parameters without carrying a 320-cell table into a test file. The
# structure the model has is the same at any number of subjects.
# ---------------------------------------------------------------------

bcm_gcm2_matrix <- function() {
  matrix(c(5, 5, 7, 7, 8, 8, 6, 8, 7, 8,
           6, 6, 7, 7, 6, 7, 5, 8, 8, 8,
           6, 8, 8, 6, 8, 8, 5, 6, 8, 8,
           2, 0, 6, 3, 7, 3, 3, 4, 4, 7,
           4, 3, 1, 5, 3, 4, 3, 4, 7, 4,
           1, 1, 2, 4, 3, 0, 1, 0, 0, 1,
           3, 5, 2, 0, 1, 0, 4, 1, 3, 1,
           1, 2, 3, 1, 0, 0, 1, 0, 0, 0),
         nrow = 8, ncol = 10, byrow = TRUE)
}

bcm_gcm2_data <- function() {
  y <- bcm_gcm2_matrix()
  data.frame(y = as.vector(y), t = 8L,
             stim = rep(seq_len(nrow(y)), times = ncol(y)),
             subj = factor(rep(seq_len(ncol(y)), each = nrow(y))))
}

test_that("GCM_2 gives every subject a gradient and a weight", {
  skip_unless_bcm(c("binomial-extras.R", "similarity.R"))
  d <- bcm_gcm2_data()
  fit <- frm(bf(y | trials(t) ~ 0 + subj, w ~ 0 + subj),
             family = bcm_gcm_family(), data = d)
  expect_length(fixef(fit)[["c"]], 10L)
  expect_length(fixef(fit)[["w"]], 10L)
  # the case study's point is that subjects differ in BOTH, so neither
  # set of estimates collapses to a common value
  expect_gt(stats::sd(exp(frm_b(fit, "c"))), 0)
  expect_gt(stats::sd(plogis(frm_b(fit, "w"))), 0)
})

test_that("GCM_2 matches its Stan program", {
  skip_unless_stan_identity()
  skip_unless_bcm(c("binomial-extras.R", "similarity.R"))
  d <- bcm_gcm2_data()
  fit <- frm(bf(y | trials(t) ~ 0 + subj, w ~ 0 + subj),
             family = bcm_gcm_family(), data = d)
  stan_lp_check(
    bcm_gcm_code(),
    data = bcm_gcm_stan_data(d, 10L, d$subj),
    fit = fit,
    pars = function(f) list(c = exp(frm_b(f, "c")),
                            w = plogis(frm_b(f, "w"))),
    const = 0)
})
