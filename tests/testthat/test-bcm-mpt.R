# Lee and Wagenmakers, Bayesian Cognitive Modeling, chapter 14:
# multinomial processing trees.
#
# Stan programs adapted from the BSD-3 ports at
# github.com/stan-dev/example-models, directory
# Bayesian_Cognitive_Modeling/CaseStudies/MPT, which holds five of them.
#
# FOUR OF THE FIVE ARE ONE MODEL. MPT_2 through MPT_5 differ only in how
# the multivariate normal over the individual differences is written:
# MPT_2 uses parameter expansion and a Wishart prior, MPT_3 an LKJ prior
# on a correlation matrix, MPT_4 its Cholesky factor, MPT_5 a
# non-centered reparameterization of the same. The README that ships
# with them says as much: these are sampler ergonomics, not models.
# Their likelihood is identical, so frmtmb fits them ONCE, as
# `1 + (1 | p | id)` on each of the three probits, and the correlation
# of that block IS the latent-trait correlation the case study reports.
#
# The family is in inst/bcm/process-trees.R.

bcm_mpt_counts <- function() {
  matrix(c(2, 4, 4, 10, 2, 1, 3, 14, 2, 2, 5, 11, 6, 0, 4, 10,
           1, 0, 4, 15, 1, 0, 2, 17, 1, 2, 4, 13, 4, 1, 6, 9,
           5, 1, 4, 10, 1, 0, 9, 10, 5, 0, 3, 12, 0, 1, 6, 13,
           1, 5, 7, 7, 1, 1, 4, 14, 2, 2, 3, 13, 2, 1, 5, 12,
           2, 0, 6, 12, 1, 0, 5, 14, 2, 1, 8, 9, 3, 0, 2, 15,
           1, 2, 3, 14),
         nrow = 21, ncol = 4, byrow = TRUE)
}

bcm_mpt_data <- function() {
  k <- bcm_mpt_counts()
  d <- data.frame(id = factor(seq_len(nrow(k))))
  d$k1 <- k[, 1]
  d$k2 <- k[, 2]
  d$k3 <- k[, 3]
  d$k4 <- k[, 4]
  d
}

# MPT_1: one tree for one aggregated respondent.
bcm_mpt1_code <- function() {
  paste(
    "data {",
    "  array[4] int<lower=0> k;",
    "}",
    "parameters {",
    "  real<lower=0, upper=1> c;",
    "  real<lower=0, upper=1> r;",
    "  real<lower=0, upper=1> u;",
    "}",
    "model {",
    "  vector[4] theta;",
    "  theta[1] = c * r;",
    "  theta[2] = (1 - c) * u ^ 2;",
    "  theta[3] = (1 - c) * 2 * u * (1 - u);",
    "  theta[4] = c * (1 - r) + (1 - c) * (1 - u) ^ 2;",
    "  target += beta_lpdf(c | 1, 1);",
    "  target += beta_lpdf(r | 1, 1);",
    "  target += beta_lpdf(u | 1, 1);",
    "  target += multinomial_lpmf(k | theta);",
    "}",
    sep = "\n")
}

# MPT_2 to MPT_5, once. The covariance of the individual differences
# rides as DATA in its Cholesky form, which is helper-rl.R's arrangement
# and keeps `log_prob` free of a constrained parameter: what the
# gradient check has to say is about the latent traits, and the outer
# gradient of the joint is not zero at a Laplace optimum anyway.
bcm_mpt_hier_code <- function() {
  paste(
    "data {",
    "  int<lower=1> S;",
    "  array[S, 4] int<lower=0> k;",
    "  matrix[3, 3] L;",
    "}",
    "parameters {",
    "  real muchat; real murhat; real muuhat;",
    "  matrix[S, 3] delta;",
    "}",
    "model {",
    "  for (i in 1 : S) {",
    "    real c = Phi(muchat + delta[i, 1]);",
    "    real r = Phi(murhat + delta[i, 2]);",
    "    real u = Phi(muuhat + delta[i, 3]);",
    "    vector[4] theta;",
    "    theta[1] = c * r;",
    "    theta[2] = (1 - c) * u ^ 2;",
    "    theta[3] = (1 - c) * 2 * u * (1 - u);",
    "    theta[4] = c * (1 - r) + (1 - c) * (1 - u) ^ 2;",
    "    target += multinomial_lpmf(k[i] | theta);",
    "    target += multi_normal_cholesky_lpdf(to_vector(delta[i]) |",
    "                                         rep_vector(0, 3), L);",
    "  }",
    "}",
    sep = "\n")
}

# ---------------------------------------------------------------------
# MPT_1: the aggregate tree.
# ---------------------------------------------------------------------

bcm_mpt1_data <- function() {
  data.frame(k1 = 45L, k2 = 24L, k3 = 97L, k4 = 254L)
}

test_that("MPT_1 estimates the three process rates", {
  skip_unless_bcm("process-trees.R")
  d <- bcm_mpt1_data()
  fit <- frm(cbind(k1, k2, k3, k4) ~ 1, family = bcm_mpt_pairs(),
             data = d)
  rt <- bcm_mpt_rates(fit)
  expect_named(rt, c("c", "r", "u"))
  expect_true(all(rt > 0 & rt < 1))
  # the tree is saturated in three of its four cells, so the fitted
  # category probabilities reproduce the observed proportions
  th <- c(rt[["c"]] * rt[["r"]],
          (1 - rt[["c"]]) * rt[["u"]]^2,
          (1 - rt[["c"]]) * 2 * rt[["u"]] * (1 - rt[["u"]]),
          rt[["c"]] * (1 - rt[["r"]]) +
            (1 - rt[["c"]]) * (1 - rt[["u"]])^2)
  k <- as.numeric(d[1, ])
  expect_equal(unname(th), k / sum(k), tolerance = 1e-5)
})

test_that("MPT_1 matches its Stan program", {
  skip_unless_stan_identity()
  skip_unless_bcm("process-trees.R")
  d <- bcm_mpt1_data()
  fit <- frm(cbind(k1, k2, k3, k4) ~ 1, family = bcm_mpt_pairs(),
             data = d)
  stan_lp_check(
    bcm_mpt1_code(),
    data = list(k = as.integer(as.numeric(d[1, ]))),
    fit = fit,
    pars = function(f) as.list(bcm_mpt_rates(f)),
    const = 0)
})

# ---------------------------------------------------------------------
# MPT_2 through MPT_5: latent traits, correlated.
# ---------------------------------------------------------------------

bcm_mpt_hier_formula <- function() {
  bf(cbind(k1, k2, k3, k4) ~ 1 + (1 | p | id),
     r ~ 1 + (1 | p | id), u ~ 1 + (1 | p | id))
}

test_that("MPT_2 to MPT_5 are one hierarchical fit", {
  skip_unless_bcm("process-trees.R")
  d <- bcm_mpt_data()
  # The fit WARNS, and the warning is the finding rather than noise; the
  # next test is what says so. Suppressed here so this one can be about
  # the estimates.
  fit <- suppressWarnings(
    frm(bcm_mpt_hier_formula(), family = bcm_mpt_pairs(), data = d))
  vc <- unname(VarCorr(fit)[[1L]])
  expect_equal(dim(vc), c(3L, 3L))
  rho <- stats::cov2cor(vc)
  # three latent traits, so three correlations, and every one of them is
  # a number the case study plots
  expect_true(all(abs(rho[upper.tri(rho)]) <= 1))
  # the group means are the book's muc, mur, muu after the probit
  mu <- pnorm(c(unname(fixef(fit)[["c"]]), unname(fixef(fit)[["r"]]),
                unname(fixef(fit)[["u"]])))
  expect_true(all(mu > 0 & mu < 1))
  expect_equal(nrow(ranef(fit)[[1L]]), 21L)
  # the random effects are POOLED and not saturated: three latent traits
  # per subject against three free counts per subject would fit every
  # subject exactly, and these do not
  dp <- eval_dpars(fit)[[single_response(fit, "an MPT fit")$resp_name]]
  cc <- as.numeric(dp[["c"]])
  rr <- as.numeric(dp[["r"]])
  uu <- as.numeric(dp[["u"]])
  th <- cbind(cc * rr, (1 - cc) * uu^2, (1 - cc) * 2 * uu * (1 - uu),
              cc * (1 - rr) + (1 - cc) * (1 - uu)^2)
  k <- bcm_mpt_counts()
  expect_gt(max(abs(th - k / rowSums(k))), 0.05)
})

test_that("the MPT covariance block has no finite standard error", {
  skip_unless_bcm("process-trees.R")
  d <- bcm_mpt_data()
  # Twenty-one subjects with three free counts each, carrying an
  # unstructured three by three covariance, is six covariance parameters
  # the data barely identify. The POINT ESTIMATES are reproducible: five
  # optimizer settings (nlminb default and tightened, restarts, optim)
  # all return the same log-likelihood and the same estimates. What is
  # not identified is their curvature, and diagnose() says so rather
  # than the fit pretending otherwise.
  #
  # The case study needs the same thing from the other side: its own
  # README records convergence trouble under Stan, and MPT_2's Wishart
  # and MPT_3's LKJ prior are what hold this block down.
  # The fit raises TWO warnings and both are part of the claim, so they
  # are collected rather than matched one at a time: expect_warning()
  # consumes the first and lets the second escape to the reporter.
  w <- character()
  fit <- withCallingHandlers(
    frm(bcm_mpt_hier_formula(), family = bcm_mpt_pairs(), data = d),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("false convergence", w)))
  expect_true(any(grepl("may not have converged", w)))
  # diagnose() prints its report even when its value is assigned, and
  # four lines of it in the check log say nothing a reader of this test
  # needs; the assertions below are the report.
  invisible(utils::capture.output(dg <- diagnose(fit)))
  expect_false(dg$pdHess)
  expect_lt(dg$min_cov_eigenvalue, 0)
  # every parameter without a finite standard error is a COVARIANCE
  # parameter; the three group means and the traits themselves are fine
  expect_gt(length(dg$bad_se), 0)
  expect_true(all(grepl("^theta", dg$bad_se)))
  # WHICH covariance coordinate carries the largest gradient is a
  # tie-break among coordinates that are all unidentified, so pinning it
  # is a bet on this machine's BLAS rather than a claim about the model.
  # It is theta_3 under every optimizer setting tried here; the claim
  # worth making is that the worst gradient is IN the covariance block.
  expect_true(dg$worst_grad %in% paste0("theta_", 1:6))
})

test_that("MPT_2 to MPT_5 match their Stan program", {
  skip_unless_stan_identity()
  skip_unless_bcm("process-trees.R")
  d <- bcm_mpt_data()
  # warns about the covariance block; the test above is what asserts it
  fit <- suppressWarnings(
    frm(bcm_mpt_hier_formula(), family = bcm_mpt_pairs(), data = d))
  stan_lp_check(
    bcm_mpt_hier_code(),
    data = list(S = nrow(d), k = bcm_mpt_counts(),
                L = frm_chol(fit)),
    fit = fit,
    pars = function(f) {
      list(muchat = unname(fixef(f)[["c"]]),
           murhat = unname(fixef(f)[["r"]]),
           muuhat = unname(fixef(f)[["u"]]),
           delta = frm_u(f))
    },
    inner = "delta",
    const = 0)
})
