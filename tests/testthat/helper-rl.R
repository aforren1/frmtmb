# Helpers for the reinforcement-learning worked example
# (vignette("reinforcement-learning")).
#
# The family itself lives in inst/rl/rw-delta.R, so the vignette and
# these tests fit exactly the same code. A test helper could not be the
# shared source: a vignette has no one relative path that reaches
# tests/testthat/ under both `R CMD build` (working directory
# vignettes/) and pkgdown (working directory the package root). inst/ is
# the one directory an installed package and a source build resolve the
# same way.

rl_source <- system.file("rl", "rw-delta.R", package = "frmtmb")
if (nzchar(rl_source)) source(rl_source, local = TRUE)

skip_unless_rl <- function() {
  testthat::skip_if_not(nzchar(rl_source),
                        "inst/rl/rw-delta.R is not installed")
}

# The Stan identity tier, gated like the brms fit tier: outside R CMD
# check BOTH variables are needed to opt in.
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")
skip_unless_rl_stan <- function() {
  skip_unless_rl()
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rstan")
  if (!identical(Sys.getenv("FRMTMB_BRMS_FIT_TESTS"), "true")) {
    testthat::skip("set FRMTMB_BRMS_FIT_TESTS=true to run the Stan tier")
  }
}

# The same model as a Stan program, written from the published
# equations. hBayesDM's own source is GPL-3 and is not read or copied.
#
# The parameters are the ones frmtmb estimates, on the scales frmtmb
# estimates them on: the two fixed effects of the learning rate, the one
# fixed effect of the softmax temperature, and the subject deviations
# themselves rather than a standardized version of them. So the map from
# frmtmb's estimates to Stan's parameters is the identity and the
# comparison carries no Jacobian. The covariance factor rides along as
# DATA, which is what keeps `log_prob` free of a constrained parameter.
rl_stan_code <- function() {
  paste(
    "data {",
    "  int<lower=1> N; int<lower=1> S; int<lower=1> T; int<lower=1> K;",
    "  array[S, T] int<lower=1> idx;",
    "  matrix[S, T] mask;",
    "  array[N] int<lower=1> subj;",
    "  array[N] int<lower=0, upper=1> choice;",
    "  vector[N] pay1; vector[N] pay2;",
    "  matrix[N, K] Xa;",
    "  matrix[2, 2] L;",
    "}",
    "parameters {",
    "  vector[K] ba; real bb; matrix[S, 2] u;",
    "}",
    "model {",
    "  vector[N] eta_a = Xa * ba;",
    "  vector[N] eta_b = rep_vector(bb, N);",
    "  vector[S] q1 = rep_vector(0, S);",
    "  vector[S] q2 = rep_vector(0, S);",
    "  for (i in 1:N) {",
    "    eta_a[i] += u[subj[i], 1];",
    "    eta_b[i] += u[subj[i], 2];",
    "  }",
    "  {",
    "    vector[N] alpha = inv_logit(eta_a);",
    "    vector[N] beta = exp(eta_b);",
    "    for (t in 1:T) {",
    "      for (s in 1:S) {",
    "        if (mask[s, t] == 1) {",
    "          int i = idx[s, t];",
    "          real e = beta[i] * (q1[s] - q2[s]);",
    "          target += bernoulli_logit_lpmf(choice[i] | e);",
    "          if (choice[i] == 1)",
    "            q1[s] += alpha[i] * (pay1[i] - q1[s]);",
    "          else",
    "            q2[s] += alpha[i] * (pay2[i] - q2[s]);",
    "        }",
    "      }",
    "    }",
    "  }",
    "  for (s in 1:S)",
    "    target += multi_normal_cholesky_lpdf(to_vector(u[s]) |",
    "                                         rep_vector(0, 2), L);",
    "}",
    sep = "\n")
}

# The design matrix is rebuilt with model.matrix() rather than read off
# the frame, so a mistake in frmtmb's own assembly would show up as a
# mismatch instead of cancelling on both sides.
rl_stan_data <- function(fit, data) {
  blk <- frame_block_of(fit$frame, "choice")
  vc <- unname(VarCorr(fit)[[1L]])
  list(N = nrow(data), S = blk[["n_subj"]], T = blk[["n_trial"]],
       K = 2L, idx = blk[["idx"]], mask = blk[["mask"]],
       subj = as.integer(factor(data$id)),
       choice = as.integer(data$choice),
       pay1 = as.numeric(data$pay1), pay2 = as.numeric(data$pay2),
       Xa = unname(stats::model.matrix(~ condition, data)),
       L = t(chol(vc)))
}

# frmtmb's full parameter vector, split the way the Stan program
# declares it. `b` is level-major (both of a subject's coefficients
# contiguous), which is what makes the reshape a plain matrix().
rl_stan_pars <- function(par) {
  list(ba = unname(par[names(par) == "beta"]),
       bb = unname(par[names(par) == "betad"]),
       u = t(matrix(unname(par[names(par) == "b"]), nrow = 2L)))
}

# The two checks, at one point.
#
#   A  Stan's log_prob at frmtmb's estimates equals frmtmb's joint log
#      density there, up to a constant this asserts is zero.
#   B  the gradient with respect to the subject effects vanishes, which
#      is what catches a wrong map: the outer gradient of the joint is
#      NOT zero at the marginal optimum, so only the inner block can be
#      asserted on.
rl_lp_check <- function(fit, data, par = NULL, tol = 1e-6,
                        tol_grad = 1e-4, check_grad = TRUE) {
  mod <- brms_stan_model(rl_stan_code())
  sdat <- rl_stan_data(fit, data)
  sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))
  par <- par %||% fit$obj$env$last.par.best
  pars <- rl_stan_pars(par)
  up <- rstan::unconstrain_pars(sf, pars)
  lp <- rstan::log_prob(sf, up, adjust_transform = FALSE, gradient = FALSE)
  ours <- -fit$obj$env$f(par)
  testthat::expect_lt(abs(lp - ours), tol * max(1, abs(ours)))
  if (check_grad) {
    gr <- rstan::grad_log_prob(sf, up, adjust_transform = FALSE)
    bump <- pars
    bump$u <- bump$u + 1
    inner <- which(rstan::unconstrain_pars(sf, bump) != up)
    testthat::expect_gt(length(inner), 0)
    testthat::expect_lt(max(abs(gr[inner])), tol_grad)
  }
  invisible(list(lp = lp, ours = ours, const = lp - ours))
}
