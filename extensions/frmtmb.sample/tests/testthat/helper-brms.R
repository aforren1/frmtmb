# Helpers for the brms cross-validation suite (test-brms-agreement.R).
#
# Two tiers. Tier 1 goes through brms's data-generating functions
# (make_standata, brmsterms, get_prior), which never touch Stan, so it
# runs whenever brms is installed and NOT_CRAN is set. Tier 2 fits Stan
# models and is opt-in only.

# Structural tier: brms must be installed, but nothing compiles.
skip_unless_brms <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("brms")
}

# Numeric tier: Stan compilation, minutes per model. Opt in with
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")
skip_unless_brms_fit <- function() {
  skip_unless_brms()
  testthat::skip_if_not_installed("rstan")
  if (!identical(Sys.getenv("FRMTMB_BRMS_FIT_TESTS"), "true")) {
    testthat::skip("set FRMTMB_BRMS_FIT_TESTS=true to run brms fit tests")
  }
}

# brms chatters about mixture ordering and dpar defaults; the design
# objects are what the tests read.
brms_standata <- function(...) {
  suppressMessages(brms::make_standata(...))
}

# Design matrices are compared by VALUE: brms names the intercept
# "Intercept" (no parentheses) so that it survives Stan's identifier
# rules, and drops matrix dimnames in places.
expect_design_equal <- function(x, y, tol = 1e-10) {
  x <- unname(as.matrix(x))
  y <- unname(as.matrix(y))
  testthat::expect_identical(dim(x), dim(y))
  testthat::expect_lt(max(abs(x - y)), tol)
}

# Column-space equality, for bases that span the same model but use a
# different parameterization (mgcv's diagonal.penalty reparameterization).
col_span_proj <- function(M) {
  M <- as.matrix(M)
  q <- qr(M)
  Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  Q %*% t(Q)
}

expect_span_equal <- function(x, y, tol = 1e-8) {
  testthat::expect_lt(max(abs(col_span_proj(x) - col_span_proj(y))), tol)
}
