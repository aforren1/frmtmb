# Heavy reference-validation file: full runs happen locally and in CI
# (NOT_CRAN=true). CRAN only needs to see the package work, not the
# validation program re-derived.
skip_on_cran()

# Confirmed defects from the v0.28 wave review.
#
# The one that changed answers is A1: spde() derived its mesh node
# ordering from the grouping variable's labels, which the finite-element
# matrices know nothing about, so integer node ids were read
# lexicographically and the field was fitted against a permuted
# precision. The rest are diagnosis and validation surfaces that either
# over-claimed, under-checked, or repeated themselves.

# The car() and spde() items of this review (A1, A5, A7 and the
# call-valued gr) moved to frmtmb.spatial with the two terms
# (extensions/frmtmb.spatial/tests/testthat/test-car-v28.R and
# test-spde.R).


# Run expr with mc.cores set, collecting every warning it raises: an
# rstan run emits its own diagnostics, and this file must not leave
# stray warnings in the suite summary.
v28_warnings_under_cores <- function(cores, expr) {
  old <- options(mc.cores = cores)
  on.exit(options(old), add = TRUE)
  msgs <- character(0)
  val <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = msgs)
}

# --------------------------------------------------- A1 spde node order







# ------------------------------------------------- A2 call-valued gr



# ------------------------------------------- A3 cores vs options(mc.cores)


# --------------------------------------------- A4 fit_error_context

test_that("only numerical failures get the undefined-likelihood cause", {
  set.seed(8)
  dd <- data.frame(x = stats::rnorm(40))
  dd$y <- 1 + dd$x + stats::rnorm(40)
  # an unknown optimizer is a typo, not an unbounded likelihood
  e1 <- tryCatch(frm(bf(y ~ x) + gaussian(), data = dd,
                     control = frmtmb_control(optimizer = "nope")),
                 error = function(e) conditionMessage(e))
  expect_match(e1, "Unknown optimizer 'nope'")
  expect_false(grepl("likelihood was undefined", e1))
  # neither is a custom optimizer that broke its return contract
  e2 <- tryCatch(frm(bf(y ~ x) + gaussian(), data = dd,
                     control = frmtmb_control(
                       optimizer = function(par, fn, gr, lower, upper,
                                            control) list(par = par))),
                 error = function(e) conditionMessage(e))
  expect_match(e2, "must return par, objective, and convergence")
  expect_false(grepl("likelihood was undefined", e2))
  # the model is still named, so the user knows which fit died
  expect_match(e1, "raised while fitting")

  # a genuinely undefined objective keeps the diagnosis and the remedies
  spec <- frm(bf(y ~ x) + gaussian(), data = dd, dry_run = "spec")
  e3 <- tryCatch(
    frmtmb:::fit_error_context(spec, NULL, FALSE, frmtmb_control(), NULL,
                               NULL, stop("NA/NaN function evaluation")),
    error = function(e) conditionMessage(e))
  expect_match(e3, "likelihood was undefined or unbounded")
  expect_match(e3, "NA/NaN function evaluation")
  expect_match(e3, "frmtmb_control\\(optimizer")
})

# ------------------------------------------------ A5 escar at rho -> 1


# ------------------------------------- A6 extreme-theta reads log sds only


test_that("a collapsed log sd is still reported, by name", {
  set.seed(12)
  dd <- data.frame(g = factor(rep(1:8, 6)), x = stats::rnorm(48))
  dd$y <- 1 + 0.5 * dd$x + stats::rnorm(48, 0, 0.7)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  # force the collapsed mode the heuristic is about
  fit$estimates$theta[1] <- -12
  dg <- diagnose(fit, quiet = TRUE)
  expect_length(dg$extreme_theta, 1L)
  # the message names the confint() row, not a bare theta index
  expect_match(names(dg$extreme_theta), "1 | g", fixed = TRUE)
  expect_match(names(dg$extreme_theta), "(Intercept)", fixed = TRUE)
  out <- utils::capture.output(diagnose(fit))
  expect_true(any(grepl("Extreme covariance parameters", out)))
  expect_true(any(grepl("1 | g", out, fixed = TRUE)))
})

# --------------------------------------------- A7 car adjacency validation



# --------------------------------------- A8 one non-finite-cov warning

test_that("a degenerate covariance warns once per fit", {
  cache <- new.env(parent = emptyenv())
  expect_warning(frmtmb:::warn_nonfinite_cov(cache), "not finite")
  expect_silent(frmtmb:::warn_nonfinite_cov(cache))
  # a fresh fit (or a refit, which replaces the cache) says it again
  expect_warning(frmtmb:::warn_nonfinite_cov(new.env(parent = emptyenv())),
                 "not finite")
  # no cache, no memory: the helper still warns when called bare
  expect_warning(frmtmb:::warn_nonfinite_cov(), "not finite")
})

test_that("repeated vcov() calls on a degenerate fit warn once", {
  set.seed(13)
  dd <- data.frame(x = stats::rnorm(30))
  dd$y <- 1 + dd$x + stats::rnorm(30)
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  # a NaN covariance is what sdreport returns on a singular Hessian;
  # inject it so the degradation path runs without a pathological model
  sdr <- sdr_of(fit)
  sdr$cov.fixed[] <- NaN
  fit$cache$sdr <- sdr
  expect_warning(vcov(fit), "not finite")
  expect_silent(vcov(fit))
  expect_silent(summary(fit))
})
