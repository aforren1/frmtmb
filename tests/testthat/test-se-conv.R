# The standard-error check on fits that carry a convergence warning
# (lane illse, item c). It used to skip them, so on a ridge whether a
# fit's standard errors were judged depended on which side of code 0
# rounding put nlminb (lane ci691). Now every fit is checked, and the
# warning says whether the stop is the maximum: a Newton step from it
# gains at most grad_tol, or more.

conv_conds <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w)
}

conv_phrase <- "Optimizer did not report convergence"

# nl `a + b`: the two enter only through their sum, an exact alias the
# check names on a converged fit
conv_alias_data <- function(seed = 1) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(100))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(100, 0, 0.5)
  d
}
conv_alias_bf <- bf(y ~ a + b + c * x, a + b + c ~ 1, nl = TRUE)

test_that("a stop at the maximum nlminb did not call converged is judged", {
  d <- conv_alias_data()
  ref <- conv_conds(frm(conv_alias_bf, data = d))
  expect_equal(ref$value$opt$convergence, 0L)
  lost0 <- sdr_of(ref$value)$se_lost
  expect_true(length(lost0) > 0L)
  # the same model started at its own maximum, with no iteration
  # allowed: nlminb reports the iteration limit at the maximum
  st <- ref$value$estimates[c("beta", "betad")]
  r <- conv_conds(frm(conv_alias_bf, data = d, start = st,
                      control = frmtmb_control(optCtrl = list(
                        iter.max = 0))))
  expect_false(r$value$opt$convergence == 0L)
  w <- grep(conv_phrase, r$warnings, fixed = TRUE, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "so the estimates are at the maximum", fixed = TRUE)
  # judged as the converged fit is: the same lost set, and the same
  # message for it (the nonlinear alias's own), beside the convergence
  # warning and no second one in other words
  expect_identical(sdr_of(r$value)$se_lost, lost0)
  expect_identical(setdiff(r$warnings, w), ref$warnings)
  expect_length(grep("Large maximum absolute gradient", r$warnings,
                     fixed = TRUE), 0L)
})

test_that("a stop short of the maximum says so and keeps its SEs", {
  # A curvature verdict at a point short of the maximum is not one about
  # the data, so such a fit keeps sdreport()'s standard errors, and the
  # warning says they are those of the stop
  # A poisson fit from an explicit start far from the maximum, so that
  # one iteration stops short whatever start frm() would choose itself
  # (from a least-squares start a gaussian y ~ x is at the maximum
  # after one iteration)
  set.seed(2)
  d <- data.frame(x = stats::rnorm(200))
  d$y <- stats::rpois(200, exp(1 + 0.5 * d$x))
  r <- conv_conds(frm(y ~ x, data = d, family = poisson(),
                      start = list(beta = c(-3, 3)),
                      control = frmtmb_control(optCtrl = list(
                        iter.max = 1), restarts = 0)))
  expect_false(r$value$opt$convergence == 0L)
  w <- grep(conv_phrase, r$warnings, fixed = TRUE, value = TRUE)
  expect_length(w, 1L)
  expect_false(grepl("so the estimates are at the maximum", w,
                     fixed = TRUE))
  expect_true(grepl("short of the maximum", w, fixed = TRUE) ||
                grepl("curvature at the stop is unusable", w, fixed = TRUE))
  # sdreport()'s covariance, not repaired or judged
  s <- sdr_of(r$value)
  expect_null(s$cov_fixed_prop)
  expect_length(s$se_lost, 0L)
  expect_true(all(is.finite(sqrt(diag(s$cov.fixed)))))
})

test_that("a converged fit without a convergence warning is unchanged", {
  # the guard's absent case: the clause belongs to the convergence
  # warning only
  d <- conv_alias_data()
  fit <- frm(y ~ x, data = d)
  r <- conv_conds(vcov(fit))
  expect_length(r$warnings, 0L)
  expect_equal(fit$opt$convergence, 0L)
})

test_that("an sd held near zero short of its maximum is not at the maximum", {
  # The sd's row is empty and its gradient is 5.6e-10, so a Newton step
  # cannot see what a larger sd gains. The headroom asks the likelihood
  # (se_sd_gain_up()); without that, this stop 92 below the maximum was
  # called the maximum (review of lane illse, round 2, R2-1)
  set.seed(1)
  n <- 300
  d <- data.frame(x = stats::rnorm(n),
                  g = factor(rep(1:20, length.out = n)))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(20, 0, 0.5)[d$g] +
    stats::rnorm(n, 0, 0.4)
  r <- conv_conds(frm(y ~ x + (1 | g), data = d,
                      start = list(theta = -15),
                      control = frmtmb_control(optCtrl = list(
                        iter.max = 3, eval.max = 1000))))
  expect_false(r$value$opt$convergence == 0L)
  w <- grep(conv_phrase, r$warnings, fixed = TRUE, value = TRUE)
  expect_length(w, 1L)
  expect_false(grepl("so the estimates are at the maximum", w,
                     fixed = TRUE))
  expect_match(w, "short of the maximum", fixed = TRUE)
  expect_length(grep("The fit may not have converged", r$warnings,
                     fixed = TRUE), 1L)
})
