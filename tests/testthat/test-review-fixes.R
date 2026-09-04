# Heavy reference-validation file: full runs happen locally and in CI
# (NOT_CRAN=true). CRAN only needs to see the package work, not the
# validation program re-derived, and this file is a main contributor
# to CRAN-condition check time.
skip_on_cran()

# v0.21 code-review fixes: cross-feature guards and the outer-parameter
# map (mi()/profile alignment across confint, bounds, and sampling).

test_that("rescor refuses cens()/trunc()/se() addition terms", {
  set.seed(401)
  d <- data.frame(y1 = rnorm(40), y2 = rnorm(40), x = rnorm(40))
  d$c1 <- 0L
  d$sdy <- 0.5
  expect_error(
    frm(mvbf(bf(y1 | cens(c1) ~ x) + gaussian(),
             bf(y2 ~ x) + gaussian(), rescor = TRUE), data = d),
    "cannot be combined with rescor")
  expect_error(
    frm(mvbf(bf(y1 | se(sdy) ~ x) + gaussian(),
             bf(y2 ~ x) + gaussian(), rescor = TRUE), data = d),
    "cannot be combined with rescor")
  expect_error(
    frm(mvbf(bf(y1 | trunc(lb = -5) ~ x) + gaussian(),
             bf(y2 ~ x) + gaussian(), rescor = TRUE), data = d),
    "cannot be combined with rescor")
})

test_that("bounds align with the outer parameters under profile = TRUE", {
  set.seed(402)
  d <- data.frame(x = rnorm(120), g = factor(rep(1:12, 10)))
  d$y <- rnorm(120, 1 + 0.7 * d$x + rnorm(12, 0, 0.5)[d$g], 1)
  # beta is inner under profile, so a bound on a beta coefficient must
  # error loudly instead of silently pinning a theta component
  expect_error(
    frm(bf(y ~ x + (1 | g)) + gaussian(), data = d,
        lower = c(x = 2),
        control = frmtmb_control(profile = TRUE)),
    "Unknown parameter")
  # a theta bound reaches the right slot
  fitb <- suppressWarnings(
    frm(bf(y ~ x + (1 | g)) + gaussian(), data = d,
        lower = c(theta_1 = 1),
        control = frmtmb_control(profile = TRUE)))
  expect_gte(fitb$estimates$theta[1], 1 - 1e-8)
  # and the slope stays at its unbounded estimate
  expect_lt(abs(unname(fitb$estimates$beta[["x"]]) - 0.7), 0.3)
})

test_that("confint() works on mi() fits (miss is inner, not outer)", {
  set.seed(403)
  n <- 120
  z <- rnorm(n)
  x <- rnorm(n, 0.5 + 0.8 * z, 0.7)
  y <- rnorm(n, 1 + 0.6 * x + 0.3 * z, 0.9)
  x[sample(n, 30)] <- NA
  dd <- data.frame(y = y, x = x, z = z)
  fit <- frm(bf(y ~ mi(x) + z) + gaussian() +
               bf(x | mi() ~ z) + gaussian(), data = dd)
  nm <- outer_par_names(fit)
  expect_identical(length(nm), length(fit$opt$par))
  expect_false(any(grepl("^miss", nm)))
  ci <- confint(fit)
  expect_identical(nrow(ci), length(fit$opt$par))
  expect_true(all(ci[, "lwr"] < ci[, "upr"]))
  # vcov(full = TRUE) labels ride the same map
  Vf <- vcov(fit, full = TRUE)
  expect_identical(rownames(Vf), nm)
})


test_that("influence machinery works with a constant dpar", {
  set.seed(406)
  d <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  d$y <- rnorm(60, 1 + 0.5 * d$x + rnorm(6, 0, 0.4)[d$g], 1)
  fit <- suppressWarnings(
    frm(bf(y ~ x + (1 | g), sigma = 1) + gaussian(), data = d))
  infl <- suppressWarnings(influence(fit, groups = "g"))
  # influence columns now align with vcov (mapped constants excluded)
  expect_identical(colnames(infl$fixed), rownames(vcov(fit)))
  cd <- cooks.distance(infl)
  expect_length(cd, 6L)
  expect_true(all(is.finite(cd)))
})

test_that("a covariate literally named sigma stays visible", {
  set.seed(407)
  d <- data.frame(sigma = rnorm(50))
  d$y <- 2 + 0.8 * d$sigma + rnorm(50)
  fit <- frm(bf(y ~ sigma) + gaussian(), data = d)
  h <- hypothesis(fit, "sigma = 0")
  expect_equal(h$estimate[1], unname(fit$estimates$beta[["sigma"]]),
               tolerance = 1e-8)
})

test_that("NA in an RE-only design variable propagates to predictions", {
  set.seed(408)
  d <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
  d$y <- rnorm(80, 1 + (0.5 + rnorm(8, 0, 0.3)[d$g]) * d$x, 1)
  fit <- frm(bf(y ~ 1 + (0 + x | g)) + gaussian(), data = d)
  nd <- data.frame(x = c(1, NA), g = factor(c(1, 2), levels = levels(d$g)))
  p <- predict(fit, newdata = nd)
  expect_true(is.finite(p[1]))
  expect_true(is.na(p[2]))
})

