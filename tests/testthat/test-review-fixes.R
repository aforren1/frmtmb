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

test_that("frm_sample(priors=) works on a fixed-effects-only GLM", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  set.seed(404)
  d <- data.frame(x = rnorm(80))
  d$y <- rpois(80, exp(0.4 + 0.5 * d$x))
  fit <- frm(bf(y ~ x) + poisson(), data = d)
  # the $b partial match used to pass random = "b" for a template that
  # only holds beta, so this errored in MakeADFun
  # short-chain R-hat/ESS warnings are not what this test checks
  ds <- suppressWarnings(
    frm_sample(fit, chains = 1, iter = 800, refresh = 0,
               priors = list(beta = prior_normal(0, 5))))
  m <- as.matrix(ds)
  expect_true("x" %in% colnames(m))
  expect_lt(abs(mean(m[, "x"]) - unname(fit$estimates$beta[["x"]])), 0.2)
})

test_that("frm_sample(laplace = TRUE) runs and labels outer draws", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  set.seed(405)
  d <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
  d$y <- rnorm(80, 1 + 0.5 * d$x + rnorm(8, 0, 0.5)[d$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d)
  # rstan warns about lp__ under laplace; not what this test checks
  ds <- suppressWarnings(
    frm_sample(fit, chains = 1, iter = 400, refresh = 0,
               laplace = TRUE))
  m <- as.matrix(ds)
  # no b columns are sampled, and theta keeps its own label instead of
  # being misattributed as b[1]
  expect_false(any(grepl("^b\\[", colnames(m))))
  expect_true(any(grepl("theta", colnames(m))))
  expect_true("x" %in% colnames(m))
  expect_lt(abs(mean(m[, "x"]) - unname(fit$estimates$beta[["x"]])), 0.25)
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

test_that("mode_inits anchors chain 1 and jitters the rest", {
  mode <- c(a = 1, b = -2, c = 0.5)
  ii <- mode_inits(mode, chains = 4, jitter = 0.25)
  expect_length(ii, 4L)
  expect_identical(ii[[1]], as.numeric(mode))
  for (k in 2:4) {
    expect_false(identical(ii[[k]], as.numeric(mode)))
    expect_lt(max(abs(ii[[k]] - as.numeric(mode))), 2)  # modest jitter
  }
  # jitter = 0 restores identical mode starts
  i0 <- mode_inits(mode, chains = 3, jitter = 0)
  expect_true(all(vapply(i0, identical, TRUE, as.numeric(mode))))
})

test_that("frm_sample runs multiple chains with jittered mode inits", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  set.seed(409)
  d <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  d$y <- rnorm(60, 1 + 0.5 * d$x + rnorm(6, 0, 0.5)[d$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d)
  ds <- suppressWarnings(
    frm_sample(fit, chains = 2, iter = 300, refresh = 0))
  a <- rstan::extract(ds$stanfit, permuted = FALSE)
  expect_identical(dim(a)[2], 2L)
  # a boundary-ish mode triggers the singular-init warning; the
  # crippled 10-iteration run may warn on its own, so collect all
  fit2 <- fit
  fit2$estimates$theta <- c(-9)
  w <- testthat::capture_warnings(
    try(frm_sample(fit2, chains = 1, iter = 10, refresh = 0),
        silent = TRUE))
  expect_true(any(grepl("extreme covariance parameter", w)))
})
