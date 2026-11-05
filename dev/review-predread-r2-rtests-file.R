# The tests that go with release fixes F1 to F4; see review-predread-r2-rtests.R.
test_that("a me() fit's mean answers on more rows than a key can hold", {
  set.seed(7981)
  n <- 200
  d <- data.frame(lat = rnorm(n), sx = 0.4, z = rnorm(n))
  d$xo <- d$lat + rnorm(n, 0, 0.4)
  d$cnt <- rpois(n, exp(0.3 + 0.5 * d$lat - 0.3 * d$z))
  fit <- frm(bf(cnt ~ me(xo, sx) + z), family = poisson(), data = d)
  big <- data.frame(xo = seq(-2, 2, length.out = 800L), sx = 0.4, z = 0)
  f <- fitted(fit, newdata = big)
  expect_equal(dim(f), c(800L, 4L))
  # the same row reads the same whatever else the grid holds
  few <- fitted(fit, newdata = big[c(1L, 400L, 800L), ])
  expect_equal(f[c(1L, 400L, 800L), "Estimate"], few[, "Estimate"], tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(f[c(1L, 400L, 800L), "Q97.5"], few[, "Q97.5"], tolerance = 1e-10,
               ignore_attr = TRUE)
})
test_that("frm_curve() reads a row of NA as NA at its default and on a contrast", {
  set.seed(7982)
  n <- 200
  d <- data.frame(x = rnorm(n), z = rnorm(n), sx = 0.4)
  d$xo <- d$x + rnorm(n, 0, 0.4)
  d$y <- 1 + d$x + rnorm(n)
  d$cat <- as.integer(cut(d$x + rlogis(n), c(-Inf, -0.5, 0.5, Inf), labels = FALSE))
  pl <- frm(bf(y ~ x + z), family = gaussian(), data = d)
  nd <- data.frame(x = c(1, NA, 0), z = 0)
  cv <- frm_curve(pl, newdata = nd, seed = 1)
  expect_equal(is.na(cv$.se), c(FALSE, TRUE, FALSE))
  # the band of the other rows is the band without that row
  ok <- frm_curve(pl, newdata = nd[-2L, ], seed = 1)
  expect_equal(cv$.se[-2L], ok$.se, tolerance = 1e-10)
  mc <- frm(bf(cat ~ me(xo, sx) + z), family = categorical(), data = d)
  nm <- data.frame(xo = c(1, NA, 0), sx = 0.4, z = 0)
  ct <- suppressWarnings(frm_curve(mc, newdata = nm, contrast = transform(nm, z = 1),
                                   simultaneous = FALSE))
  expect_true(all(is.na(ct$.se[c(2L, 5L)])))
  expect_true(all(is.finite(ct$.se[c(1L, 3L, 4L, 6L)])))
})
test_that("a middle category's Q on a new me() row is the draw's own quantile", {
  set.seed(7971)
  n <- 400
  d <- data.frame(lat = rnorm(n), sx = 0.4)
  d$xo <- d$lat + rnorm(n, 0, 0.4)
  d$ord <- as.integer(cut(0.4 + 0.7 * d$lat + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf),
                          labels = FALSE))
  fit <- frm(bf(ord ~ me(xo, sx)), family = cumulative(), data = d)
  f <- fitted(fit, newdata = data.frame(xo = 1, sx = 1.2))
  cf <- fixef(fit)[, "Estimate"]
  b <- cf[grep("me", names(cf))]
  th <- cf[grep("Intercept", names(cf))]
  set.seed(1)
  l <- rnorm(1e5, 1, 1.2)
  # the draw alone, with the coefficients at their estimates
  p <- cbind(stats::plogis(th[1L] - b * l),
             stats::plogis(th[2L] - b * l) - stats::plogis(th[1L] - b * l),
             stats::plogis(th[3L] - b * l) - stats::plogis(th[2L] - b * l),
             1 - stats::plogis(th[3L] - b * l))
  q <- apply(p, 2L, stats::quantile, c(0.025, 0.975), names = FALSE)
  for (k in 2:3) {
    expect_lt(abs(f[1L, "Q97.5", k] - q[2L, k]), 0.02)
    expect_lt(abs(f[1L, "Q2.5", k] - q[1L, k]), 0.02)
  }
  # the mean stays inside its interval, and the order holds
  for (k in 1:4) {
    expect_lt(f[1L, "Q2.5", k], f[1L, "Q97.5", k])
  }
})
test_that("the warning says the standard error can be off in either direction", {
  set.seed(7983)
  d <- data.frame(x = rnorm(120), g = factor(rep(1:12, 10)))
  d$y <- 1 + d$x + rep(rnorm(12), 10) + rnorm(120)
  fit <- suppressWarnings(frm(y ~ x + (1 | g), data = d))
  expect_warning(marginaleffects::get_vcov(fit), "too wide or too narrow")
  expect_warning(marginaleffects::get_vcov(fit), "standard errors from the population-level")
})
