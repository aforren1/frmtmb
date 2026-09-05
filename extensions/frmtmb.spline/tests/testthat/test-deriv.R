## Derivatives and features, checked against functions whose derivative
## and whose peak are known in closed form.

sp_bell_fit <- function(n = 500, seed = 12) {
  set.seed(seed)
  d <- data.frame(t = sort(stats::runif(n)))
  # a bell-shaped velocity profile: the minimum-jerk speed profile
  d$v <- 30 * d$t^2 * (1 - d$t)^2
  d$y <- d$v + stats::rnorm(n, 0, 0.15)
  list(d = d,
       fit = frmtmb::frm(frmtmb::bf(y ~ s(t, k = 15)),
                         family = stats::gaussian(), data = d))
}

test_that("the central difference reproduces a known derivative", {
  o <- sp_bell_fit()
  g <- data.frame(t = seq(0.1, 0.9, length.out = 21))
  d1 <- frm_curve_deriv(o$fit, var = "t", order = 1, newdata = g,
                        simultaneous = FALSE)
  truth <- 60 * g$t * (1 - g$t) * (1 - 2 * g$t)
  # within two pointwise standard errors everywhere: this is a
  # statistical statement about a fitted curve, not a numerical one
  expect_true(all(abs(d1$.estimate - truth) < 2.5 * d1$.se))
  d2 <- frm_curve_deriv(o$fit, var = "t", order = 2, newdata = g,
                        simultaneous = FALSE)
  truth2 <- 60 * (1 - 6 * g$t + 6 * g$t^2)
  expect_true(all(abs(d2$.estimate - truth2) < 2.5 * d2$.se))
})

test_that("the default step size is per order and scales with the range", {
  o <- sp_bell_fit()
  g <- data.frame(t = seq(0.1, 0.9, length.out = 11))
  a <- frm_curve_deriv(o$fit, var = "t", order = 1, newdata = g,
                       simultaneous = FALSE)
  b <- frm_curve_deriv(o$fit, var = "t", order = 2, newdata = g,
                       simultaneous = FALSE)
  expect_equal(attr(a, "eps"), 1e-6 * diff(range(g$t)))
  expect_equal(attr(b, "eps"), 1e-4 * diff(range(g$t)))
  expect_gt(attr(b, "eps"), attr(a, "eps"))
})

test_that("the step size that gratia fixes at 1e-7 would ruin order 2 here", {
  o <- sp_bell_fit()
  g <- data.frame(t = seq(0.1, 0.9, length.out = 11))
  good <- frm_curve_deriv(o$fit, var = "t", order = 2, newdata = g,
                          simultaneous = FALSE)
  bad <- frm_curve_deriv(o$fit, var = "t", order = 2, newdata = g,
                         eps = 1e-7, simultaneous = FALSE)
  # Cancellation at 1e-7 costs an ABSOLUTE amount, about 2 here, which
  # does not shrink when the curve does. That is why the same step size
  # is merely untidy on this curve, whose second derivative spans 360,
  # and ruinous on a flatter one: test-gratia.R measures 1.76 on a curve
  # whose second derivative is of order 5, which is 35 percent.
  truth2 <- 60 * (1 - 6 * g$t + 6 * g$t^2)
  expect_lt(max(abs(good$.estimate - truth2)), 40)
  expect_gt(max(abs(bad$.estimate - good$.estimate)), 1)
})

test_that("a derivative reuses a curve's grid and predictor", {
  o <- sp_bell_fit()
  g <- data.frame(t = seq(0.1, 0.9, length.out = 15))
  cv <- frm_curve(o$fit, newdata = g, simultaneous = FALSE)
  a <- frm_curve_deriv(cv, var = "t", order = 1, simultaneous = FALSE)
  b <- frm_curve_deriv(o$fit, var = "t", order = 1, newdata = g,
                       simultaneous = FALSE)
  expect_equal(a$.estimate, b$.estimate)
  expect_equal(a$.se, b$.se)
})

test_that("derivative refusals name the argument", {
  o <- sp_bell_fit()
  g <- data.frame(t = seq(0.1, 0.9, length.out = 11))
  expect_error(frm_curve_deriv(o$fit, var = "t", order = 3, newdata = g),
               "must be 1")
  expect_error(frm_curve_deriv(o$fit, var = "nope", newdata = g),
               "not a column of the grid")
  expect_error(frm_curve_deriv(o$fit, var = "t"), "`newdata` is required")
  expect_error(frm_curve_deriv(o$fit, var = "t", newdata = g, eps = -1),
               "one positive finite number")
  gf <- g
  gf$f <- factor("a")
  expect_error(frm_curve_deriv(o$fit, var = "f", newdata = gf),
               "NUMERIC covariate")
})

test_that("the peak of a known bell profile is recovered with an interval", {
  o <- sp_bell_fit()
  g <- data.frame(t = seq(0.05, 0.95, length.out = 41))
  ft <- frm_curve_feature(o$fit, var = "t", type = "maximum", newdata = g)
  expect_s3_class(ft, "frmtmb_feature")
  expect_equal(nrow(ft), 1L)
  # 30 t^2 (1-t)^2 peaks at t = 0.5, height 30/16 = 1.875
  expect_lt(abs(ft$.estimate - 0.5), 2.5 * ft$.se)
  expect_lt(abs(ft$.value - 1.875), 2.5 * ft$.value_se)
  expect_true(ft$.lower_ci < 0.5 && ft$.upper_ci > 0.5)
  # the peak HEIGHT standard error is the pointwise one at t*, because
  # the location term of the delta method drops out at a stationary point
  cv <- frm_curve(o$fit, newdata = data.frame(t = ft$.estimate + c(0, 1e-9)),
                  simultaneous = FALSE)
  expect_equal(ft$.value_se, cv$.se[1L], tolerance = 1e-5)
})

test_that("a crossing of a level is located with an interval", {
  o <- sp_bell_fit()
  g <- data.frame(t = seq(0.05, 0.95, length.out = 41))
  fc <- frm_curve_feature(o$fit, var = "t", type = "crossing", at = 1,
                          newdata = g)
  expect_equal(nrow(fc), 2L)
  # 30 t^2 (1-t)^2 = 1 where t(1-t) = sqrt(1/30), so t = 0.24030 and
  # t = 0.75970. The tolerance is absolute rather than a multiple of
  # the standard error, because what separates the fitted crossing from
  # the true one here is the penalty's bias and not sampling noise: the
  # standard error is 0.0034 and it is a statement about the fitted
  # curve, which is what the delta method describes.
  expect_lt(abs(fc$.estimate[1L] - 0.24030), 0.03)
  expect_lt(abs(fc$.estimate[2L] - 0.75970), 0.03)
  expect_true(all(abs(fc$.value - 1) < 1e-6))
  expect_true(all(fc$.se > 0))
})

test_that("no root in the window is an answer, not an error", {
  o <- sp_bell_fit()
  g <- data.frame(t = seq(0.05, 0.95, length.out = 41))
  fc <- frm_curve_feature(o$fit, var = "t", type = "crossing", at = 99,
                          newdata = g)
  expect_equal(nrow(fc), 0L)
  expect_output(print(fc), "no crossing in this window")
  # a minimum of a curve that only ever rises then falls is not there
  fm <- frm_curve_feature(o$fit, var = "t", type = "minimum",
                          newdata = data.frame(t = seq(0.3, 0.7, length.out = 21)))
  expect_equal(nrow(fm), 0L)
})

test_that("maximum and minimum select different stationary points", {
  set.seed(21)
  n <- 600
  d <- data.frame(t = sort(stats::runif(n)))
  d$y <- sin(4 * pi * d$t) + stats::rnorm(n, 0, 0.25)
  fit <- frmtmb::frm(frmtmb::bf(y ~ s(t, k = 20)),
                     family = stats::gaussian(), data = d)
  g <- data.frame(t = seq(0.05, 0.95, length.out = 61))
  mx <- frm_curve_feature(fit, var = "t", type = "maximum", newdata = g)
  mn <- frm_curve_feature(fit, var = "t", type = "minimum", newdata = g)
  ex <- frm_curve_feature(fit, var = "t", type = "extremum", newdata = g)
  # sin(4 pi t) peaks at 0.125 and 0.625, troughs at 0.375 and 0.875
  expect_equal(nrow(mx), 2L)
  expect_equal(nrow(mn), 2L)
  expect_equal(nrow(ex), nrow(mx) + nrow(mn))
  expect_lt(abs(mx$.estimate[1L] - 0.125), 0.03)
  expect_lt(abs(mn$.estimate[1L] - 0.375), 0.03)
  expect_true(all(mx$.value > 0))
  expect_true(all(mn$.value < 0))
})

test_that("a feature needs a grid to scan", {
  o <- sp_bell_fit()
  expect_error(frm_curve_feature(o$fit, var = "t",
                                 newdata = data.frame(t = c(0.2, 0.8))),
               "at least three points")
  expect_error(frm_curve_feature(o$fit, var = "t", at = NA,
                                 newdata = data.frame(t = seq(0, 1, 0.1))),
               "one finite number")
})
