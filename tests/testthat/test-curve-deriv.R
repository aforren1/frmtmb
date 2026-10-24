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
  # Cancellation at 1e-7 costs an ABSOLUTE amount that does not shrink
  # when the curve does: about 2 on Windows and 0.5 on Linux, since it
  # is machine rounding divided by 1e-14 and the rounding differs by
  # platform. So the claim is relative, not a number: the fixed step is
  # thousands of times noisier than the measured one, whose own
  # numerical error is read off by halving it. That is why the same
  # step size is merely untidy on this curve, whose second derivative
  # spans 360, and ruinous on a flatter one: test-gratia.R measures 1.76
  # on a curve whose second derivative is of order 5, 35 percent.
  truth2 <- 60 * (1 - 6 * g$t + 6 * g$t^2)
  expect_lt(max(abs(good$.estimate - truth2)), 40)
  half <- frm_curve_deriv(o$fit, var = "t", order = 2, newdata = g,
                          eps = attr(good, "eps") / 2, simultaneous = FALSE)
  good_err <- max(abs(good$.estimate - half$.estimate))
  bad_err <- max(abs(bad$.estimate - good$.estimate))
  expect_gt(bad_err, 1000 * good_err)
  expect_gt(bad_err, 0.05)
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

test_that("an exact gp()'s curvature past its data keeps its kriging part", {
  skip_on_cran()
  # dev/reviews/2026-10-05-gpby.md, B1: the stencil applied to the extra
  # covariance divided its rounding by eps^4 at order 2, and the .se past
  # the data came out 0, 0, 0, 0.875. The extra part is now the kernel's
  # closed form, written out here by hand from the fitted theta
  set.seed(5)
  x <- sort(stats::runif(50, 0, 4))
  dg <- data.frame(x = x, y = sin(1.3 * x) + stats::rnorm(50, 0, 0.15))
  fg <- frmtmb::frm(frmtmb::bf(y ~ gp(x)), family = stats::gaussian(),
                    data = dg)
  th <- fg$estimates$theta
  s2 <- exp(2 * th[1])
  l <- exp(th[2])
  pos <- fg$frame$linpreds[["y.mu"]]$gps[[1]]$positions[, 1]
  K <- exp(-outer(pos, pos, "-")^2 / (2 * l^2)) + diag(1e-6, length(pos))
  xs <- c(3.0, 4.2, 4.8, 5.5)
  r <- outer(xs, pos, "-")
  k0 <- exp(-r^2 / (2 * l^2))
  k1 <- -r / l^2 * k0
  k2 <- (r^2 / l^4 - 1 / l^2) * k0
  v1 <- s2 * (1 / l^2 - rowSums((k1 %*% solve(K)) * k1))
  v2 <- s2 * (3 / l^4 - rowSums((k2 %*% solve(K)) * k2))
  nd <- data.frame(x = xs)
  # the curve's own standard error carries it, does not depend on the
  # design step, and is at least the extra part
  d2 <- frm_curve_deriv(fg, var = "x", order = 2, newdata = nd,
                        simultaneous = FALSE)
  d2b <- frm_curve_deriv(fg, var = "x", order = 2, newdata = nd,
                         simultaneous = FALSE, eps = 4 * attr(d2, "eps"))
  expect_true(all(d2$.se^2 >= v2))
  expect_lt(max(abs(d2$.se / d2b$.se - 1)), 1e-4)
  d1 <- frm_curve_deriv(fg, var = "x", order = 1, newdata = nd,
                        simultaneous = FALSE)
  expect_true(all(d1$.se^2 >= v1))
  # and a simultaneous band over it is an ordinary one
  g2 <- frm_curve_deriv(fg, var = "x", order = 2,
                        newdata = data.frame(x = seq(3.5, 5.5,
                                                     length.out = 21)),
                        nsim = 2000, seed = 1)
  expect_true(all(g2$.se > 0))
  expect_lt(g2$.crit_sim[1], 2 * stats::qnorm(0.975))
  # and the extra part is the closed form. Both subtract the kriging
  # quadratic form from the derivative's prior variance, so inside the
  # data (x = 3, v1 near 1e-6) they agree to rounding of the prior
  # variance, not of v1 itself
  e1 <- diag(frmtmb::frm_extra_cov_deriv(fg, nd, var = "x", order = 1))
  e2 <- diag(frmtmb::frm_extra_cov_deriv(fg, nd, var = "x", order = 2))
  expect_lt(max(abs(e1 - v1)), 1e-8 * s2 / l^2)
  expect_lt(max(abs(e2 - v2)), 1e-8 * 3 * s2 / l^4)
})

test_that("an unseen level's draw carries no curvature a line does not", {
  skip_on_cran()
  # B1 again: y ~ t + (1 + t | g) is linear in t, so the new level's
  # draw has no second derivative and the order-2 standard error is the
  # coefficients' alone. The stencil on the extra covariance gave .se
  # from 0 to 1.495
  set.seed(23)
  ng <- 25
  d <- data.frame(g = factor(rep(seq_len(ng), each = 12)),
                  t = rep(seq(0, 1, length.out = 12), ng))
  u0 <- stats::rnorm(ng, 0, 0.5)
  u1 <- stats::rnorm(ng, 0, 0.3)
  d$y <- -1 + 2 * d$t + u0[d$g] + u1[d$g] * d$t +
    stats::rnorm(nrow(d), 0, 0.2)
  fit <- frmtmb::frm(frmtmb::bf(y ~ t + (1 + t | g)),
                     family = stats::gaussian(), data = d)
  tt <- seq(0, 1, length.out = 21)
  nd <- data.frame(t = tt, g = factor("n1", levels = c(levels(d$g), "n1")))
  # order 1 at the level: the slope's own variance, from VarCorr()
  d1 <- frm_curve_deriv(fit, var = "t", order = 1, newdata = nd,
                        re_formula = NULL, allow_new_levels = TRUE,
                        simultaneous = FALSE)
  Sg <- frmtmb::VarCorr(fit)$g$cov[, "Estimate", ]
  ref <- sqrt(stats::vcov(fit)["t", "t"] + Sg[2, 2])
  expect_lt(max(abs(d1$.se / ref - 1)), 1e-6)
  # order 2: coefficients and draw alike are a line in t, so the whole
  # .se is rounding, far below the slope's
  d2 <- frm_curve_deriv(fit, var = "t", order = 2, newdata = nd,
                        re_formula = NULL, allow_new_levels = TRUE,
                        simultaneous = FALSE)
  expect_lt(max(d2$.se), 1e-6 * min(d1$.se))
  E2 <- frmtmb::frm_extra_cov_deriv(fit, nd, var = "t", order = 2,
                                    re_formula = NULL,
                                    allow_new_levels = TRUE)
  expect_lt(max(abs(E2)), 1e-12 * min(d1$.se^2))
})
