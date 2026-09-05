## Against mgcv and gratia, which is where the claims about this
## package's method are settled.
##
## Two different comparisons live here and they must not be confused.
## The MODEL comparison asks whether frmtmb's ML fit is mgcv's, which
## the core case-studies vignette already establishes and which is
## re-checked here because everything else rests on it. The ALGORITHM
## comparison asks whether this package's max-deviation simulation is
## gratia's, and it is run on gratia's own inputs so that a difference
## in the fitted covariance cannot be mistaken for a difference in the
## method.

sp_gratia_data <- function(n = 300, seed = 7) {
  set.seed(seed)
  d <- data.frame(x = sort(stats::runif(n)))
  d$y <- 2 * sin(pi * d$x) + 0.6 * d$x + stats::rnorm(n, 0, 0.4)
  d
}

test_that("frmtmb reproduces the mgcv ML fit this package reads curves off", {
  skip_if_not_installed("mgcv")
  d <- sp_gratia_data()
  fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 10)),
                     family = stats::gaussian(), data = d)
  gm <- mgcv::gam(y ~ s(x, k = 10), data = d, method = "ML")
  # frmtmb's marginal log likelihood IS mgcv's ML score
  expect_equal(as.numeric(stats::logLik(fit)),
               as.numeric(-gm$gcv.ubre), tolerance = 1e-6)
  g <- data.frame(x = seq(0, 1, length.out = 40))
  cv <- frm_curve(fit, newdata = g, simultaneous = FALSE)
  expect_equal(cv$.estimate,
               as.numeric(mgcv::predict.gam(gm, newdata = g)),
               tolerance = 1e-4)
})

test_that("the max-deviation simulation is gratia's, on gratia's own inputs", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("gratia")
  skip_on_cran()
  d <- sp_gratia_data()
  gm <- mgcv::gam(y ~ s(x, k = 10), data = d, method = "ML")
  g <- data.frame(x = seq(0, 1, length.out = 60))
  sm <- gm$smooth[[1L]]
  cs <- sm$first.para:sm$last.para
  Cg <- mgcv::PredictMat(sm, g)
  Xp <- mgcv::predict.gam(gm, newdata = g, type = "lpmatrix")
  # gratia divides a SMOOTH-ONLY deviation by smooth_estimates()'s .se,
  # which is the FULL predictor's standard error. Handing this package
  # the same two objects is what makes the critical values comparable;
  # standardizing the smooth-only deviation by its own standard error
  # instead gives a critical value 8 percent larger, and a band that is
  # just as exact.
  Ssm <- Cg %*% gm$Vp[cs, cs] %*% t(Cg)
  div <- sqrt(rowSums((Xp %*% gm$Vp) * Xp))
  nsim <- 100000L
  mine <- sp_sim_crit(Ssm, div, nsim, 0.95, seed = 11)
  ci <- stats::confint(gm, parm = "s(x)", type = "simultaneous", data = g,
                       nsim = nsim, level = 0.95)
  theirs <- unique(round((ci$.upper_ci - ci$.estimate) / ci$.se, 8))
  expect_length(theirs, 1L)
  expect_true(is.finite(mine$mcse) && mine$mcse > 0)
  # inside three Monte Carlo standard errors, which is the only claim
  # two independent simulations of the same quantile can support
  expect_lt(abs(mine$crit - theirs), 3 * mine$mcse)
})

test_that("gratia's .se really is the full predictor's, which is why the divisor is an argument", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("gratia")
  d <- sp_gratia_data()
  gm <- mgcv::gam(y ~ s(x, k = 10), data = d, method = "ML")
  g <- data.frame(x = seq(0, 1, length.out = 30))
  se_g <- gratia::smooth_estimates(gm, select = "s(x)", data = g)$.se
  Xp <- mgcv::predict.gam(gm, newdata = g, type = "lpmatrix")
  full <- sqrt(rowSums((Xp %*% gm$Vp) * Xp))
  sm <- gm$smooth[[1L]]
  cs <- sm$first.para:sm$last.para
  Cg <- mgcv::PredictMat(sm, g)
  only <- sqrt(rowSums((Cg %*% gm$Vp[cs, cs]) * Cg))
  expect_equal(unname(se_g), unname(full), tolerance = 1e-8)
  expect_gt(max(abs(se_g / only - 1)), 0.05)
})

test_that("first derivatives agree with gratia's to the fit difference", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("gratia")
  d <- sp_gratia_data()
  fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 10)),
                     family = stats::gaussian(), data = d)
  gm <- mgcv::gam(y ~ s(x, k = 10), data = d, method = "ML")
  g <- data.frame(x = seq(0.05, 0.95, length.out = 25))
  d1 <- frm_curve_deriv(fit, var = "x", order = 1, newdata = g,
                        simultaneous = FALSE)
  gd <- gratia::derivatives(gm, select = "s(x)", data = g, order = 1,
                            type = "central", eps = 1e-7)
  # the two fitted curves themselves differ by about 6e-07, and
  # differentiating divides a difference of that size by 2 eps
  expect_lt(max(abs(d1$.estimate - gd$.derivative)), 1e-4)
  # the standard errors differ by a few percent, and for a stated
  # reason: mgcv's Vp is conditional on the smoothing parameter and the
  # joint precision this package inverts is not
  expect_lt(max(abs(d1$.se / gd$.se - 1)), 0.15)
  expect_gt(max(abs(d1$.se / gd$.se - 1)), 0.005)
})

test_that("the second derivative is the one gratia's fixed eps gets wrong", {
  skip_if_not_installed("mgcv")
  d <- sp_gratia_data()
  fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 10)),
                     family = stats::gaussian(), data = d)
  gm <- mgcv::gam(y ~ s(x, k = 10), data = d, method = "ML")
  g <- data.frame(x = seq(0.05, 0.95, length.out = 25))
  d2 <- frm_curve_deriv(fit, var = "x", order = 2, newdata = g,
                        simultaneous = FALSE)
  # mgcv's own lpmatrix at the same stencil is the reference: it uses a
  # different package's basis code and the same step size
  stencil <- function(e) {
    gl <- g; gl$x <- g$x - e
    gh <- g; gh$x <- g$x + e
    (mgcv::predict.gam(gm, newdata = gh, type = "lpmatrix") -
       2 * mgcv::predict.gam(gm, newdata = g, type = "lpmatrix") +
       mgcv::predict.gam(gm, newdata = gl, type = "lpmatrix")) / e^2
  }
  ref <- as.numeric(stencil(1e-4 * diff(range(g$x))) %*% stats::coef(gm))
  expect_lt(max(abs(d2$.estimate - ref)), 1e-3)
  # and the same calculation at gratia's fixed 1e-7 is wrong by orders
  # of magnitude, which is the measurement the default eps rests on
  bad <- as.numeric(stencil(1e-7) %*% stats::coef(gm))
  expect_gt(max(abs(bad - ref)), 0.5)
})
