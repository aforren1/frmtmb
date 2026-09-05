## The whole point of the change of variable is that the likelihood
## differentiates. These are the tests that say so: automatic derivatives
## through the PDE solve, the shift and the renormalization, against
## numerical ones, at parameter points including a collapsing boundary.
##
## numDeriv is given a wider step than it takes by default. Each density
## value is the end of a long chain of floating-point work, so at the
## default step the difference between two evaluations is mostly
## cancellation noise, and Richardson extrapolation then extrapolates the
## noise: at one point below it disagrees with the automatic derivative by
## 4e-4 relative while a plain central difference at 1e-3 agrees with it to
## 6e-7. The wider step is the honest comparison, not a loosened one.
gd_nd <- function(f, v) {
  numDeriv::grad(f, v, method.args = list(d = 0.01, r = 6))
}

test_that("AD gradients through the solver match numDeriv", {
  skip_if_not_installed("numDeriv")
  skip_if_not_installed("RTMB")
  comp <- gddm(drift = list(gddm_drift_constant(), gddm_drift_leak()),
               bound = gddm_bound_exponential())[["gddm"]]$comp
  ctl <- list(dt = 0.02, ny = 101L, t_max = 2, nt = 100L,
              renormalize = TRUE, max_ndt = 0.3, wmax = 17L,
              dpars = c("mu", "leak", "bs", "tau", "bias", "ndt"))
  d <- list(ncond = 1L, first = 1L, cov = matrix(0), gindex = 1L)
  dens <- function(v) {
    dp <- list(mu = v[1], leak = v[2], bs = exp(v[3]), tau = exp(v[4]),
               bias = 1 / (1 + exp(-v[5])), ndt = v[6])
    gd_densities(dp, comp, ctl, d)
  }
  points <- list(
    ## perfect integration, a bound that barely moves
    c(1.5, 0, log(2.5), log(20), 0, 0.2),
    ## leaky integration, a bound collapsing on the timescale of the window
    c(2.0, 1.2, log(3.0), log(1.0), 0, 0.25),
    ## unstable integration, a faster collapse, a biased start
    c(0.8, -0.7, log(2.0), log(0.8), 0.4, 0.15))
  for (v0 in points) {
    ## Read the density where a fit would read it: at the times carrying
    ## the most mass. Fixed indices would sometimes land in the far tail,
    ## where a defective density underflows and the log is not the thing
    ## under test.
    p0 <- as.numeric(dens(v0))
    idx <- order(p0, decreasing = TRUE)[c(1, 5, 20, 60, 120)]
    expect_true(all(p0[idx] > 0))
    ## a functional shaped like a log likelihood
    f <- function(v) sum(log(gd_densities(
      list(mu = v[1], leak = v[2], bs = exp(v[3]), tau = exp(v[4]),
           bias = 1 / (1 + exp(-v[5])), ndt = v[6]),
      comp, ctl, d)[idx]))
    tp <- RTMB::MakeTape(f, v0)
    ad <- as.numeric(tp$jacfun()(v0))
    ## 1e-4 relative is numDeriv's own floor on this function, not a
    ## slack tolerance: no single finite-difference step serves all six
    ## components at once, and whichever step is chosen, the component it
    ## suits least is the one that disagrees.
    expect_equal(ad, gd_nd(f, v0), tolerance = 1e-4)
    ## and the taped value is the plain value
    expect_equal(as.numeric(tp(v0)), f(v0), tolerance = 1e-10)
  }
})

test_that("the hand-written tridiagonal adjoint is the true derivative", {
  skip_if_not_installed("numDeriv")
  skip_if_not_installed("RTMB")
  ## gd_tri_solve is one atomic tape node with a derivative written by
  ## hand, so nothing else in the package would notice if it were wrong.
  n <- 7L
  mk <- function(p) c(p[1] * (1:(n - 1)) / n, 2 + p[2] * (1:n) / n,
                      p[3] * ((n - 1):1) / n, p[4] * cos(1:n))
  g <- function(p) {
    y <- gd_tri_solve(mk(p))
    sum(y * y * (1:n))
  }
  p0 <- c(0.3, 0.5, -0.4, 1.1)
  tp <- RTMB::MakeTape(g, p0)
  expect_equal(as.numeric(tp$jacfun()(p0)), gd_nd(g, p0), tolerance = 1e-7)
  ## and it survives a second derivative, which is what standard errors need
  expect_true(all(is.finite(as.numeric(tp$jacfun()$jacfun()(p0)))))
})

test_that("zero coherence has a finite gradient in the exponent", {
  skip_if_not_installed("numDeriv")
  skip_if_not_installed("RTMB")
  ## d/d(alpha) of C^alpha is C^alpha log(C), which is NaN at C = 0. A
  ## coherence design normally contains a zero condition, so this is not a
  ## corner case, it is the first row of the design.
  comp <- gddm(drift = gddm_drift_coherence(cmax = 0.512))[["gddm"]]$comp
  ctl <- list(dt = 0.02, ny = 101L, t_max = 2, nt = 100L,
              renormalize = TRUE, max_ndt = 0.3, wmax = 17L,
              dpars = c("mu", "alpha", "bs", "bias", "ndt"))
  ## the zero condition alongside two nonzero ones
  d <- list(ncond = 3L, first = 1:3, cov = matrix(c(0, 0.128, 0.512)),
            gindex = 1:3)
  f <- function(v) {
    dp <- list(mu = rep(v[1], 3), alpha = rep(exp(v[2]), 3),
               bs = rep(exp(v[3]), 3), bias = rep(0.5, 3),
               ndt = rep(v[4], 3))
    p <- gd_densities(dp, comp, ctl, d)
    sum(log(p[idx]))
  }
  v0 <- c(6, log(0.8), log(2.5), 0.2)
  ## read each condition where it carries mass, not at fixed offsets
  p0 <- as.numeric({
    dp <- list(mu = rep(v0[1], 3), alpha = rep(exp(v0[2]), 3),
               bs = rep(exp(v0[3]), 3), bias = rep(0.5, 3),
               ndt = rep(v0[4], 3))
    gd_densities(dp, comp, ctl, d)
  })
  idx <- unlist(lapply(0:5, function(b) {
    o <- b * 101L + seq_len(101L)
    o[order(p0[o], decreasing = TRUE)[c(1, 8)]]
  }))
  expect_true(all(p0[idx] > 0))
  tp <- RTMB::MakeTape(f, v0)
  ad <- as.numeric(tp$jacfun()(v0))
  expect_true(all(is.finite(ad)))
  expect_equal(ad, gd_nd(f, v0), tolerance = 1e-6)

  ## the zero condition alone carries no dependence on the exponent at
  ## all, which is why it can be dropped from the tape rather than
  ## regularized: the drift is identically zero whatever alpha is
  d0 <- list(ncond = 1L, first = 1L, cov = matrix(0), gindex = 1L)
  g <- function(v) {
    dp <- list(mu = v[1], alpha = exp(v[2]), bs = exp(v[3]), bias = 0.5,
               ndt = v[4])
    p <- gd_densities(dp, comp, ctl, d0)
    sum(log(p[i0]))
  }
  p00 <- as.numeric(gd_densities(
    list(mu = v0[1], alpha = exp(v0[2]), bs = exp(v0[3]), bias = 0.5,
         ndt = v0[4]), comp, ctl, d0))
  i0 <- order(p00, decreasing = TRUE)[c(1, 10, 40)]
  ga <- as.numeric(RTMB::MakeTape(g, v0)$jacfun()(v0))
  expect_true(is.finite(ga[2]))
  expect_equal(ga[2], 0)

  ## and a nonzero coherence does depend on it
  d1 <- list(ncond = 1L, first = 1L, cov = matrix(0.128), gindex = 1L)
  g1 <- function(v) {
    dp <- list(mu = v[1], alpha = exp(v[2]), bs = exp(v[3]), bias = 0.5,
               ndt = v[4])
    p <- gd_densities(dp, comp, ctl, d1)
    sum(log(p[i1]))
  }
  p11 <- as.numeric(gd_densities(
    list(mu = v0[1], alpha = exp(v0[2]), bs = exp(v0[3]), bias = 0.5,
         ndt = v0[4]), comp, ctl, d1))
  i1 <- order(p11, decreasing = TRUE)[c(1, 10, 40)]
  expect_gt(abs(as.numeric(RTMB::MakeTape(g1, v0)$jacfun()(v0))[2]), 1e-3)
})

test_that("a negative coherence flips the drift rather than failing", {
  ## stimulus coding: the sign is data, so the exponent applies to the
  ## magnitude and the sign is carried through unchanged
  comp <- gddm(drift = gddm_drift_coherence(cmax = 0.512))[["gddm"]]$comp
  p <- list(mu = 4, alpha = 0.8)
  x <- c(-1, 0, 1)
  pos <- comp$drift$fn(x, 0, p, 0.128)
  neg <- comp$drift$fn(x, 0, p, -0.128)
  expect_equal(as.numeric(pos), -as.numeric(neg))
  expect_true(all(is.finite(as.numeric(neg))))
})

test_that("a fitted model's standard errors are finite", {
  skip_on_cran()
  set.seed(21)
  d <- gddm_simulate(600, mu = 2.5, bs = 2.5, ndt = 0.25,
                     control = gddm_control(t_max = 2))
  d$cond <- 1L
  fit <- frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5),
             family = gddm(control = gddm_control(t_max = 2, dt = 0.02,
                                                  ny = 101L)),
             data = d)
  se <- summary(fit)
  expect_s3_class(fit, "frmtmb_fit")
  expect_true(all(is.finite(unlist(fixef(fit)))))
  expect_true(is.finite(as.numeric(logLik(fit))))
})
