## The generalized family has no closed-form density, so the solver is
## checked against the one case that does have one: with a constant drift
## and boundaries that do not move, the generalized model IS the Wiener
## model, and the flux through each wall must reproduce this package's own
## Navarro-Fuss density.

gd_ctl <- function(dt = 0.01, ny = 201L, t_max = 3, renormalize = FALSE,
                   max_ndt = 0.4) {
  nt <- as.integer(round(t_max / dt))
  list(dt = dt, ny = as.integer(ny), t_max = nt * dt, nt = nt,
       renormalize = renormalize, max_ndt = max_ndt,
       wmax = as.integer(ceiling(max_ndt / dt)) + 2L)
}

## the family's own lookup: a response time is data, so the two nodes it
## falls between and the weight between them are constants
gd_look <- function(v, tt, dt) {
  s <- tt / dt
  k <- floor(s)
  w <- s - k
  (1 - w) * v[k + 1] + w * v[k + 2]
}

gd_wiener_cases <- list(
  list(mu = 1, bs = 2, bias = 0.5),
  list(mu = -0.8, bs = 1.5, bias = 0.5),
  list(mu = 2, bs = 1.2, bias = 0.35),
  list(mu = 0, bs = 1.8, bias = 0.5),
  list(mu = 3, bs = 2.5, bias = 0.6))

## Worst log-density error against the analytic family, over both walls and
## every case, for decision times from `from_t` on. The window is an
## absolute time, not a multiple of the step, because comparing grids at a
## fixed multiple would compare them over different windows: a fixed grid
## cannot resolve the leading edge, where the exact density climbs through
## orders of magnitude inside one cell, so a finer grid measured closer to
## zero can look worse while being better everywhere they overlap.
gd_worst_vs_wiener <- function(dt, ny, from_t) {
  comp <- gddm()[["gddm"]]$comp
  ctl <- gd_ctl(dt = dt, ny = ny)
  worst <- 0
  for (cs in gd_wiener_cases) {
    p <- c(cs, list(ndt = 0))
    s <- gd_solve(p, 0, comp, ctl)
    tt <- seq(from_t, 2, by = 0.0025)
    for (up in c(1, 0)) {
      ex <- exp(ddm_lpdf_both(tt, p$mu, p$bs, p$bias, up))
      num <- gd_look(as.numeric(if (up == 1) s$up else s$lo), tt, dt)
      k <- num > 0 & ex > 0
      worst <- max(worst, max(abs(log(num[k]) - log(ex[k]))))
    }
  }
  worst
}

test_that("with a constant drift and fixed bounds the solver is the Wiener density", {
  # At the shipped grid, over decision times from 0.2 s on. Measured, not
  # guessed: this is the accuracy the help page claims.
  expect_lt(gd_worst_vs_wiener(0.01, 201L, 0.2), 0.01)
})

test_that("the density degrades gracefully on a coarser grid", {
  # compared over one window, so that the grids are compared and not the
  # windows they are read over
  coarse <- gd_worst_vs_wiener(0.02, 101L, 0.2)
  default <- gd_worst_vs_wiener(0.01, 201L, 0.2)
  fine <- gd_worst_vs_wiener(0.005, 401L, 0.2)
  expect_lt(coarse, 0.1)
  expect_lt(default, coarse)
  expect_lt(fine, default)
})

test_that("refining the spatial grid alone does not make the density worse", {
  # A delta start excites every frequency the grid carries and
  # Crank-Nicolson does not damp the highest of them, so a start spread
  # over one cell rings harder as ny grows. The shipped start is spread
  # over two, and this is the test that pins why.
  at_ny <- vapply(c(101L, 201L, 401L),
                  function(ny) gd_worst_vs_wiener(0.01, ny, 0.2), numeric(1))
  expect_lt(max(at_ny), 0.02)
})

test_that("the leading edge is where the density is not to be trusted", {
  # Stated rather than hidden: an implicit scheme spreads a little mass
  # everywhere at once, so at decision times of a few steps the density is
  # far larger than the truth. It is small in absolute terms, and a lapse
  # component floors it, but a fit should not be asked to read it.
  near <- gd_worst_vs_wiener(0.01, 201L, 0.03)
  far <- gd_worst_vs_wiener(0.01, 201L, 0.2)
  expect_gt(near, 10 * far)
})

test_that("the solver conserves probability against the analytic family", {
  comp <- gddm()[["gddm"]]$comp
  ctl <- gd_ctl(dt = 0.005, ny = 201L, t_max = 4)
  for (cs in gd_wiener_cases) {
    p <- c(cs, list(ndt = 0))
    s <- gd_solve(p, 0, comp, ctl)
    tt <- seq(0, ctl$t_max, by = ctl$dt)
    mass <- (sum(s$up) + sum(s$lo)) * ctl$dt
    # everything that has not absorbed by t_max is the only shortfall
    expect_gt(mass, 0.98)
    expect_lt(mass, 1 + 1e-6)
  }
})

test_that("the wall flux is never negative for a well-resolved model", {
  comp <- gddm(drift = list(gddm_drift_constant(), gddm_drift_leak()),
               bound = gddm_bound_exponential())[["gddm"]]$comp
  ctl <- gd_ctl(dt = 0.01, ny = 201L, t_max = 2)
  s <- gd_solve(list(mu = 1.5, leak = 0.5, bs = 3, tau = 1, bias = 0.5,
                     ndt = 0), 0, comp, ctl)
  expect_true(all(as.numeric(s$up) >= 0))
  expect_true(all(as.numeric(s$lo) >= 0))
})

test_that("the non-decision-time shift conserves mass and moves the density", {
  comp <- gddm()[["gddm"]]$comp
  ctl <- gd_ctl(dt = 0.01, ny = 101L, t_max = 3)
  s <- gd_solve(list(mu = 1.5, bs = 2, bias = 0.5, ndt = 0), 0, comp, ctl)
  v <- as.numeric(s$up)
  # The kernel is a partition of unity, so every admissible shift keeps
  # the mass except what leaves the end of the window. Zero is included
  # deliberately: the kernel is centered, so at a shift of under two steps
  # part of it reaches backward, and a version that dropped those taps
  # lost a sixth of the mass exactly where the non-decision time is small.
  for (ndt in c(0, 0.005, 0.05, 0.137, 0.3)) {
    sh <- as.numeric(gd_shift(v, ndt, ctl$dt, ctl$wmax))
    # nothing is created, and the only loss is what the shift pushes past
    # the end of the window
    expect_lte(sum(sh), sum(v) * (1 + 1e-9))
    expect_gt(sum(sh), sum(v) * 0.99)
    # the weights are non-negative, so the shift invents no new negatives
    expect_gte(min(sh), min(0, min(v)) - 1e-12)
    # the shifted density peaks one non-decision time later
    expect_equal(which.max(sh) - which.max(v), round(ndt / ctl$dt),
                 tolerance = 1)
  }
  # a shift is a smoothing, not an interpolation, so it does not sharpen
  expect_lt(max(as.numeric(gd_shift(v, 0.1, ctl$dt, ctl$wmax))), max(v))
})

test_that("a collapsing bound really does collapse, and bounds stay positive", {
  ctl <- gd_ctl(t_max = 2)
  ex <- gddm_bound_exponential()$fn
  b0 <- ex(0, list(bs = 3, tau = 1), ctl)
  b1 <- ex(1, list(bs = 3, tau = 1), ctl)
  expect_equal(b0$B, 1.5)
  expect_equal(b1$B, 1.5 * exp(-1))
  expect_equal(b0$dlogB, -1)

  li <- gddm_bound_linear()$fn
  l0 <- li(0, list(bs = 3, kappa = 0.9), ctl)
  l1 <- li(ctl$t_max, list(bs = 3, kappa = 0.9), ctl)
  expect_equal(l0$B, 1.5)
  # kappa is the fraction lost by the end of the window, so the bound is
  # strictly positive across it for any kappa the logit link can reach
  expect_equal(l1$B, 1.5 * 0.1)
  expect_gt(l1$B, 0)

  cn <- gddm_bound_constant()$fn
  expect_equal(cn(0, list(bs = 2), ctl)$B, cn(5, list(bs = 2), ctl)$B)
  expect_equal(cn(1, list(bs = 2), ctl)$dlogB, 0)
})

test_that("both starting distributions carry unit mass and sit where asked", {
  h <- 2 / 202
  yg <- seq(-1 + h, 1 - h, length.out = 201L)
  for (bias in c(0.3, 0.5, 0.72)) {
    r <- gddm_start_point()$fn(yg, h, list(bias = bias))
    expect_equal(sum(r) * h, 1, tolerance = 1e-12)
    expect_equal(sum(yg * r) * h, 2 * bias - 1, tolerance = 1e-6)
  }
  r <- gddm_start_uniform()$fn(yg, h, list(bias = 0.5, sz = 0.3))
  expect_equal(sum(r) * h, 1, tolerance = 1e-12)
  expect_equal(sum(yg * r) * h, 0, tolerance = 1e-6)
  # a uniform start is wider than a point start
  rp <- gddm_start_point()$fn(yg, h, list(bias = 0.5))
  expect_gt(sum(yg^2 * r) * h, sum(yg^2 * rp) * h)
})

test_that("the tridiagonal solve is a solve", {
  n <- 9L
  set.seed(4)
  lo <- runif(n - 1, -0.3, 0.3)
  up <- runif(n - 1, -0.3, 0.3)
  di <- runif(n, 2, 3)
  rhs <- rnorm(n)
  A <- diag(di, n)
  for (i in seq_len(n - 1)) {
    A[i + 1, i] <- lo[i]
    A[i, i + 1] <- up[i]
  }
  expect_equal(gd_tri_f(c(lo, di, up, rhs)), as.numeric(solve(A, rhs)),
               tolerance = 1e-12)
})

## The floor. A parameter value whose density the grid cannot represent
## at some observed response time has to give the optimizer a number, not
## a NaN: NaN is not a value a line search can use, and inside mixture()
## it takes every other component with it through the log-sum-exp.

test_that("the floor is a floor, not a smooth maximum", {
  # The smooth maximum (x + lo + |x - lo|)/2 is wrong at these
  # magnitudes: the +lo and -lo both vanish into the rounding of x and
  # it collapses to the positive part, returning exactly 0 for a
  # negative x, which is the value the floor exists to prevent.
  lo <- 1e-300
  smooth_max <- function(x, lo) 0.5 * (x + lo + abs(x - lo))
  expect_identical(smooth_max(-1e-17, lo), 0)
  expect_identical(gd_floor(-1e-17, lo), lo)

  # inert on anything a density legitimately takes
  for (x in c(1e-8, 1e-3, 0.5, 2.5, 1e3)) {
    expect_identical(gd_floor(x, lo), x)
  }
  # and never returns something log() cannot take
  expect_true(all(gd_floor(c(-1, -1e-320, 0, 1e-320), lo) > 0))
})

test_that("a row the grid cannot represent is flat, not undefined", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("numDeriv")
  # log(gd_floor(x)) below the floor: finite, very negative, and with a
  # derivative of exactly zero, which is what makes it flat rather than
  # undefined. -Inf would differentiate to NaN, the same failure again.
  f <- function(v) log(gd_floor(v[1L], gd_dens_floor))
  below <- -1e-20
  expect_true(is.finite(f(below)))
  expect_lt(f(below), -600)
  g <- as.numeric(RTMB::MakeTape(f, below)$jacfun()(below))
  expect_identical(g, 0)
  expect_false(anyNA(g))
  # above the floor it is the ordinary log, with the ordinary derivative
  expect_equal(f(0.25), log(0.25))
  expect_equal(as.numeric(RTMB::MakeTape(f, 0.25)$jacfun()(0.25)), 1 / 0.25)
})

test_that("a floored row is invisible to a mixture's log-sum-exp", {
  # The property that matters is not that exp() is literally zero -- the
  # floor is a positive number, so it is not -- but that the row cannot
  # move a component beside it. Adding it to any ordinary log density in
  # log space returns that density bit for bit.
  lp <- log(gd_floor(-1e-20, gd_dens_floor))
  expect_true(is.finite(lp))
  expect_lt(lp, -600)
  logspace_add <- function(a, b) {
    m <- pmax(a, b)
    m + log1p(exp(pmin(a, b) - m))
  }
  for (other in c(-1.5, 2.25, -12)) {
    expect_identical(logspace_add(other, lp), other)
  }
  # at a log density of exactly 0 there is no magnitude to absorb it, and
  # what is left is the floor itself, which is as small as a double gets
  # without being zero
  expect_lt(abs(logspace_add(0, lp)), 1e-299)
  # and it is a number, which is the whole point: NaN is not
  expect_false(is.nan(lp))
})

test_that("extreme parameters give a bad number, never no number", {
  # The optimizer overshoots into these. A collapsing bound driven to
  # zero makes the diffusion coefficient 0.5/B^2 infinite, and an Inf on
  # both the diagonal and the off-diagonals of the tridiagonal system
  # comes back from the sweep as NaN. NaN is not a value a line search
  # can use, so the boundary is floored and every one of these returns a
  # finite number instead.
  comp <- gddm(bound = gddm_bound_exponential())[["gddm"]]$comp
  ctl <- list(dt = 0.02, ny = 51L, t_max = 2, nt = 100L,
              renormalize = TRUE, max_ndt = 0.3, wmax = 17L,
              dpars = c("mu", "bs", "tau", "bias", "ndt"),
              tridiagonal = "recorded")
  d <- list(ncond = 1L, first = 1L, cov = matrix(0), gindex = 1L)
  extremes <- list(
    list(bs = 2.5, tau = 1e-3), list(bs = 2.5, tau = 1e-8),
    list(bs = 1e-6, tau = 1), list(bs = 1e6, tau = 1),
    list(bs = 1e-4, tau = 1e-4))
  for (e in extremes) {
    v <- as.numeric(gd_densities(
      list(mu = 1, bs = e$bs, tau = e$tau, bias = 0.5, ndt = 0.2),
      comp, ctl, d))
    expect_false(anyNA(v))
    expect_true(all(is.finite(v)))
    # and once logged through the floor it is still a number
    lp <- log(gd_floor(v, gd_dens_floor))
    expect_false(anyNA(lp))
  }

  # a linear bound approaches zero the same way, as its collapse
  # fraction approaches one
  compl <- gddm(bound = gddm_bound_linear())[["gddm"]]$comp
  ctll <- ctl
  ctll$dpars <- c("mu", "bs", "kappa", "bias", "ndt")
  v <- as.numeric(gd_densities(
    list(mu = 1, bs = 2, kappa = 1 - 1e-12, bias = 0.5, ndt = 0.2),
    compl, ctll, d))
  expect_false(anyNA(v))
})
