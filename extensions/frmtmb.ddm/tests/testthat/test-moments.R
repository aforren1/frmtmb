## The closed-form boundary probability and conditional mean decision
## time, against quadrature of RWiener's density. These are what
## fitted(), predict(type = "response") and residuals(type = "response")
## return, so a wrong constant here is a wrong fitted value.

qparts <- function(v, a, w, up) {
  resp <- if (up == 1) "upper" else "lower"
  tau <- 1e-9
  # dwiener warns that q and resp differ in length on every vectorized
  # call; it does the right thing, and the warning would drown the run
  f <- function(t) {
    suppressWarnings(RWiener::dwiener(t + tau, a, tau, w, v, resp = resp))
  }
  P <- stats::integrate(f, 0, 100, subdivisions = 4000L,
                        rel.tol = 1e-11)$value
  M <- stats::integrate(function(t) t * f(t), 0, 100,
                        subdivisions = 4000L, rel.tol = 1e-11)$value
  c(P = P, mean = M / P)
}

test_that("the boundary probability is the gambler's-ruin answer", {
  skip_if_not_installed("RWiener")
  gr <- expand.grid(v = c(-4, -2, -0.5, 0.5, 2, 4), a = c(0.7, 1.4, 2.5),
                    w = c(0.3, 0.5, 0.7))
  e <- mapply(function(v, a, w) {
    r <- qparts(v, a, w, 1)[["P"]]
    abs(ddm_p_upper(v, a, w) - r) / r
  }, gr$v, gr$a, gr$w)
  expect_lt(max(e), 1e-11)
  # at zero drift the start point IS the probability, and the expm1
  # spelling has to reach that exactly rather than through 0/0
  expect_equal(ddm_p_upper(0, 1.4, 0.3), 0.3)
  expect_equal(ddm_p_upper(0, 2.0, 0.5), 0.5)
  # the two boundaries exhaust the probability
  expect_equal(ddm_p_upper(1.2, 1.4, 0.4) +
                 ddm_p_upper(-1.2, 1.4, 0.6), 1)
})

test_that("the conditional mean decision time matches quadrature", {
  skip_if_not_installed("RWiener")
  gr <- expand.grid(v = c(-4, -2, -0.5, 0, 0.5, 2, 4), a = c(0.7, 1.4, 2.5),
                    w = c(0.3, 0.5, 0.7), up = c(0, 1))
  e <- mapply(function(v, a, w, up) {
    r <- qparts(v, a, w, up)[["mean"]]
    abs(ddm_cond_mean_dt(v, a, w, up) - r) / r
  }, gr$v, gr$a, gr$w, gr$up)
  expect_true(all(is.finite(e)))
  expect_lt(max(e), 1e-8)
})

test_that("the zero-drift crossover costs less than 1e-7", {
  skip_if_not_installed("RWiener")
  # The direct formulas cancel like eps / m^2 and the limits are wrong
  # like m^2, where m = 2 v a. The band boundary is where they meet.
  # Sweeping through it is the only way to know the seam is not a step.
  vs <- c(0, 1e-8, 1e-6, 1e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 1e-1)
  e <- vapply(vs, function(v) {
    r <- qparts(v, 1.4, 0.4, 1)[["mean"]]
    abs(ddm_cond_mean_dt(v, 1.4, 0.4, 1) - r) / r
  }, 0)
  expect_true(all(is.finite(e)))
  expect_lt(max(e), 1e-7)
  # and the seam is continuous to the accuracy the band was chosen for:
  # the step across it is the limit's own m^2 error, not a discontinuity
  m <- vapply(c(0.99, 1.01) * ddm_drift_tol / (2 * 1.4),
              function(v) ddm_cond_mean_dt(v, 1.4, 0.4, 1), 0)
  expect_equal(m[1], m[2], tolerance = 1e-7)
})

test_that("the zero-drift limits are the known polynomials", {
  a <- 1.7; w <- 0.35
  expect_equal(ddm_cond_mean_dt(0, a, w, 1), a^2 * (1 - w^2) / 3)
  expect_equal(ddm_cond_mean_dt(0, a, w, 0), a^2 * w * (2 - w) / 3)
})

test_that("the conditional mean obeys the boundary reflection", {
  # E[T | lower] at (v, w) is E[T | upper] at (-v, 1 - w): the same
  # reflection the density uses, so the moments and the density cannot
  # drift apart
  for (p in list(c(1.3, 1.4, 0.35), c(-0.7, 2.1, 0.6), c(0, 1.0, 0.5))) {
    expect_equal(ddm_cond_mean_dt(p[1], p[2], p[3], 0),
                 ddm_cond_mean_dt(-p[1], p[2], 1 - p[3], 1))
  }
})

test_that("an unbiased start gives both boundaries the same mean time", {
  # The classic drift-diffusion result: at w = 0.5 the conditional mean
  # decision time does not depend on which boundary was reached, for any
  # drift. It is worth pinning because it is the reason a fit with
  # bias fixed at 0.5 has fitted values that ignore the boundary, which
  # looks like a bug in the family until you know it is a theorem.
  for (v in c(-2, -0.5, 0, 0.5, 2, 5)) {
    expect_equal(ddm_cond_mean_dt(v, 1.4, 0.5, 1),
                 ddm_cond_mean_dt(v, 1.4, 0.5, 0))
  }
  # away from 0.5 the two differ, and by the reflection
  expect_gt(abs(ddm_cond_mean_dt(0.9, 1.4, 0.3, 1) -
                  ddm_cond_mean_dt(0.9, 1.4, 0.3, 0)), 0.1)
})

test_that("post$mean_fn is the mean of the sim slot", {
  skip_if_not_installed("RWiener")
  # An independent check of the closed form: the simulator draws by
  # inverse transform through RWiener's quantile function and shares no
  # code with the moment formulas.
  set.seed(21)
  dp <- list(mu = 1.1, bs = 1.4, ndt = 0.25, bias = 0.45)
  for (up in c(0, 1)) {
    at <- list(vint1 = rep(up, 4000))
    draws <- ddm_sim_rt(dp, at, 4000)
    m <- ddm_mean_rt(dp, list(vint1 = up))
    expect_true(all(is.finite(draws)))
    expect_true(all(draws > dp$ndt))
    # 4000 draws, so the standard error of the mean is the yardstick
    expect_lt(abs(mean(draws) - m), 4 * stats::sd(draws) / sqrt(4000))
  }
})

test_that("mean_fn and the simulator both refuse a missing boundary", {
  dp <- list(mu = 1, bs = 1.4, ndt = 0.2, bias = 0.5)
  expect_error(ddm_mean_rt(dp, list()), "conditional on the boundary")
  expect_error(ddm_sim_rt(dp, list(), 5), "needs the boundary")
})
