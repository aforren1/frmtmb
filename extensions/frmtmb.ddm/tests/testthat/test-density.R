## The density is the whole package. RWiener's C implementation is the
## external reference: it chooses its truncation adaptively, which is
## exactly what a tape cannot do, so agreeing with it over a wide range
## of normalized times is the claim being tested.

## RWiener refuses tau = 0, so the reference is taken a hair above it.
ref_lower <- function(t, v, a, w, tau = 1e-9) {
  log(RWiener::dwiener(t + tau, a, tau, w, v, resp = "lower"))
}
ref_upper <- function(t, v, a, w, tau = 1e-9) {
  log(RWiener::dwiener(t + tau, a, tau, w, v, resp = "upper"))
}

test_that("the blended density matches RWiener over both series regimes", {
  skip_if_not_installed("RWiener")
  gr <- expand.grid(t = c(0.01, 0.05, 0.15, 0.4, 0.8, 1.5, 2.5, 4, 8),
                    a = c(0.4, 0.8, 1.4, 2.2, 3.5),
                    w = c(0.2, 0.35, 0.5, 0.7, 0.85),
                    v = c(-3, -1.5, 0, 1.5, 3))
  # the grid spans normalized times from below 1e-3 to 50, which is both
  # series regimes and the band between them several times over
  u <- gr$t / gr$a^2
  expect_lt(min(u), 1e-3)
  expect_gt(max(u), 40)

  lo <- mapply(function(t, a, w, v) {
    r <- ref_lower(t, v, a, w)
    abs(ddm_lpdf_lower(t, v, a, w) - r) / abs(r)
  }, gr$t, gr$a, gr$w, gr$v)
  expect_true(all(is.finite(lo)))
  expect_lt(max(lo), 1e-11)

  # the upper boundary is the reflection, not a second series
  up <- mapply(function(t, a, w, v) {
    r <- ref_upper(t, v, a, w)
    abs(ddm_lpdf_both(t, v, a, w, 1) - r) / abs(r)
  }, gr$t, gr$a, gr$w, gr$v)
  expect_true(all(is.finite(up)))
  expect_lt(max(up), 1e-11)

  # up = 0 is the unreflected density, exactly
  expect_equal(ddm_lpdf_both(0.5, 1.2, 1.4, 0.4, 0),
               ddm_lpdf_lower(0.5, 1.2, 1.4, 0.4))
})

test_that("each series is accurate where the weight gives it any weight", {
  skip_if_not_installed("RWiener")
  # This is the design claim: the blend is only correct because both
  # series are accurate wherever the logistic weight is not saturated.
  # The band where the weight exceeds 1e-10 must sit strictly inside the
  # band where both series hold.
  us <- 10^seq(-2, 0.8, length.out = 40)
  lam <- 0.5 * (1 + tanh((log(us) - log(ddm_u0)) / ddm_us))
  live <- us[lam > 1e-10 & lam < 1 - 1e-10]
  expect_gt(length(live), 5)

  for (w in c(0.15, 0.5, 0.85)) {
    for (u in live) {
      r <- ref_lower(u, 0, 1, w)
      expect_lt(abs(ddm_log_gs(u, w) - r) / abs(r), 1e-12)
      expect_lt(abs(ddm_log_gl(u, w) - r) / abs(r), 1e-12)
    }
  }
})

test_that("the weight saturates, so a garbage series cannot leak in", {
  # Outside the band one series is nonsense. The weight has to be
  # exactly 0 or exactly 1 there, not merely small, and the density has
  # to stay finite either way.
  lam <- function(u) 0.5 * (1 + tanh((log(u) - log(ddm_u0)) / ddm_us))
  expect_equal(lam(1e-4), 0)
  expect_equal(lam(1e4), 1)
  # and the individual series stay finite even where they are wrong,
  # because 0 * NaN would still be NaN
  for (u in c(1e-6, 1e-3, 1, 1e3, 1e6)) {
    expect_true(is.finite(ddm_log_gs(u, 0.4)))
    expect_true(is.finite(ddm_log_gl(u, 0.4)))
  }
})

test_that("a single series would not do", {
  skip_if_not_installed("RWiener")
  # The motivation for paying for two series. At a normalized time of
  # 16 the small-time series is not slightly off, it is wrong, and more
  # terms do not rescue it: past the crossover its error is
  # cancellation, not truncation.
  r <- ref_lower(16, 0, 1, 0.45)
  expect_gt(abs(ddm_log_gs(16, 0.45) - r) / abs(r), 1e-2)
  expect_gt(abs(ddm_log_gs(16, 0.45, K = 40L) - r) / abs(r), 1e-2)
  # the blend is right there
  expect_lt(abs(ddm_lpdf_lower(16, 0, 1, 0.45) - r) / abs(r), 1e-12)
})

test_that("the density is zero below the non-decision time", {
  # Outside the support the value is not finite, which is what the
  # optimizer reads as a rejected step. This is the same contract
  # core's shifted_lognormal states for its own ndt.
  expect_false(is.finite(suppressWarnings(
    ddm_lpdf_lower(-0.1, 1, 1.4, 0.5))))
  expect_false(is.finite(suppressWarnings(
    ddm_lpdf_lower(0, 1, 1.4, 0.5))))

  # Strictly inside the support the log density stays a NUMBER, however
  # far down: as the decision time goes to zero the -w^2 / 2u term
  # outruns the -1.5 log(u) term, so the value falls off a cliff without
  # ever becoming undefined. That is the difference that matters,
  # because the bounded ndt link keeps the fit strictly inside, and a
  # finite gradient there is what steers the optimizer back.
  expect_false(is.nan(ddm_lpdf_lower(1e-8, 1, 1.4, 0.5)))
  expect_true(is.finite(ddm_lpdf_lower(1e-8, 1, 1.4, 0.5)))
  expect_lt(ddm_lpdf_lower(1e-8, 1, 1.4, 0.5), -1e6)
  expect_lt(ddm_lpdf_lower(1e-3, 1, 1.4, 0.5),
            ddm_lpdf_lower(0.5, 1, 1.4, 0.5) - 100)
  # monotone descent into the cliff, so there is no spurious mode there
  tt <- 10^seq(-6, -1, length.out = 20)
  expect_true(all(diff(ddm_lpdf_lower(tt, 1, 1.4, 0.5)) > 0))
})

test_that("the density integrates to the boundary probability", {
  # an internal consistency check that needs no external reference
  for (p in list(c(0.8, 1.4, 0.4), c(-1.5, 0.9, 0.6), c(0, 2.0, 0.5))) {
    v <- p[1]; a <- p[2]; w <- p[3]
    for (up in c(0, 1)) {
      m <- stats::integrate(
        function(t) exp(ddm_lpdf_both(t, v, a, w, up)),
        0, 200, subdivisions = 2000L, rel.tol = 1e-10)$value
      pu <- ddm_p_upper(v, a, w)
      expect_equal(m, if (up == 1) pu else 1 - pu, tolerance = 1e-7)
    }
  }
})
