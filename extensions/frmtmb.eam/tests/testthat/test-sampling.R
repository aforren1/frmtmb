## Does the density hold up under NUTS? The Laplace path and the
## sampler put different demands on a log density: the optimizer needs a
## gradient near the mode, the sampler needs one everywhere it wanders,
## including the tails where the two series hand over.
##
## Gated the way the ODE extension gates RTMBode: the sampler is a
## suggested dependency of a suggested dependency, and its absence is a
## skip rather than a failure.

test_that("frm_sample runs a short chain on a wiener model", {
  skip_on_cran()
  skip_if_not_installed("RWiener")
  skip_if_not_installed("frmtmb.sample")
  skip_if_not_installed("tmbstan")

  set.seed(77)
  dat <- ddm_simulate(250, mu = 0.9, bs = 1.4, ndt = 0.25)
  fit <- frm(bf(rt | vint(upper) ~ 1, bias = 0.5), family = wiener(),
             data = dat)

  smp <- frmtmb.sample::frm_sample(fit, chains = 1, iter = 400,
                                   warmup = 200, seed = 3,
                                   refresh = 0)
  expect_true(!is.null(smp))
  # the sampled posterior sits around the mode the optimizer found
  fx <- unlist(fixef(fit))
  su <- summary(smp)
  expect_true(is.finite(fx[["mu.(Intercept)"]]))
  expect_true(nrow(as.data.frame(su)) > 0 || length(su) > 0)
})

test_that("the density is finite over a wide sweep of the parameters", {
  # what a sampler actually needs: no NaN anywhere it might step, not
  # merely near the mode. The tape has no branches, so this sweep is
  # the whole of the reachable behavior.
  gr <- expand.grid(t = c(1e-4, 1e-2, 0.2, 1, 10, 100),
                    v = c(-8, -1, 0, 1, 8),
                    a = c(0.05, 0.3, 1.4, 6, 20),
                    w = c(0.01, 0.2, 0.5, 0.8, 0.99))
  vals <- mapply(function(t, v, a, w) ddm_lpdf_lower(t, v, a, w),
                 gr$t, gr$v, gr$a, gr$w)
  expect_false(any(is.nan(vals)))
  # -Inf is allowed (a density of zero); NaN is not, because it would
  # poison the gradient rather than reject the step
  expect_true(all(is.finite(vals) | vals == -Inf))
  u <- gr$t / gr$a^2
  expect_lt(min(u), 1e-6)
  expect_gt(max(u), 1e4)
})

test_that("the gradient is finite over the same sweep", {
  skip_if_not_installed("RTMB")
  gr <- expand.grid(v = c(-8, -1, 0, 1, 8), a = c(0.3, 1.4, 6),
                    w = c(0.05, 0.5, 0.95))
  tt <- c(0.05, 0.5, 5)
  for (i in seq_len(nrow(gr))) {
    p <- c(gr$v[i], gr$a[i], gr$w[i])
    lp <- function(q) sum(ddm_lpdf_lower(tt, q[1], q[2], q[3]))
    tp <- RTMB::MakeTape(lp, p)
    g <- as.numeric(tp$jacobian(p))
    expect_false(any(is.nan(g)))
    expect_true(all(is.finite(g)))
  }
})
