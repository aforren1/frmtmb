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
  # Installed is not the same as working. A tmbstan built against
  # StanHeaders >= 2.39 samples a standard normal instead of the
  # model, and frm_sample() refuses it, so without this the block
  # ERRORS rather than skipping. No CI job installs both frmtmb.sample
  # and tmbstan here, but a developer machine does.
  skip_if(frmtmb.sample:::tmbstan_build_broken(),
          paste("this tmbstan samples a standard normal instead of",
                "the model; see dev/prior-dropping-investigation.md"))

  set.seed(77)
  dat <- ddm_simulate(250, mu = 0.9, bs = 1.4, ndt = 0.25)
  fit <- frm(bf(rt | vint(upper) ~ 1, bias = 0.5), family = wiener(),
             data = dat)

  # One chain of 200 draws: whether its ESS or R-hat crosses rstan's
  # warning line is rounding (both ESS warnings with OpenBLAS 0.3.26,
  # none with the reference BLAS; dev/ciharden-findings.md), and this
  # block asserts the chain runs and sits at the mode, not its mixing
  smp <- allow_warnings(
    frmtmb.sample::frm_sample(fit, chains = 1, iter = 400, warmup = 200,
                              seed = 3, refresh = 0),
    c("Effective Samples Size", "The largest R-hat is"))
  expect_true(!is.null(smp))
  # the sampled posterior sits around the mode the optimizer found
  fx <- unlist(fixef_by_dpar(fit))
  su <- summary(smp)
  expect_true(is.finite(fx[["mu.(Intercept)"]]))
  expect_true(nrow(as.data.frame(su)) > 0 || length(su) > 0)
  # With its ESS and R-hat warnings allowed, this is what says the chain
  # works: no divergent transition, and each parameter's posterior mean
  # within 3 posterior sds of the mode on the sampler's scale (measured
  # at most 0.152 sd with the reference BLAS and 0.096 with OpenBLAS
  # 0.3.32, dev/ciharden-rev-eamchain.R)
  sf <- smp$stanfit
  sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
  expect_identical(sum(vapply(sp, function(m) sum(m[, "divergent__"]), 0)),
                   0)
  est <- fit$opt$par
  s <- rstan::summary(sf)$summary
  z <- (s[seq_along(est), "mean"] - est) / s[seq_along(est), "sd"]
  expect_true(all(is.finite(z)))
  expect_lt(max(abs(z)), 3)
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
