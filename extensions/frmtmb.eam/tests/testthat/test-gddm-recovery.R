## Does the estimator find the parameters it was given, and does
## renormalizing the defective density matter? Both questions are answered
## by fitting, so both are slow, and both are skipped on CRAN.

gd_truth <- list(mu = 6, leak = 1, bs = 3, tau = 1.2, ndt = 0.25)
gd_coh <- c(0, 0.128, 0.512)
gd_drift <- function() list(gddm_drift_constant(), gddm_drift_leak())

## A deliberately coarse grid, so the suite finishes. It is coarser than
## the shipped default and the drift rate is the parameter that notices.
gd_rctl <- function(renormalize = TRUE) {
  gddm_control(t_max = 2, dt = 0.02, ny = 101L, renormalize = renormalize)
}

## Draws come from the model's own density rather than from a forward
## simulation, whose first passages are late by however much the step
## misses excursions between monitoring times. Any gap left is the
## estimator's.
gd_rsim <- function(n_per, seed) {
  set.seed(seed)
  d <- do.call(rbind, lapply(gd_coh, function(cc)
    gddm_simulate(n_per, mu = gd_truth$mu * (cc / 0.512),
                  leak = gd_truth$leak, bs = gd_truth$bs, tau = gd_truth$tau,
                  ndt = gd_truth$ndt, drift = gd_drift(),
                  bound = gddm_bound_exponential(), control = gd_rctl())))
  d$coh <- rep(gd_coh, each = n_per)
  d$cond <- gddm_conditions(d, coh)
  d
}

gd_rfit <- function(d, renormalize = TRUE) {
  f <- frm(bf(rt | vint(upper, cond) ~ 0 + factor(coh), bias = 0.5),
           family = gddm(drift = gd_drift(), bound = gddm_bound_exponential(),
                         control = gd_rctl(renormalize)),
           data = d)
  e <- unlist(fixef(f))
  c(leak = unname(e[["leak.(Intercept)"]]),
    bs = exp(unname(e[["bs.(Intercept)"]])),
    tau = exp(unname(e[["tau.(Intercept)"]])),
    ndt = min(d$rt) / (1 + exp(-unname(e[["ndt.(Intercept)"]]))),
    mu_hi = unname(e[["mu.factor(coh)0.512"]]))
}

test_that("the estimator recovers the parameters it was given", {
  skip_on_cran()
  skip_on_ci()
  R <- 12L
  est <- t(vapply(seq_len(R), function(i) {
    r <- try(gd_rfit(gd_rsim(400L, 500L + i)), silent = TRUE)
    if (inherits(r, "try-error")) rep(NA_real_, 5L) else r
  }, numeric(5L)))
  ok <- stats::complete.cases(est)
  expect_gte(sum(ok), R - 1L)

  tru <- c(leak = gd_truth$leak, bs = gd_truth$bs, tau = gd_truth$tau,
           ndt = gd_truth$ndt, mu_hi = gd_truth$mu)
  m <- colMeans(est[ok, , drop = FALSE])
  mcse <- apply(est[ok, , drop = FALSE], 2, stats::sd) / sqrt(sum(ok))
  z <- (m - tru) / mcse

  ## Boundary height, its time constant, the leak and the non-decision
  ## time are the parameters a generalized model is fitted FOR, so they
  ## are the ones the test is strict about. The seeds are fixed, so this
  ## is deterministic; the bar is a little above the largest deviation
  ## measured, because the replicate-to-replicate spread on the leak is
  ## heavy enough that a standard error from twelve replicates understates
  ## it, and a bar set exactly at the measurement would be brittle.
  expect_lt(max(abs(z[c("leak", "bs", "tau", "ndt")])), 2.5)

  ## The drift rate is checked on a relative scale rather than in
  ## Monte Carlo standard errors. Across independent blocks of replicates
  ## its bias moved between about minus one and plus three percent, in
  ## both directions and at two grids, so there is no bias to pin, only a
  ## band to stay inside.
  expect_lt(abs(m[["mu_hi"]] / tru[["mu_hi"]] - 1), 0.06)
})

test_that("renormalizing the defective density is not optional", {
  skip_on_cran()
  ## The discretized solve loses probability mass, and how much it loses
  ## depends on the parameters: a configuration that absorbs faster loses
  ## less. A likelihood that does not divide it out therefore pays a
  ## hidden bonus for fast absorption, and the bonus lands exactly on the
  ## two parameters a generalized model exists to estimate.
  d <- gd_rsim(400L, 101L)
  a <- gd_rfit(d, renormalize = TRUE)
  b <- gd_rfit(d, renormalize = FALSE)

  ## the leak more than doubles
  expect_gt(b[["leak"]] / a[["leak"]], 2)
  ## and the boundary separation shrinks by at least a tenth
  expect_lt(b[["bs"]] / a[["bs"]], 0.9)
  ## in the direction that says "absorb sooner": more leak, closer walls
  expect_gt(b[["leak"]], a[["leak"]])
  expect_lt(b[["bs"]], a[["bs"]])

  ## and renormalizing is the one that is close to the truth on the leak
  expect_lt(abs(a[["leak"]] - gd_truth$leak),
            abs(b[["leak"]] - gd_truth$leak))
})

test_that("a fitted model supports the surface a user reaches for next", {
  skip_on_cran()
  d <- gd_rsim(200L, 77L)
  fit <- frm(bf(rt | vint(upper, cond) ~ 0 + factor(coh), bias = 0.5),
             family = gddm(drift = gd_drift(),
                           bound = gddm_bound_exponential(),
                           control = gd_rctl()),
             data = d)
  expect_s3_class(fit, "frmtmb_fit")
  expect_no_error(summary(fit))
  expect_true(is.finite(AIC(fit)))

  ## fitted() is the conditional mean response time for the row's own wall
  ft <- fitted(fit)
  expect_length(ft, nrow(d))
  expect_true(all(is.finite(ft)))
  expect_true(all(ft > 0))
  expect_equal(mean(ft), mean(d$rt), tolerance = 0.15)
  expect_equal(predict(fit, type = "response"), ft, tolerance = 1e-8)
  expect_no_error(predict(fit, type = "link"))

  ## the two walls have different mean response times under a leak
  expect_gt(stats::sd(ft), 0)

  ## simulate() draws response times of the right shape
  s <- simulate(fit, nsim = 1L)
  expect_equal(nrow(as.data.frame(s)), nrow(d))
})

test_that("the generalized family reproduces the analytic one end to end", {
  skip_on_cran()
  ## Same data, two families. With a constant drift and fixed bounds they
  ## are the same model, so their estimates should agree to about the
  ## solver's accuracy rather than to machine precision.
  skip_if_not_installed("RWiener")
  set.seed(31)
  d <- ddm_simulate(1200, mu = 1.2, bs = 2, ndt = 0.3, bias = 0.5)
  d$cond <- 1L
  aw <- frm(bf(rt | vint(upper) ~ 1, bias = 0.5), family = wiener(),
            data = d)
  ag <- frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5),
            family = gddm(control = gddm_control(t_max = max(d$rt) + 0.5,
                                                 dt = 0.01, ny = 201L)),
            data = d)
  ew <- unlist(fixef(aw))
  eg <- unlist(fixef(ag))
  expect_equal(unname(eg[["mu.(Intercept)"]]),
               unname(ew[["mu.(Intercept)"]]), tolerance = 0.05)
  expect_equal(exp(unname(eg[["bs.(Intercept)"]])),
               exp(unname(ew[["bs.(Intercept)"]])), tolerance = 0.05)
  ub <- min(d$rt)
  nd <- function(e) ub / (1 + exp(-unname(e[["ndt.(Intercept)"]])))
  expect_equal(nd(eg), nd(ew), tolerance = 0.05)
})

test_that("gddm_floored counts the rows the grid could not represent", {
  skip_on_cran()
  ## Clean data on an adequate grid: nothing should be at the floor, and
  ## the count is the evidence for that rather than an absence of
  ## warnings, which is what the floor removed.
  d <- gd_rsim(200L, 55L)
  fit <- frm(bf(rt | vint(upper, cond) ~ 0 + factor(coh), bias = 0.5),
             family = gddm(drift = gd_drift(),
                           bound = gddm_bound_exponential(),
                           control = gd_rctl()),
             data = d)
  n <- gddm_floored(fit)
  expect_type(n, "integer")
  expect_length(n, 1L)
  expect_identical(as.integer(n), 0L)
  expect_identical(attr(n, "n_obs"), nrow(d))
  expect_identical(attr(n, "rows"), integer(0))

  ## and it refuses a fit from another family rather than counting
  ## something that has no floor
  skip_if_not_installed("RWiener")
  set.seed(41)
  dw <- ddm_simulate(200, mu = 1, bs = 1.5, ndt = 0.3, bias = 0.5)
  fw <- frm(bf(rt | vint(upper) ~ 1, bias = 0.5), family = wiener(),
            data = dw)
  expect_error(gddm_floored(fw), "this is a wiener model")
})
