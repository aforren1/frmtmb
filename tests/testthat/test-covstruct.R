test_that("ar1 matches glmmTMB", {
  skip_if_not_installed("glmmTMB")
  dd <- sim_ar1_data()
  fit <- frm(bf(y ~ 1 + ar1(tim + 0 | g)) + gaussian(), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ 1 + ar1(tim + 0 | g), data = dd,
                          REML = FALSE)
  expect_loglik_equal(fit, ref, tol = 1e-6)
  V <- VarCorr(fit)[[1]]
  Vg <- glmmTMB::VarCorr(ref)$cond$g
  expect_lt(abs(sqrt(V[1, 1]) - attr(Vg, "stddev")[1]), 1e-3)
  # estimated autocorrelation
  rho_frm <- V[1, 2] / V[1, 1]
  rho_ref <- attr(Vg, "correlation")[1, 2]
  expect_lt(abs(rho_frm - rho_ref), 1e-3)
})

test_that("cs matches glmmTMB", {
  skip_if_not_installed("glmmTMB")
  # genuinely compound-symmetric heteroscedastic data: shared group effect
  # plus time-specific noise scales
  set.seed(52)
  n_g <- 80; n_t <- 4
  g <- factor(rep(seq_len(n_g), each = n_t))
  tim <- factor(rep(seq_len(n_t), n_g))
  sd_t <- c(0.4, 0.6, 0.8, 1.0)
  u <- rnorm(n_g, 0, 0.8)
  dd <- data.frame(
    y = 1 + u[g] + rnorm(n_g * n_t, 0, sd_t[as.integer(tim)]),
    g = g, tim = tim
  )
  fit <- suppressWarnings(
    frm(bf(y ~ 1 + cs(tim + 0 | g)) + gaussian(), data = dd)
  )
  ref <- suppressWarnings(
    glmmTMB::glmmTMB(y ~ 1 + cs(tim + 0 | g), data = dd, REML = FALSE)
  )
  expect_loglik_equal(fit, ref, tol = 1e-5)
})

test_that("homdiag matches glmmTMB", {
  skip_if_not_installed("glmmTMB")
  dd <- sim_ar1_data(seed = 53, rho = 0)
  fit <- frm(bf(y ~ 1 + homdiag(tim + 0 | g)) + gaussian(), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ 1 + homdiag(tim + 0 | g),
                          data = dd, REML = FALSE)
  expect_loglik_equal(fit, ref, tol = 1e-6)
})

test_that("ar1 validation errors are clear", {
  dd <- sim_ar1_data(n_g = 5)
  expect_error(frm(bf(y ~ 1 + ar1(tim | g)) + gaussian(), data = dd),
               "ar1\\(times \\+ 0 \\| g\\)")
  dd$xnum <- rnorm(nrow(dd))
  expect_error(frm(bf(y ~ 1 + ar1(xnum + 0 | g)) + gaussian(), data = dd),
               "at least 2")
})
