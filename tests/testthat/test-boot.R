# Heavy reference-validation file: full runs happen locally and in CI
# (NOT_CRAN=true). CRAN only needs to see the package work, not the
# validation program re-derived, and this file is a main contributor
# to CRAN-condition check time.
skip_on_cran()

# v0.13: frm_bootstrap and hypothesis(method = "profile"/"boot"),
# including the natural-scale random-effect names.

sim_lmm_boot <- function(seed = 61, n_g = 15, n_per = 12) {
  set.seed(seed)
  dd <- data.frame(x = rnorm(n_g * n_per),
                   g = factor(rep(seq_len(n_g), each = n_per)))
  dd$y <- rnorm(nrow(dd), 1 + 0.5 * dd$x + rnorm(n_g, 0, 0.8)[dd$g], 1)
  dd
}

test_that("frm_bootstrap recovers the sampling distribution", {
  dd <- sim_lmm_boot()
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

  bs <- frm_bootstrap(fit, nsim = 40, seed = 1)
  expect_s3_class(bs, "frmtmb_boot")
  # default FUN covers every dpar's fixed effects, sigma included
  expect_equal(dim(bs$t), c(40L, 3L))
  expect_true(all(bs$converged))
  expect_named(bs$t0, c("mu.(Intercept)", "mu.x", "sigma.(Intercept)"))
  # bootstrap mean near the estimate, bootstrap SE near the Wald SE
  expect_lt(abs(mean(bs$t[, 2]) - bs$t0[[2]]), 0.1)
  se_wald <- sqrt(vcov(fit)["x", "x"])
  expect_lt(abs(stats::sd(bs$t[, 2]) - se_wald), 0.5 * se_wald)

  ci <- confint(bs)
  expect_equal(dim(ci), c(3L, 3L))
  expect_true(all(ci[, "lwr"] < ci[, "est"] & ci[, "est"] < ci[, "upr"]))
  expect_output(print(bs), "Parametric bootstrap")
})

test_that("hypothesis exposes natural-scale RE names", {
  dd <- sim_lmm_boot()
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

  h <- hypothesis(fit, "sd_g__Intercept")
  expect_equal(h$estimate, sqrt(VarCorr(fit)[[1]][1, 1]),
               tolerance = 1e-8)
  # Wald se agrees with confint_varcorr's delta method (log-scale
  # transform makes the intervals differ; the raw SEs should be close)
  cv <- confint_varcorr(fit)
  se_log <- (log(cv$upr[1]) - log(cv$lwr[1])) / (2 * stats::qnorm(0.975))
  expect_equal(h$se / h$estimate, se_log, tolerance = 0.05)

  icc <- "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)"
  hw <- hypothesis(fit, icc)
  vc <- VarCorr(fit)[[1]][1, 1]
  expect_equal(hw$estimate, vc / (vc + sigma(fit)^2), tolerance = 1e-8)

  hb <- hypothesis(fit, icc, method = "boot", nsim = 40, seed = 2)
  expect_true(hb$lwr > 0 && hb$upr < 1)
  expect_true(hb$lwr < hw$estimate && hw$estimate < hb$upr)

  # correlation name on a correlated-slopes fit
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  fs <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
            data = sleepstudy)
  hc <- hypothesis(fs, "cor_Subject__Intercept__Days")
  C <- stats::cov2cor(VarCorr(fs)[[1]])
  expect_equal(hc$estimate, C[1, 2], tolerance = 1e-8)

  # unknown names error with the available list
  expect_error(hypothesis(fit, "sd_h__Intercept"), "Available names")
})

test_that("hypothesis method = 'profile' handles linear contrasts", {
  set.seed(9)
  dd <- data.frame(x1 = rnorm(150), x2 = rnorm(150))
  dd$y <- rnorm(150, 1 + 0.6 * dd$x1 + 0.4 * dd$x2, 1)
  fit <- frm(bf(y ~ x1 + x2) + gaussian(), data = dd)

  # single coefficient: profile lincomb matches confint's uniroot
  # (interpolated from the tmbprofile curve, hence the tolerance)
  hp <- hypothesis(fit, "x1", method = "profile")
  cu <- confint(fit, parm = "x1", method = "uniroot")
  expect_equal(unname(c(hp$lwr, hp$upr)), unname(cu[1, 1:2]),
               tolerance = 1e-2)
  # the profile curve is stored and centered on the estimate
  pr <- attr(hp, "profiles")[[1]]
  expect_s3_class(pr, "data.frame")
  expect_true(min(pr[[1]]) < hp$estimate && hp$estimate < max(pr[[1]]))

  # contrast: profile interval close to Wald here (near-quadratic
  # likelihood), and the constant offset is handled
  hw <- hypothesis(fit, "x1 - x2 = 0.1")
  hp2 <- hypothesis(fit, "x1 - x2 = 0.1", method = "profile")
  expect_equal(hp2$estimate, hw$estimate, tolerance = 1e-10)
  expect_equal(unname(c(hp2$lwr, hp2$upr)), unname(c(hw$lwr, hw$upr)),
               tolerance = 0.02)

  # nonlinear expressions and REML fits are rejected
  expect_error(hypothesis(fit, "exp(x1)", method = "profile"),
               "not linear")
  fr <- frm(bf(y ~ x1 + x2) + gaussian(), data = dd, REML = TRUE)
  expect_error(hypothesis(fr, "x1", method = "profile"), "ML fit")
  # wald on an REML fit still works, including RE-free sd-less names
  expect_equal(hypothesis(fr, "x1")$estimate, fixef(fr)$mu[["x1"]],
               tolerance = 1e-10)
})

test_that("hypothesis method = 'boot' handles nonlinear expressions", {
  set.seed(13)
  dd <- data.frame(x = rnorm(120))
  dd$y <- rnorm(120, 1 + 0.5 * dd$x, 1)
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)

  hb <- hypothesis(fit, c("exp(x)", "x"), method = "boot", nsim = 60,
                   seed = 3)
  expect_s3_class(hb, "frmtmb_hypothesis")
  expect_equal(nrow(hb), 2L)
  # the draws ride along for interrogation and are coupled across rows
  d <- attr(hb, "draws")
  expect_equal(dim(d), c(60L, 2L))
  expect_equal(colnames(d), c("exp(x)", "x"))
  expect_equal(d[, 1], exp(d[, 2]), tolerance = 1e-10)
  expect_equal(hb$estimate[1], exp(fixef(fit)$mu[["x"]]),
               tolerance = 1e-10)
  # one shared bootstrap run: the draws are coupled, exp(x) row must be
  # consistent with the x row
  expect_true(hb$lwr[1] > 0)
  expect_true(hb$p[2] < 0.1)
  hw <- hypothesis(fit, "exp(x)")
  expect_lt(abs(hb$se[1] - hw$se[1]) / hw$se[1], 0.6)
})

test_that("backend controls pass through hypothesis's dots", {
  set.seed(23)
  dd <- sim_lmm_boot(seed = 81, n_g = 10, n_per = 8)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

  # profile: a wider ytol extends the stored curve
  h1 <- hypothesis(fit, "x", method = "profile")
  h2 <- hypothesis(fit, "x", method = "profile", ytol = 8)
  r1 <- diff(range(attr(h1, "profiles")[[1]][[1]]))
  r2 <- diff(range(attr(h2, "profiles")[[1]][[1]]))
  expect_gt(r2, r1)

  # boot: re.form reaches frm_bootstrap (conditional bootstrap runs)
  hb <- hypothesis(fit, "x", method = "boot", nsim = 10, seed = 5,
                   re.form = NULL)
  expect_equal(dim(attr(hb, "draws")), c(10L, 1L))

  # wald: stray arguments warn instead of vanishing
  expect_warning(hypothesis(fit, "x", ytol = 8), "unused by method")
})

test_that("hypothesis objects print and plot for every method", {
  set.seed(17)
  dd <- data.frame(x = rnorm(100))
  dd$y <- rnorm(100, 1 + 0.5 * dd$x, 1)
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)

  hw <- hypothesis(fit, c("x", "exp(x)"))
  hp <- hypothesis(fit, "x", method = "profile")
  hb <- hypothesis(fit, "x", method = "boot", nsim = 20, seed = 4)
  expect_output(print(hw), "method = wald")
  expect_output(print(hb), "bootstrap draws: 20")

  tmp <- file.path(tempdir(), "frmtmb-hyp-plots.pdf")
  grDevices::pdf(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  })
  expect_no_error(plot(hw))
  expect_no_error(plot(hp))
  expect_no_error(plot(hb))
})

test_that("hypothesis wald with RE names works under REML", {
  dd <- sim_lmm_boot(seed = 71)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd, REML = TRUE)
  h <- hypothesis(fit, "sd_g__Intercept")
  expect_equal(h$estimate, sqrt(VarCorr(fit)[[1]][1, 1]),
               tolerance = 1e-8)
  expect_true(is.finite(h$se) && h$se > 0)
})
