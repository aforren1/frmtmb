test_that("sleepstudy LMM matches lmer and glmmTMB (ML)", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("glmmTMB")
  data(sleepstudy, package = "lme4")

  fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
                data = sleepstudy)

  ref_lmer <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy,
                         REML = FALSE)
  expect_loglik_equal(fit, ref_lmer, tol = 1e-6)
  # within 1e-3 of the coefficients' standard errors (near 6.8 and 1.5):
  # an absolute 1e-4 asked more than either optimizer's own test, and
  # the two agree to 1.1e-4 at log-likelihoods equal to 1.2e-8
  # (dev/optrap-findings.md)
  se_mu <- summary(fit)$coefficients$mu[, "Std. Error"]
  expect_lt(max(abs(unname(fixef_by_dpar(fit)$mu) -
                      unname(lme4::fixef(ref_lmer))) / se_mu), 1e-3)

  ref_tmb <- glmmTMB::glmmTMB(Reaction ~ Days + (Days | Subject),
                              sleepstudy, REML = FALSE)
  expect_loglik_equal(fit, ref_tmb, tol = 1e-6)

  # as a ratio: glmmTMB's own optimum is 7.7e-4 from lme4's in the
  # intercept, and this fit's 1.1e-4, so the two SEs differ by more than
  # an absolute 1e-4 on an SE near 6.8 (dev/optrap-findings.md)
  se_f <- summary(fit)$coefficients$mu[, "Std. Error"]
  se_g <- summary(ref_tmb)$coefficients$cond[, "Std. Error"]
  expect_vector_equal(se_f / se_g, rep(1, length(se_g)), tol = 1e-4)

  vc_f <- varcorr_matrices(fit)[[1]]
  vc_g <- glmmTMB::VarCorr(ref_tmb)$cond$Subject
  expect_vector_equal(sqrt(diag(vc_f)), attr(vc_g, "stddev"), tol = 1e-3)
})

test_that("sleepstudy REML matches lmer REML", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")

  fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
                data = sleepstudy, REML = TRUE)
  ref <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy,
                    REML = TRUE)
  expect_loglik_equal(fit, ref, tol = 1e-6)
  expect_vector_equal(fixef_by_dpar(fit)$mu, lme4::fixef(ref), tol = 1e-4)
  # REML fixed-effect SEs
  se_f <- summary(fit)$coefficients$mu[, "Std. Error"]
  se_l <- coef(summary(ref))[, "Std. Error"]
  expect_vector_equal(se_f, se_l, tol = 1e-4)
})

test_that("uncorrelated slopes (x || g) match glmmTMB", {
  skip_if_not_installed("glmmTMB")
  data(sleepstudy, package = "lme4")

  fit <- frm(bf(Reaction ~ Days + (Days || Subject)) + gaussian(),
                data = sleepstudy)
  ref <- glmmTMB::glmmTMB(Reaction ~ Days + (Days || Subject),
                          sleepstudy, REML = FALSE)
  expect_loglik_equal(fit, ref, tol = 1e-6)
})

test_that("linear model without random effects matches lm", {
  set.seed(61)   # unseeded, this failed ~1 in 3 runs
  dd <- data.frame(y = rnorm(60), x = rnorm(60))
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  ref <- stats::lm(y ~ x, dd)
  expect_loglik_equal(fit, ref, tol = 1e-6)
  # a 1e-5 coefficient displacement costs ~1e-9 logLik here, below any
  # optimizer stopping rule, so the likelihood gate above is the sharp
  # one and the coefficient gate allows the flat-region slack
  expect_vector_equal(fixef_by_dpar(fit)$mu, coef(ref), tol = 5e-5)
})
