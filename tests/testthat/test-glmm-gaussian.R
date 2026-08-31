test_that("sleepstudy LMM matches lmer and glmmTMB (ML)", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("glmmTMB")
  data(sleepstudy, package = "lme4")

  fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
                data = sleepstudy)

  ref_lmer <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy,
                         REML = FALSE)
  expect_loglik_equal(fit, ref_lmer, tol = 1e-6)
  expect_vector_equal(fixef(fit)$mu, lme4::fixef(ref_lmer), tol = 1e-4)

  ref_tmb <- glmmTMB::glmmTMB(Reaction ~ Days + (Days | Subject),
                              sleepstudy, REML = FALSE)
  expect_loglik_equal(fit, ref_tmb, tol = 1e-6)

  se_f <- summary(fit)$coefficients$mu[, "Std. Error"]
  se_g <- summary(ref_tmb)$coefficients$cond[, "Std. Error"]
  expect_vector_equal(se_f, se_g, tol = 1e-4)

  vc_f <- VarCorr(fit)[[1]]
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
  expect_vector_equal(fixef(fit)$mu, lme4::fixef(ref), tol = 1e-4)
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
  dd <- data.frame(y = rnorm(60), x = rnorm(60))
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  ref <- stats::lm(y ~ x, dd)
  expect_loglik_equal(fit, ref, tol = 1e-6)
  expect_vector_equal(fixef(fit)$mu, coef(ref), tol = 1e-5)
})
