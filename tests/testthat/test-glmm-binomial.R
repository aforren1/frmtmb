test_that("cbpp binomial GLMM with trials() matches glmmTMB and glmer", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("glmmTMB")
  data(cbpp, package = "lme4")

  fit <- frm(bf(incidence | trials(size) ~ period + (1 | herd)) +
                  binomial(),
                data = cbpp)
  ref <- glmmTMB::glmmTMB(cbind(incidence, size - incidence) ~
                            period + (1 | herd),
                          data = cbpp, family = binomial)
  expect_loglik_equal(fit, ref, tol = 1e-6)
  expect_vector_equal(fixef(fit)$mu, unname(glmmTMB::fixef(ref)$cond),
                      tol = 1e-4)

  ref2 <- lme4::glmer(cbind(incidence, size - incidence) ~
                        period + (1 | herd),
                      data = cbpp, family = binomial)
  expect_vector_equal(fixef(fit)$mu, lme4::fixef(ref2), tol = 1e-3)
})

test_that("bernoulli responses (0/1 and factor) match glm", {
  set.seed(42)
  dd <- data.frame(x = rnorm(200))
  dd$y <- rbinom(200, 1, plogis(-0.5 + 0.8 * dd$x))
  fit <- frm(bf(y ~ x) + binomial(), data = dd)
  ref <- stats::glm(y ~ x, family = binomial, data = dd)
  expect_loglik_equal(fit, ref, tol = 1e-6)
  expect_vector_equal(fixef(fit)$mu, coef(ref), tol = 1e-5)

  dd$yf <- factor(ifelse(dd$y == 1, "yes", "no"), levels = c("no", "yes"))
  fit2 <- frm(bf(yf ~ x) + binomial(), data = dd)
  expect_vector_equal(fixef(fit2)$mu, coef(ref), tol = 1e-5)
})
