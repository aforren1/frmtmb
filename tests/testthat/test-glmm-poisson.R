test_that("poisson GLMM with correlated slopes matches glmmTMB", {
  skip_if_not_installed("glmmTMB")
  dd <- sim_pois_glmm()

  fit <- frm(bf(y ~ x + (x | g)) + poisson(), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ x + (x | g), dd, family = poisson)

  expect_loglik_equal(fit, ref, tol = 1e-6)
  expect_vector_equal(fixef(fit)$mu, unname(glmmTMB::fixef(ref)$cond),
                      tol = 1e-4)
  se_f <- summary(fit)$coefficients$mu[, "Std. Error"]
  se_g <- summary(ref)$coefficients$cond[, "Std. Error"]
  expect_vector_equal(se_f, se_g, tol = 1e-4)
})

test_that("poisson GLM (no random effects) matches glm", {
  dd <- sim_pois_glmm(n_g = 1, n_per = 200)
  fit <- frm(bf(y ~ x) + poisson(), data = dd)
  ref <- stats::glm(y ~ x, family = poisson, data = dd)
  expect_loglik_equal(fit, ref, tol = 1e-6)
  expect_vector_equal(fixef(fit)$mu, coef(ref), tol = 1e-5)
})

test_that("weights() matches weighted glm", {
  dd <- sim_pois_glmm(n_g = 1, n_per = 150)
  dd$w <- runif(nrow(dd), 0.5, 2)
  fit <- frm(bf(y | weights(w) ~ x) + poisson(), data = dd)
  ref <- suppressWarnings(
    stats::glm(y ~ x, family = poisson, data = dd, weights = w)
  )
  expect_vector_equal(fixef(fit)$mu, coef(ref), tol = 1e-5)
})

test_that("offset() is applied", {
  dd <- sim_pois_glmm(n_g = 1, n_per = 150)
  dd$expo <- runif(nrow(dd), 0.5, 2)
  fit <- frm(bf(y ~ x + offset(log(expo))) + poisson(), data = dd)
  ref <- stats::glm(y ~ x + offset(log(expo)), family = poisson, data = dd)
  expect_vector_equal(fixef(fit)$mu, coef(ref), tol = 1e-5)
})
