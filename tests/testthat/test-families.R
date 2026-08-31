test_that("Gamma matches glmmTMB", {
  skip_if_not_installed("glmmTMB")
  set.seed(21)
  n <- 400
  x <- rnorm(n)
  mu <- exp(0.5 + 0.3 * x)
  y <- rgamma(n, shape = 2, scale = mu / 2)
  dd <- data.frame(y, x)
  fit <- frm(bf(y ~ x) + Gamma(link = "log"), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ x, family = Gamma(link = "log"), data = dd)
  expect_loglik_equal(fit, ref, tol = 1e-6)
  expect_vector_equal(fixef(fit)$mu, unname(glmmTMB::fixef(ref)$cond),
                      tol = 1e-4)
})

test_that("lognormal equals lm on the log scale plus the Jacobian", {
  set.seed(22)
  n <- 300
  x <- rnorm(n)
  y <- exp(rnorm(n, 0.5 + 0.7 * x, 0.4))
  dd <- data.frame(y, x)
  fit <- frm(bf(y ~ x) + lognormal(), data = dd)
  ref <- lm(log(y) ~ x, dd)
  expect_lt(abs(as.numeric(logLik(fit)) -
                  (as.numeric(logLik(ref)) - sum(log(y)))), 1e-6)
  expect_vector_equal(fixef(fit)$mu, coef(ref), tol = 1e-5)
})

test_that("student matches MASS::fitdistr", {
  skip_if_not_installed("MASS")
  set.seed(23)
  y <- 3 + 1.5 * rt(800, df = 5)
  dd <- data.frame(y = y)
  fit <- frm(bf(y ~ 1) + student(), data = dd)
  ref <- suppressWarnings(MASS::fitdistr(y, "t"))
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))), 1e-4)
  expect_lt(abs(fixef(fit)$mu[[1]] - ref$estimate[["m"]]), 1e-3)
  est_sigma <- exp(fixef(fit)$sigma[[1]])
  expect_lt(abs(est_sigma - ref$estimate[["s"]]), 1e-3)
  est_nu <- 1 + exp(fixef(fit)$nu[[1]])
  expect_lt(abs(est_nu - ref$estimate[["df"]]) / ref$estimate[["df"]], 1e-2)
})

test_that("negbinomial GLM matches MASS::glm.nb", {
  skip_if_not_installed("MASS")
  set.seed(24)
  n <- 500
  x <- rnorm(n)
  y <- rnbinom(n, size = 1.5, mu = exp(0.5 + 0.4 * x))
  dd <- data.frame(y, x)
  fit <- frm(bf(y ~ x) + negbinomial(), data = dd)
  ref <- MASS::glm.nb(y ~ x, data = dd)
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))), 1e-5)
  expect_vector_equal(fixef(fit)$mu, coef(ref), tol = 1e-4)
  expect_lt(abs(exp(fixef(fit)$shape[[1]]) - ref$theta), 1e-2)
})

test_that("nbinom1 matches glmmTMB", {
  skip_if_not_installed("glmmTMB")
  set.seed(25)
  n <- 500
  x <- rnorm(n)
  mu <- exp(1 + 0.3 * x)
  phi <- 1.5
  y <- rnbinom(n, size = mu / phi, mu = mu)
  dd <- data.frame(y, x)
  fit <- frm(bf(y ~ x) + nbinom1(), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ x, family = glmmTMB::nbinom1, data = dd)
  expect_loglik_equal(fit, ref, tol = 1e-5)
  expect_vector_equal(fixef(fit)$mu, unname(glmmTMB::fixef(ref)$cond),
                      tol = 1e-4)
})

test_that("beta matches glmmTMB beta_family", {
  skip_if_not_installed("glmmTMB")
  set.seed(26)
  n <- 400
  x <- rnorm(n)
  mu <- plogis(0.3 + 0.6 * x)
  phi <- 8
  y <- rbeta(n, mu * phi, (1 - mu) * phi)
  dd <- data.frame(y, x)
  fit <- frm(bf(y ~ x) + Beta(), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ x, family = glmmTMB::beta_family(),
                          data = dd)
  expect_loglik_equal(fit, ref, tol = 1e-6)
  expect_vector_equal(fixef(fit)$mu, unname(glmmTMB::fixef(ref)$cond),
                      tol = 1e-4)
  expect_lt(abs(exp(fixef(fit)$phi[[1]]) - glmmTMB::sigma(ref)), 1e-2)
})

test_that("tweedie matches glmmTMB", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("mgcv")
  set.seed(27)
  n <- 400
  x <- rnorm(n)
  mu <- exp(0.8 + 0.3 * x)
  y <- mgcv::rTweedie(mu, p = 1.5, phi = 1.2)
  dd <- data.frame(y, x)
  fit <- frm(bf(y ~ x) + tweedie(), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ x, family = glmmTMB::tweedie(), data = dd)
  expect_loglik_equal(fit, ref, tol = 1e-4)
  expect_vector_equal(fixef(fit)$mu, unname(glmmTMB::fixef(ref)$cond),
                      tol = 1e-3)
})

test_that("compois fits and recovers the mean model", {
  set.seed(28)
  n <- 120
  x <- rnorm(n)
  y <- rpois(n, exp(0.6 + 0.4 * x))   # nu = 1 reduces to poisson
  dd <- data.frame(y, x)
  fit <- frm(bf(y ~ x) + compois(), data = dd)
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_vector_equal(fixef(fit)$mu, c(0.6, 0.4), tol = 0.3)
  # poisson nested in compois: likelihood at optimum can't be worse
  ref <- glm(y ~ x, family = poisson, data = dd)
  expect_gt(as.numeric(logLik(fit)), as.numeric(logLik(ref)) - 1e-6)
})

test_that("simulate() round-trips through family simulators", {
  set.seed(29)
  dd <- data.frame(x = rnorm(300), g = factor(rep(1:10, 30)))
  dd$y <- rpois(300, exp(0.5 + 0.3 * dd$x + rnorm(10, 0, 0.4)[dd$g]))
  fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)

  s1 <- simulate(fit, nsim = 3, seed = 1)
  s2 <- simulate(fit, nsim = 3, seed = 1)
  expect_identical(s1, s2)
  expect_identical(dim(s1), c(300L, 3L))
  expect_true(all(s1$sim_1 >= 0))

  # marginal simulation redraws random effects
  m1 <- simulate(fit, nsim = 1, seed = 2, re.form = NA)
  c1 <- simulate(fit, nsim = 1, seed = 2)
  expect_false(identical(m1, c1))
})
