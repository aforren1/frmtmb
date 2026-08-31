sim_zi_data <- function(seed = 71, n = 500) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  mu <- exp(0.5 + 0.4 * x)
  zi <- plogis(-0.5 + 0.8 * z)
  data.frame(y = rbinom(n, 1, 1 - zi) * rpois(n, mu), x = x, z = z)
}

test_that("zero_inflated_poisson with zi formula matches glmmTMB", {
  skip_if_not_installed("glmmTMB")
  dd <- sim_zi_data()
  fit <- suppressWarnings(
    frm(bf(y ~ x, zi ~ z) + zero_inflated_poisson(), data = dd)
  )
  ref <- glmmTMB::glmmTMB(y ~ x, ziformula = ~z, family = poisson,
                          data = dd)
  expect_loglik_equal(fit, ref, tol = 1e-6)
  expect_vector_equal(fixef(fit)$mu, unname(glmmTMB::fixef(ref)$cond),
                      tol = 1e-3)
  expect_vector_equal(fixef(fit)$zi, unname(glmmTMB::fixef(ref)$zi),
                      tol = 1e-3)
})

test_that("zi model with random effects in mu matches glmmTMB", {
  skip_if_not_installed("glmmTMB")
  set.seed(72)
  n <- 600
  g <- factor(rep(1:20, 30))
  x <- rnorm(n)
  mu <- exp(0.5 + 0.3 * x + rnorm(20, 0, 0.5)[g])
  dd <- data.frame(y = rbinom(n, 1, 0.75) * rpois(n, mu), x = x, g = g)
  fit <- suppressWarnings(
    frm(bf(y ~ x + (1 | g), zi ~ 1) + zero_inflated_poisson(), data = dd)
  )
  ref <- glmmTMB::glmmTMB(y ~ x + (1 | g), ziformula = ~1,
                          family = poisson, data = dd)
  expect_loglik_equal(fit, ref, tol = 1e-5)
})

test_that("zero_inflated_negbinomial matches glmmTMB", {
  skip_if_not_installed("glmmTMB")
  set.seed(73)
  n <- 600
  x <- rnorm(n)
  mu <- exp(0.8 + 0.4 * x)
  dd <- data.frame(y = rbinom(n, 1, 0.8) *
                     rnbinom(n, size = 1.5, mu = mu), x = x)
  fit <- suppressWarnings(
    frm(bf(y ~ x, zi ~ 1) + zero_inflated_negbinomial(), data = dd)
  )
  ref <- glmmTMB::glmmTMB(y ~ x, ziformula = ~1,
                          family = glmmTMB::nbinom2, data = dd)
  expect_loglik_equal(fit, ref, tol = 1e-5)
  expect_vector_equal(fixef(fit)$mu, unname(glmmTMB::fixef(ref)$cond),
                      tol = 1e-3)
})

test_that("hurdle_poisson matches glmmTMB truncated_poisson hurdle", {
  skip_if_not_installed("glmmTMB")
  set.seed(74)
  n <- 500
  x <- rnorm(n)
  mu <- exp(0.7 + 0.3 * x)
  y <- rpois(n, mu)
  keep0 <- rbinom(n, 1, 0.3)   # 30% structural zeros
  y[keep0 == 1] <- 0
  # resample the sampling zeros into positives for a clean hurdle
  while (any(y == 0 & keep0 == 0)) {
    i <- y == 0 & keep0 == 0
    y[i] <- rpois(sum(i), mu[i])
  }
  dd <- data.frame(y = y, x = x)
  fit <- suppressWarnings(
    frm(bf(y ~ x, hu ~ 1) + hurdle_poisson(), data = dd)
  )
  ref <- glmmTMB::glmmTMB(y ~ x, ziformula = ~1,
                          family = glmmTMB::truncated_poisson, data = dd)
  expect_loglik_equal(fit, ref, tol = 1e-5)
})

test_that("zi simulation respects the mixture", {
  dd <- sim_zi_data(seed = 75)
  fit <- suppressWarnings(
    frm(bf(y ~ x, zi ~ z) + zero_inflated_poisson(), data = dd)
  )
  s <- simulate(fit, nsim = 20, seed = 1)
  expect_lt(abs(mean(as.matrix(s) == 0) - mean(dd$y == 0)), 0.05)
})
