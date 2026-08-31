# Comparison helpers against reference packages. Tolerances are absolute.
expect_loglik_equal <- function(fit, ref, tol = 1e-6) {
  expect_lt(abs(as.numeric(stats::logLik(fit)) -
                  as.numeric(stats::logLik(ref))), tol)
}

expect_vector_equal <- function(x, y, tol) {
  expect_equal(length(x), length(y))
  expect_lt(max(abs(unname(x) - unname(y))), tol)
}

# Grouped AR(1) data shared across tests.
sim_ar1_data <- function(seed = 51, n_g = 60, n_t = 6, rho = 0.6,
                         sd_re = 1, sd_res = 0.5) {
  set.seed(seed)
  g <- factor(rep(seq_len(n_g), each = n_t))
  tim <- factor(rep(seq_len(n_t), n_g))
  u <- replicate(n_g, {
    e <- rnorm(n_t, 0, sd_re)
    for (t in 2:n_t) e[t] <- rho * e[t - 1] + sqrt(1 - rho^2) * rnorm(1, 0, sd_re)
    e
  })
  data.frame(y = 1 + as.vector(u) + rnorm(n_g * n_t, 0, sd_res),
             g = g, tim = tim)
}

# Simulated poisson GLMM data shared across tests.
sim_pois_glmm <- function(seed = 101, n_g = 30, n_per = 20) {
  set.seed(seed)
  g <- factor(rep(seq_len(n_g), each = n_per))
  x <- stats::rnorm(n_g * n_per)
  b0 <- stats::rnorm(n_g, 0, 0.5)
  b1 <- stats::rnorm(n_g, 0, 0.3)
  eta <- 0.3 + 0.4 * x + b0[g] + b1[g] * x
  data.frame(y = stats::rpois(length(eta), exp(eta)), x = x, g = g)
}
