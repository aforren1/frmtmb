test_that("emmeans marginal means match glmmTMB", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("glmmTMB")
  set.seed(131)
  n <- 300
  dd <- data.frame(
    f = factor(sample(c("a", "b", "c"), n, replace = TRUE)),
    x = rnorm(n),
    g = factor(rep(1:15, length.out = n))
  )
  dd$y <- rnorm(n, 1 + c(a = 0, b = 1, c = -1)[dd$f] + 0.5 * dd$x +
                  rnorm(15, 0, 0.7)[dd$g], 1)

  fit <- frm(bf(y ~ f + x + (1 | g)) + gaussian(), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ f + x + (1 | g), data = dd, REML = FALSE)

  em_f <- as.data.frame(emmeans::emmeans(fit, "f"))
  em_g <- as.data.frame(emmeans::emmeans(ref, "f"))
  expect_vector_equal(em_f$emmean, em_g$emmean, tol = 1e-4)
  expect_vector_equal(em_f$SE, em_g$SE, tol = 1e-4)

  ct_f <- as.data.frame(emmeans::contrast(emmeans::emmeans(fit, "f"),
                                          "pairwise"))
  expect_length(ct_f$estimate, 3)
})

test_that("as_tmbstan hands the objective to NUTS", {
  skip_if_not_installed("tmbstan")
  set.seed(132)
  dd <- data.frame(x = rnorm(80))
  dd$y <- rnorm(80, 1 + 0.5 * dd$x, 1)
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  sf <- suppressWarnings(as_tmbstan(fit, chains = 1, iter = 400,
                                    refresh = 0, seed = 1))
  expect_s4_class(sf, "stanfit")
  skip_if_not_installed("rstan")
  dr <- rstan::extract(sf, "beta")$beta
  expect_lt(abs(mean(dr[, 1]) - fixef(fit)$mu[[1]]), 0.2)
})
