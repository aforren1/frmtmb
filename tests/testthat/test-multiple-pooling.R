# Pooled variance components and hypothesis tests for frm_multiple,
# plus mice interop. Validation leans on the degenerate identity:
# m identical imputations have zero between-variance, so every pooled
# quantity must collapse to its single-fit counterpart.

test_that("identical imputations reduce pooling to the single fit", {
  set.seed(61)
  n_g <- 15
  n_per <- 8
  g <- factor(rep(seq_len(n_g), each = n_per))
  x <- rnorm(n_g * n_per)
  b0 <- rnorm(n_g, 0, 0.7)
  b1 <- rnorm(n_g, 0, 0.4)
  dd <- data.frame(
    y = 1 + 0.5 * x + b0[g] + b1[g] * x + rnorm(n_g * n_per, 0, 0.6),
    x = x, g = g
  )
  form <- bf(y ~ x + (1 + x | g)) + gaussian()
  mfit <- frm_multiple(form, data = list(dd, dd, dd))
  f0 <- frm(form, data = dd)
  cv <- confint_varcorr(f0)
  pvc <- mfit$pooled_varcorr

  expect_equal(pvc$grp, cv$block)
  expect_equal(pvc$term, cv$term)
  expect_equal(pvc$type, cv$type)
  # between-variance is exactly zero, so estimates match the single fit
  expect_vector_equal(pvc$estimate, cv$estimate, tol = 1e-8)
  expect_true(all(pvc$fmi < 1e-8))
  # intervals differ only by the t-vs-z quantile: on the transformed
  # scale the half-width ratio is qt(.975, df) / qnorm(.975) exactly
  tr <- function(type, v) {
    out <- numeric(length(v))
    i <- type == "cor"
    out[i] <- atanh(v[i])
    out[!i] <- log(v[!i])
    out
  }
  hw_p <- tr(pvc$type, pvc$upr) - tr(pvc$type, pvc$estimate)
  hw_s <- tr(cv$type, cv$upr) - tr(cv$type, cv$estimate)
  expect_vector_equal(hw_p / hw_s,
                      stats::qt(0.975, pvc$df) / stats::qnorm(0.975),
                      tol = 1e-6)
  # hence the pooled (t) interval contains the single-fit (z) interval
  expect_true(all(pvc$lwr <= cv$lwr & pvc$upr >= cv$upr))

  hyps <- c("x = 0", "sd_g__Intercept - 0.5", "cor_g__Intercept__x")
  hp <- hypothesis(mfit, hyps)
  h0 <- hypothesis(f0, hyps)
  expect_s3_class(hp, "frmtmb_hypothesis")
  expect_vector_equal(hp$estimate, h0$estimate, tol = 1e-8)
  expect_vector_equal(hp$se, h0$se, tol = 1e-8)
  expect_vector_equal(hp$t, h0$z, tol = 1e-6)
  # t reference is heavier-tailed than the single fit's normal
  expect_true(all(hp$p >= h0$p - 1e-12))
  expect_true(all(hp$lwr <= h0$lwr & hp$upr >= h0$upr))

  expect_output(print(mfit), "variance components")
})

test_that("pooling across distinct imputations is Rubin-consistent", {
  set.seed(62)
  n_g <- 20
  n_per <- 6
  g <- factor(rep(seq_len(n_g), each = n_per))
  x <- rnorm(n_g * n_per)
  b0 <- rnorm(n_g, 0, 0.8)
  y <- 1 + 0.5 * x + b0[g] + rnorm(n_g * n_per, 0, 0.7)
  x_mis <- x
  x_mis[sample(n_g * n_per, 25)] <- NA
  imps <- lapply(1:4, function(i) {
    xi <- x_mis
    xi[is.na(xi)] <- rnorm(sum(is.na(xi)), mean(x_mis, na.rm = TRUE),
                           stats::sd(x_mis, na.rm = TRUE))
    data.frame(y = y, x = xi, g = g)
  })
  mfit <- frm_multiple(bf(y ~ x + (1 | g)) + gaussian(), data = imps)
  pvc <- mfit$pooled_varcorr
  expect_equal(nrow(pvc), 1L)
  expect_equal(pvc$type, "sd")
  # log-scale mean: the pooled sd is the geometric mean of the
  # per-fit sds, so it lies inside their range
  ests <- vapply(mfit$fits,
                 function(f) confint_varcorr(f)$estimate[1], numeric(1))
  expect_gte(pvc$estimate, min(ests))
  expect_lte(pvc$estimate, max(ests))
  expect_true(pvc$fmi >= 0 && pvc$fmi <= 1)
  expect_true(pvc$lwr < pvc$estimate && pvc$estimate < pvc$upr)

  h <- "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)"
  hp <- hypothesis(mfit, h)
  per <- lapply(mfit$fits, function(f) hypothesis(f, h))
  # qbar is the plain mean; the total variance adds a nonnegative
  # between component to the mean within variance
  expect_vector_equal(hp$estimate,
                      mean(vapply(per, `[[`, numeric(1), "estimate")),
                      tol = 1e-8)
  ubar <- mean(vapply(per, `[[`, numeric(1), "se")^2)
  expect_gte(hp$se^2, ubar - 1e-12)
  expect_true(hp$df > 0 && hp$p >= 0 && hp$p <= 1)
})

test_that("mids input matches mice::pool over equivalent lm fits", {
  skip_if_not_installed("mice")
  set.seed(73)
  n <- 200
  x <- rnorm(n)
  z <- rnorm(n, 0.4 * x, 1)
  y <- rnorm(n, 1 + 0.5 * x + 0.3 * z, 1)
  dd <- data.frame(y = y, x = x, z = z)
  # MAR: missingness driven by the always-observed response
  dd$x[runif(n) < plogis(-1.2 + 0.5 * y)] <- NA
  dd$z[runif(n) < 0.15] <- NA
  imp <- mice::mice(dd, m = 5, seed = 11, printFlag = FALSE)
  mfit <- frm_multiple(bf(y ~ x + z) + gaussian(), data = imp)
  expect_length(mfit$fits, 5L)
  ref <- summary(mice::pool(with(imp, lm(y ~ x + z))))
  est <- mfit$pooled[1:3, ]   # mu coefficients, lm term order
  # Both sides run the same Rubin formulas, and gaussian ML beta-hat
  # is OLS, so pooled estimates agree essentially exactly. SEs differ
  # by the residual-variance convention (ML divides by n, lm by n - p,
  # scaling the within part by ~ n / (n - p)), so compare loosely.
  expect_vector_equal(est$estimate, ref$estimate, tol = 1e-4)
  expect_lt(max(abs(est$se / ref$std.error - 1)), 0.02)
})
