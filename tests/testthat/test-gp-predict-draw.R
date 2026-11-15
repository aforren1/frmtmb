# predict() on a gp() field draws its coefficients with the other estimates,
# and the field at a new position given the drawn field, as brms's draws
# carry them (R/predict-brms.R predict_b_drawer()). The brms table is in
# dev/gpfield-findings.md, 1919.

gpd_data <- function() {
  set.seed(9200)
  n <- 30
  d <- data.frame(x = sort(stats::runif(n, 0, 10)))
  d$y <- 2 * sin(d$x) + stats::rnorm(n, 0, 0.2)
  d
}

test_that("predict() draws a gp() field past its positions", {
  d <- gpd_data()
  fit <- frm(bf(y ~ gp(x)), data = d, family = gaussian())
  sdgp <- exp(fit$estimates$theta[1L])
  sigma <- exp(fit$estimates$betad[1L])
  nd <- data.frame(x = c(60, d$x[5L]))
  set.seed(1)
  p <- predict(fit, newdata = nd, summary = FALSE, ndraws = 2000)
  # far from every position the field is its prior: sd near
  # sqrt(sdgp^2 + sigma^2) and not sigma alone, which the mode gave
  s_far <- stats::sd(p[, 1L])
  prior <- sqrt(sdgp^2 + sigma^2)
  expect_gt(s_far, 0.8 * prior)
  expect_gt(s_far, 2 * sigma)
  # at a position the fit saw, the field is pinned by the data: the spread
  # is about the noise's. A field read through each replicate's own sd at
  # the mode's z spread 1.93, ten times the noise
  expect_lt(stats::sd(p[, 2L]), 0.5 * s_far)
  expect_lt(stats::sd(p[, 2L]), 3 * sigma)
  # a population term: re_formula = NA draws it too
  set.seed(1)
  p0 <- predict(fit, newdata = nd, summary = FALSE, ndraws = 2000,
                re_formula = NA)
  expect_gt(stats::sd(p0[, 1L]), 0.8 * prior)
})

test_that("an approximate field's coefficients are drawn", {
  d <- gpd_data()
  fit <- frm(bf(y ~ gp(x, k = 15)), data = d, family = gaussian())
  sigma <- exp(fit$estimates$betad[1L])
  set.seed(2)
  p <- predict(fit, newdata = data.frame(x = 5.5), summary = FALSE,
               ndraws = 2000)
  set.seed(2)
  q <- predict(fit, newdata = data.frame(x = 5.5), summary = FALSE,
               ndraws = 2000, propagate_error = FALSE)
  # inside the positions the field and the intercept are drawn together:
  # a little more than the noise alone. Holding the field at its mode while
  # the intercept, which it can absorb, was drawn on its own read 78 times
  # the noise's variance
  expect_gt(stats::var(p[, 1L]), stats::var(q[, 1L]))
  expect_lt(stats::var(p[, 1L]), 3 * stats::var(q[, 1L]))
})
