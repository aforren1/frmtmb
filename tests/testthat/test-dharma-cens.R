# dharma_residuals() on a cens() fit. A censored row observes an event,
# not a value, so its residual is a uniform draw between the predictive
# distribution function at the two ends of the event (Dunn and Smyth's
# randomized quantile residual, extended to censoring). Through 0.81.0
# the function scored a censored row as the exact value it recorded,
# against the latent draws: a right-censored row read F(y), the bottom
# of its interval, every time, and nothing was said.
#
# The assertions are exact: each residual is checked against the
# interval its own row's draws give, with the draws the function itself
# makes (simulate() with the same seed), so no tolerance is involved.

cens_gauss <- function(seed = 41, n = 300) {
  set.seed(seed)
  x <- stats::rnorm(n)
  ys <- 1 + 0.5 * x + stats::rnorm(n)
  # row-varying censoring times on both sides, and intervals
  rc <- stats::runif(n, 0.8, 2.5)
  lc <- stats::runif(n, -1, 0.3)
  cen <- ifelse(ys > rc, 1, ifelse(ys < lc, -1, 0))
  y <- ifelse(cen == 1, rc, ifelse(cen == -1, lc, ys))
  hi <- y
  iv <- which(cen == 0)[seq(1, sum(cen == 0), by = 6)]
  cen[iv] <- 2
  y[iv] <- floor(ys[iv] * 2) / 2
  hi[iv] <- y[iv] + 0.5
  data.frame(x = x, y = y, cen = cen, hi = hi)
}

# the PIT interval of each row from its own draws
row_interval <- function(sims, lo_end, hi_end) {
  lo <- rowMeans(sims <= lo_end)
  hi <- rowMeans(sims <= hi_end)
  list(lo = lo, hi = hi)
}

test_that("a censored gaussian row's residual is drawn over its event", {
  skip_if_not_installed("DHARMa")
  d <- cens_gauss()
  expect_true(all(c(-1, 1, 2) %in% d$cen))
  fit <- frm(bf(y | cens(cen, hi) ~ x), family = gaussian(), data = d)
  res <- dharma_residuals(fit, nsim = 200, seed = 7)
  r <- res$scaledResiduals
  sims <- as.matrix(simulate(fit, nsim = 200, seed = 7))
  # right: Y > y
  i <- which(d$cen == 1)
  iv <- row_interval(sims[i, ], d$y[i], Inf)
  ok <- iv$lo < 1
  expect_gt(sum(ok), 10L)
  expect_true(all(r[i][ok] > iv$lo[ok] & r[i][ok] <= 1))
  # left: Y <= y
  i <- which(d$cen == -1)
  iv <- row_interval(sims[i, ], -Inf, d$y[i])
  ok <- iv$hi > 0
  expect_gt(sum(ok), 10L)
  expect_true(all(r[i][ok] >= 0 & r[i][ok] < iv$hi[ok]))
  # interval: y < Y <= y2
  i <- which(d$cen == 2)
  iv <- row_interval(sims[i, ], d$y[i], d$hi[i])
  ok <- iv$lo < iv$hi
  expect_gt(sum(ok), 10L)
  expect_true(all(r[i][ok] > iv$lo[ok] & r[i][ok] < iv$hi[ok]))
  # an uncensored row keeps F(y), as on a fit without cens()
  i <- which(d$cen == 0)
  expect_identical(r[i], unname(rowMeans(sims[i, ] < d$y[i])))
  # the observation of a censored row is a point of its event
  i <- which(d$cen == 1)
  expect_true(all(res$observedResponse[i] > d$y[i]))
})

test_that("a right-censored count is read on the lattice", {
  skip_if_not_installed("DHARMa")
  set.seed(43)
  n <- 300
  x <- stats::rnorm(n)
  ys <- stats::rpois(n, exp(1 + 0.4 * x))
  k <- sample(1:4, n, TRUE)
  cen <- as.numeric(ys > k)
  y <- ifelse(cen == 1, k, ys)
  d <- data.frame(x = x, y = y, cen = cen)
  fit <- frm(bf(y | cens(cen) ~ x), family = poisson(), data = d)
  res <- dharma_residuals(fit, nsim = 200, seed = 9)
  r <- res$scaledResiduals
  sims <- as.matrix(simulate(fit, nsim = 200, seed = 9))
  # "censored at k" is Y > k: a draw AT k is below the event, so the
  # residual lies above F(k), not between F(k - 1) and F(k)
  i <- which(d$cen == 1)
  lo <- rowMeans(sims[i, ] <= d$y[i])
  ok <- lo < 1
  expect_gt(sum(ok), 10L)
  expect_true(all(r[i][ok] > lo[ok] & r[i][ok] <= 1))
  # an uncensored count keeps its ordinary randomized residual
  i <- which(d$cen == 0)
  lo0 <- rowMeans(sims[i, ] < d$y[i])
  hi0 <- rowMeans(sims[i, ] <= d$y[i])
  expect_true(all(r[i] >= lo0 & r[i] <= hi0))
})

test_that("a fit without a censored row is untouched", {
  skip_if_not_installed("DHARMa")
  d <- cens_gauss()
  d$cen <- 0
  d$y <- d$hi
  a <- frm(bf(y | cens(cen) ~ x), family = gaussian(), data = d)
  ra <- dharma_residuals(a, nsim = 100, seed = 3)
  sims <- as.matrix(simulate(a, nsim = 100, seed = 3))
  expect_identical(ra$scaledResiduals, unname(rowMeans(sims < d$y)))
  expect_identical(unname(ra$observedResponse), d$y)
  expect_identical(unname(ra$simulatedResponse), unname(sims))
  # the absent case of the refusal below: no censored row, so DHARMa's
  # own options pass through
  expect_s3_class(dharma_residuals(a, nsim = 50, seed = 3,
                                   method = "traditional"), "DHARMa")
})

test_that("a censored fit refuses DHARMa options that break its residual", {
  skip_if_not_installed("DHARMa")
  d <- cens_gauss()
  fit <- frm(bf(y | cens(cen, hi) ~ x), family = gaussian(), data = d)
  # passed through, both gave residuals far from uniform on a correct
  # censored model (KS p 2.3e-25 and 2.3e-14 where PIT gives 0.86)
  expect_error(dharma_residuals(fit, nsim = 50, seed = 1,
                                method = "traditional"),
               "method = \"traditional\") is refused on a cens() fit",
               fixed = TRUE)
  expect_error(dharma_residuals(fit, nsim = 50, seed = 1, method = "trad"),
               "is refused on a cens() fit", fixed = TRUE)
  expect_error(dharma_residuals(fit, nsim = 50, seed = 1,
                                rotation = "estimated"),
               "rotation = ) is refused on a cens() fit", fixed = TRUE)
  expect_s3_class(dharma_residuals(fit, nsim = 50, seed = 1,
                                   method = "PIT"), "DHARMa")
})
