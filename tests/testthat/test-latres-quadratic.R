# The warning that a latent residual's Laplace approximation is biased
# (frmtmb_latres_laplace) read its first-order correction by differences
# of each row's log density. On a lognormal, whose log density is
# quadratic in the linear predictor, the approximation is exact and those
# differences read rounding: it fired on the weight-0 spelling of an
# unread mi() response and not on the mi() spelling, estimates equal to
# 7 digits (dev/quadseam2-findings.md, items 5615 and 5620).

latres_ln_data <- function(seed) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(6L), g = factor(seq_len(10L)))
  d <- d[order(d$g, d$t), ]
  d$x <- stats::rnorm(nrow(d))
  e <- as.vector(replicate(10L, stats::arima.sim(list(ar = 0.6), 6L,
                                                 sd = 0.5)))
  d$y <- exp(0.6 + 0.4 * d$x + e + stats::rnorm(nrow(d), 0, 0.3))
  d$y[sort(sample(nrow(d), 8))] <- NA
  d
}

latres_count <- function(expr) {
  n <- 0L
  v <- withCallingHandlers(expr, frmtmb_latres_laplace = function(w) {
    n <<- n + 1L
    invokeRestart("muffleWarning")
  })
  list(value = v, n = n)
}

test_that("an exact Laplace approximation is not called biased", {
  skip_on_cran()
  d <- latres_ln_data(40401)
  dw <- d
  dw$w <- ifelse(is.na(d$y), 0, 1)
  dw$y[is.na(dw$y)] <- 1
  r1 <- latres_count(frm(bf(y | mi() ~ x + ar(t, g, cov = TRUE)),
                         family = lognormal(), data = d))
  r2 <- latres_count(frm(bf(y | weights(w) ~ x + ar(t, g, cov = TRUE)),
                         family = lognormal(), data = dw))
  expect_identical(c(r1$n, r2$n), c(0L, 0L))
  expect_true(latres_quadratic(r2$value,
                               frame_latres(r2$value$frame)[[1L]]))
})

test_that("a bernoulli latent residual still gets the warning", {
  # p near 0.5 has a third derivative near 0 and a fourth that is not
  skip_on_cran()
  set.seed(502)
  d <- expand.grid(t = seq_len(10L), g = factor(seq_len(20L)))
  d <- d[order(d$g, d$t), ]
  d$x <- stats::rnorm(nrow(d))
  e <- as.vector(replicate(20L, stats::arima.sim(list(ar = 0.5), 10L,
                                                 sd = 0.8)))
  d$yp <- stats::rpois(nrow(d), exp(0.3 * d$x + e))
  d$yb <- stats::rbinom(nrow(d), 1, stats::plogis(0.3 * d$x + e))
  r <- latres_count(frm(bf(yb ~ x + ar(t, g, cov = TRUE)),
                        family = bernoulli(), data = d))
  expect_identical(r$n, 1L)
  expect_false(latres_quadratic(r$value, frame_latres(r$value$frame)[[1L]]))
})

test_that("the quadratic gate does not read a vanishing step as exact", {
  # a bernoulli in 4 groups of 5 whose sderr stops at exp(-10.12): the
  # latent's marginal sd is 7.6e-5, and differences at that step read
  # every row as quadratic, dropping the warning 0.82.0 gave (review of
  # lane quadseam2, F2; dev/reviews/quadseam2-review-latres.R, FEWBERN 3)
  skip_on_cran()
  set.seed(703)
  d <- expand.grid(t = seq_len(5L), g = factor(seq_len(4L)))
  d <- d[order(d$g, d$t), ]
  d$x <- stats::rnorm(nrow(d))
  d$e <- as.vector(replicate(4L, stats::arima.sim(list(ar = 0.5), 5L,
                                                  sd = 0.8)))
  d$yp <- stats::rpois(nrow(d), exp(0.3 * d$x + d$e))
  d$yb <- stats::rbinom(nrow(d), 1, stats::plogis(0.3 * d$x + d$e))
  r <- latres_count(frm(bf(yb ~ x + ar(t, g, cov = TRUE)),
                        family = bernoulli(), data = d))
  expect_identical(r$n, 1L)
  expect_false(latres_quadratic(r$value, frame_latres(r$value$frame)[[1L]]))
})
