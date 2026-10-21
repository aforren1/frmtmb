# frm_sample() from a formula, with every row owning a random effect in
# mu and weights summing to n + 1 or more: under brms's default prior on
# sigma, finite at 0, the posterior does not integrate. Measured before
# the refusal: 2 chains of 1000 gave sigma between 5.7e-05 and 1.1e-04
# with a largest Rhat of 2.5. frm() refuses the fit route already.

wus_data <- function(seed = 20261021) {
  set.seed(seed)
  n <- 60
  d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n),
                  w = runif(n, 0.5, 2))
  d$y <- 0.5 + 0.4 * d$x + rnorm(n, 0, 0.5) + rnorm(n, 0, 0.8)
  d
}

test_that("the formula route refuses an improper posterior", {
  skip_sampler()
  d <- wus_data()
  expect_gt(sum(d$w) - nrow(d), 1)
  allow_warnings(
    expect_error(
      frm_sample(bf(y | weights(w) ~ x + (1 | obs)), family = gaussian(),
                 data = d, chains = 1, iter = 200, refresh = 0, seed = 1),
      "posterior of this model is improper"),
    "gives every observation its own random effect")
})

test_that("the formula route samples when the weights sum below n", {
  skip_on_cran()
  skip_sampler()
  d <- wus_data()
  d$w <- d$w * 0.5
  expect_lt(sum(d$w), nrow(d))
  ds <- suppressWarnings(suppressMessages(
    frm_sample(bf(y | weights(w) ~ x + (1 | obs)), family = gaussian(),
               data = d, chains = 1, iter = 200, refresh = 0, seed = 1)))
  expect_s3_class(ds, "frmtmb_draws")
})
