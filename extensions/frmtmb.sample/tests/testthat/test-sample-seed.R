# frm_sample(seed = ) fixes every chain, as brm(seed = ) does. Chain 1
# starts at the mode and later chains at the mode plus a jitter drawn in
# R; through frmtmb.sample 0.19.0 the jitter came from the session's
# stream, so two calls with seed = 3 and 2 chains differed unless
# set.seed() came first (dev/dblbar-rev2-seed2.R, lane dblbar).

skip_on_cran()

sps_fit <- function() {
  set.seed(20261014)
  d <- data.frame(x = stats::rnorm(60), g = factor(rep(1:6, 10)))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(6, 0, 0.5)[d$g] + stats::rnorm(60)
  frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
}

sps_draws <- function(fit, ...) {
  suppressWarnings(suppressMessages(as.matrix(
    frm_sample(fit, chains = 2, iter = 200, refresh = 0, ...))))
}

test_that("two calls with the same seed give the same chains", {
  skip_sampler()
  fit <- sps_fit()
  # the session's stream moves between the calls, which is the case a
  # seed has to survive
  set.seed(1)
  a <- sps_draws(fit, seed = 3)
  stats::runif(5)
  b <- sps_draws(fit, seed = 3)
  expect_identical(a, b)
  # and a different seed gives different chains, so the identity above
  # is the seed's and not a sampler that ignores it
  c3 <- sps_draws(fit, seed = 4)
  expect_false(identical(a, c3))
})

test_that("the seeded jitter leaves the session's random stream alone", {
  # The whole call still moves it: rstan::sampling() draws each chain's
  # permutation of its kept draws with sample.int() whatever the seed
  # (dev/sampleprior-rngprobe.R), and so does brm(seed = ) through the
  # same function. The jitter is frm_sample()'s own draw, and it does not
  set.seed(7)
  before <- .Random.seed
  v <- frmtmb.sample:::with_sample_seed(3, stats::rnorm(2))
  expect_identical(.Random.seed, before)
  set.seed(3)
  expect_identical(v, stats::rnorm(2))
  # and with no stream yet, none is left behind
  rm(".Random.seed", envir = globalenv())
  invisible(frmtmb.sample:::with_sample_seed(3, stats::rnorm(1)))
  expect_false(exists(".Random.seed", envir = globalenv(),
                      inherits = FALSE))
  set.seed(7)
})

test_that("without a seed the jitter is the session's, as brms leaves it", {
  skip_sampler()
  fit <- sps_fit()
  # chain 1 is the mode and the jitter of chain 2 is all that differs
  # between two inits; set.seed() decides it when no seed is given
  ini <- function() {
    frmtmb.sample:::with_sample_seed(NULL, frmtmb.sample:::mode_inits(
      c(1, 2, 3), 2, 0.25))
  }
  set.seed(11)
  i1 <- ini()
  set.seed(11)
  i2 <- ini()
  set.seed(12)
  i3 <- ini()
  expect_identical(i1, i2)
  expect_false(identical(i1[[2]], i3[[2]]))
  expect_identical(i1[[1]], c(1, 2, 3))
  # a seed of NA is brms's "no seed" and reads the session's stream too
  set.seed(11)
  expect_identical(frmtmb.sample:::with_sample_seed(NA, stats::runif(2)),
                   {
                     set.seed(11)
                     stats::runif(2)
                   })
  expect_error(frmtmb.sample:::with_sample_seed("a", 1),
               "one finite number")
})
