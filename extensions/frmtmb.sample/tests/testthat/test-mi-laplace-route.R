# frm_sample(laplace = TRUE) on a missing value no density reads: the formula
# and the fit take one route, which leaves the row out as frm() does. The same
# model with the rows deleted is the reference; the error of the old route was
# many posterior standard deviations (0.910 against 0.359 on sigma).

skip_on_cran()

route_data <- function() {
  set.seed(10600)
  n <- 200
  x <- rnorm(n)
  d <- data.frame(x = x, y = 0.3 + 0.5 * x + rnorm(n, 0, 0.5) + rexp(n, 1 / 2))
  m <- sample(n, 60)
  d$y[m] <- NA
  list(d = d, m = m)
}
route_prior <- function() {
  set_prior("normal(0, 5)", class = "b") +
    set_prior("normal(0, 5)", class = "Intercept") +
    set_prior("exponential(1)", class = "sigma") +
    set_prior("exponential(0.5)", class = "beta")
}
route_means <- function(s) {
  dm <- as.matrix(as_draws_matrix(s))
  colMeans(dm[, c("sigma", "beta"), drop = FALSE])
}
route_go <- function(x, ...) {
  suppressWarnings(suppressMessages(frm_sample(
    x, prior = route_prior(), chains = 2, iter = 800, seed = 1, refresh = 0,
    cores = 1, ...)))
}

test_that("the formula's laplace = TRUE reads the rows the fit reads", {
  skip_sampler()
  r <- route_data()
  a <- route_go(bf(y | mi() ~ x), family = exgaussian(), data = r$d,
                laplace = TRUE)
  b <- route_go(frm(bf(y | mi() ~ x), family = exgaussian(), data = r$d),
                laplace = TRUE)
  ref <- route_go(bf(y ~ x), family = exgaussian(), data = r$d[-r$m, ])
  dm <- as.matrix(as_draws_matrix(ref))
  sds <- apply(dm[, c("sigma", "beta"), drop = FALSE], 2L, stats::sd)
  # in the reference's own posterior sd: the two chains' means, at 800 draws,
  # are within a third of that; the old route sat at 5 of them
  expect_lt(max(abs(route_means(a) - route_means(ref)) / sds), 1 / 3)
  expect_lt(max(abs(route_means(b) - route_means(ref)) / sds), 1 / 3)
})

test_that("a read value with no mode is refused on the formula's route too", {
  set.seed(12200)
  n <- 80
  d <- data.frame(x = rnorm(n))
  d$y <- 0.3 * d$x + rnorm(n) - rnorm(n)
  d$z <- 0.5 + 0.6 * d$y + rnorm(n)
  d$y[1:8] <- NA
  expect_error(
    frm_sample(bf(y | mi() ~ x) + asym_laplace() + bf(z ~ mi(y)) + gaussian(),
               data = d, laplace = TRUE, chains = 1, iter = 10, refresh = 0),
    "no such mode")
})
