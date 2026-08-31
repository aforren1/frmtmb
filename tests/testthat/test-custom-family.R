my_poisson <- function() {
  custom_family(
    "my_poisson",
    dpars = "mu",
    links = list(mu = "log"),
    lpdf = function(y, dpars, aterms) {
      y * log(dpars$mu) - dpars$mu - lgamma(y + 1)
    },
    init_dpars = list(mu = function(y, aterms) mean(y) + 0.1),
    type = "discrete",
    post = list(mean_fn = function(dpars, aterms) dpars$mu),
    sim = function(dpars, aterms, n) stats::rpois(n, dpars$mu)
  )
}

test_that("a hand-written custom family matches the built-in", {
  set.seed(121)
  dd <- data.frame(x = rnorm(300), g = factor(rep(1:15, 20)))
  dd$y <- rpois(300, exp(0.5 + 0.3 * dd$x + rnorm(15, 0, 0.4)[dd$g]))

  fit_c <- frm(bf(y ~ x + (1 | g)) + my_poisson(), data = dd)
  fit_b <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
  expect_lt(abs(as.numeric(logLik(fit_c)) - as.numeric(logLik(fit_b))),
            1e-8)
  expect_vector_equal(fixef(fit_c)$mu, fixef(fit_b)$mu, tol = 1e-6)
  # the whole post-processing stack works on the custom family
  expect_equal(fitted(fit_c), fitted(fit_b), tolerance = 1e-6)
  s <- simulate(fit_c, nsim = 1, seed = 1)
  expect_true(all(s$sim_1 >= 0))
})

test_that("check_custom_family passes a correct lpdf", {
  expect_true(check_custom_family(
    my_poisson(), y = rpois(50, 3),
    dpars = list(mu = rep(2.5, 50))
  ))
})

test_that("check_custom_family catches tape-unsafe code", {
  bad <- custom_family(
    "bad",
    dpars = "mu",
    links = list(mu = "log"),
    lpdf = function(y, dpars, aterms) {
      # base::matrix strips the advector class: values silently wrong
      m <- matrix(dpars$mu, ncol = 1)
      y * log(m[, 1]) - m[, 1] - lgamma(y + 1)
    }
  )
  expect_error(check_custom_family(bad, y = rpois(50, 3),
                                   dpars = list(mu = rep(2.5, 50))),
               "tape|differently|gradient")
})
