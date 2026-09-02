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

test_that("a custom lpdf needs no ADoverload boilerplate of its own", {
  # numeric-first c() and pad[i] <- on advectors, no ADoverload lines:
  # frmtmb_family() splices the overloads onto the body
  cf <- custom_family(
    "gauss_bare", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      pad <- c(0, dpars$mu)
      pad[1] <- pad[2]
      RTMB::dnorm(y, pad[-1], dpars$sigma, log = TRUE)
    })
  # the wrap is visible on the stored function, and a function that
  # binds the overloads itself is left untouched
  expect_true("ADoverload" %in% all.names(body(cf$lpdf)))
  own <- function(y, dpars, aterms) {
    "c" <- RTMB::ADoverload("c")
    RTMB::dnorm(y, dpars$mu, dpars$sigma, log = TRUE)
  }
  cf2 <- custom_family("gauss_own", dpars = c("mu", "sigma"),
                       links = list(mu = "identity", sigma = "log"),
                       lpdf = own)
  expect_identical(cf2$lpdf, own)
  set.seed(2)
  dd <- data.frame(x = stats::rnorm(120))
  dd$y <- stats::rnorm(120, 1 + 0.5 * dd$x, 1.3)
  fc <- frm(bf(y ~ x), family = cf, data = dd)
  ref <- frm(bf(y ~ x), family = gaussian(), data = dd)
  expect_equal(as.numeric(logLik(fc)), as.numeric(logLik(ref)),
               tolerance = 1e-6)
  expect_equal(unname(unlist(fixef(fc)$mu)),
               unname(unlist(fixef(ref)$mu)), tolerance = 1e-4)
})
