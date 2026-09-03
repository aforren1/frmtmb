# MAP / regularized ML via frm(prior=).

test_that("MAP priors shrink estimates and penalize the objective", {
  set.seed(501)
  dd <- data.frame(x = rnorm(120), g = factor(rep(1:12, 10)))
  dd$y <- rnorm(120, 1 + 0.5 * dd$x + rnorm(12, 0, 0.7)[dd$g], 1)

  ml <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  map <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
             prior = set_prior("normal(0, 0.1)", class = "b"))

  # slope shrunk toward zero, ML untouched
  expect_lt(abs(fixef(map)$mu[["x"]]), abs(fixef(ml)$mu[["x"]]))
  # penalized objective identity at the MAP solution
  nlp <- -stats::dnorm(fixef(map)$mu[["x"]], 0, 0.1, log = TRUE)
  raw_nll <- ml$obj$fn(map$opt$par)
  expect_lt(abs((-as.numeric(logLik(map))) - (raw_nll + nlp)), 1e-6)
  expect_output(print(map), "MAP")

  # near-flat priors reproduce ML
  map2 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
              prior = set_prior("normal(0, 1000)", class = "b"))
  expect_vector_equal(fixef(map2)$mu, fixef(ml)$mu, tol = 1e-3)
})

test_that("an sd prior regularizes a singular variance component", {
  set.seed(502)
  # 4 groups, tiny group effect: ML collapses the variance component
  dd <- data.frame(x = rnorm(80), g = factor(rep(1:4, 20)))
  dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(4, 0, 0.05)[dd$g], 1)
  ml <- suppressWarnings(frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd))
  sd_ml <- sqrt(VarCorr(ml)[[1]][1, 1])

  map <- suppressWarnings(
    frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
        prior = set_prior("exponential(2)", class = "sd"))
  )
  sd_map <- sqrt(VarCorr(map)[[1]][1, 1])
  # regularized away from the boundary, and finite standard errors
  expect_gt(sd_map, sd_ml - 1e-8)
  d <- diagnose(map, quiet = TRUE)
  expect_length(d$bad_se, 0)
})

test_that("a MAP fit's priors carry into frm_sample by default", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  set.seed(503)
  dd <- data.frame(x = rnorm(100))
  dd$y <- rnorm(100, 1 + 0.5 * dd$x, 1)
  map <- frm(bf(y ~ x) + gaussian(), data = dd,
             prior = set_prior("normal(0, 0.05)", class = "b",
                                coef = "x"))
  ds <- suppressWarnings(frm_sample(map, chains = 1, iter = 500,
                                    refresh = 0, seed = 1))
  # judged against the chain's own spread: a seeded chain is not
  # platform-deterministic, and this asserts wiring, not mixing
  mx <- as.matrix(ds)[, "x"]
  if (sampler_gates_on()) {
    expect_lt(abs(mean(mx)), 5 * stats::sd(mx) + 1e-8)
  }
})