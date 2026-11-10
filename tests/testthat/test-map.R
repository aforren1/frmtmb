# MAP / regularized ML via frm(prior=).

test_that("MAP priors shrink estimates and penalize the objective", {
  set.seed(501)
  dd <- data.frame(x = rnorm(120), g = factor(rep(1:12, 10)))
  dd$y <- rnorm(120, 1 + 0.5 * dd$x + rnorm(12, 0, 0.7)[dd$g], 1)

  ml <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  map <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
             prior = set_prior("normal(0, 0.1)", class = "b"))

  # slope shrunk toward zero, ML untouched
  expect_lt(abs(fixef_by_dpar(map)$mu[["x"]]), abs(fixef_by_dpar(ml)$mu[["x"]]))
  # penalized objective identity at the MAP solution
  nlp <- -stats::dnorm(fixef_by_dpar(map)$mu[["x"]], 0, 0.1, log = TRUE)
  raw_nll <- ml$obj$fn(map$opt$par)
  expect_lt(abs((-as.numeric(logLik_penalized(map))) - (raw_nll + nlp)) /
              abs(raw_nll), 1e-8)
  # and logLik() is the data part alone
  expect_lt(abs(-as.numeric(logLik(map)) - raw_nll) / abs(raw_nll), 1e-8)
  expect_output(print(map), "MAP")

  # near-flat priors reproduce ML
  map2 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
              prior = set_prior("normal(0, 1000)", class = "b"))
  expect_vector_equal(fixef_by_dpar(map2)$mu, fixef_by_dpar(ml)$mu, tol = 1e-3)
})

test_that("an sd prior regularizes a singular variance component", {
  set.seed(502)
  # 4 groups, tiny group effect: ML collapses the variance component
  dd <- data.frame(x = rnorm(80), g = factor(rep(1:4, 20)))
  dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(4, 0, 0.05)[dd$g], 1)
  ml <- suppressWarnings(frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd))
  sd_ml <- sqrt(varcorr_matrices(ml)[[1]][1, 1])

  map <- suppressWarnings(
    frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
        prior = set_prior("exponential(2)", class = "sd"))
  )
  sd_map <- sqrt(varcorr_matrices(map)[[1]][1, 1])
  # regularized away from the boundary, and finite standard errors
  expect_gt(sd_map, sd_ml - 1e-8)
  d <- diagnose(map, quiet = TRUE)
  expect_length(d$bad_se, 0)
})

test_that("hypothesis() on a MAP fit does not call it maximum likelihood", {
  # Through 0.80.0 print() said "a maximum-likelihood fit has no
  # posterior" on a fit with a prior, which has one. The columns stay NA
  # on every fit (?hypothesis, "The returned object" says why); the
  # printed reason is now true of both kinds of fit.
  set.seed(394)
  dd <- data.frame(x = rnorm(120), g = factor(rep(1:12, 10)))
  dd$y <- rnorm(120, 1 + 0.15 * dd$x + rnorm(12, 0, 0.6)[dd$g], 1)
  map <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
             prior = set_prior("normal(0, 1)", class = "b") +
               set_prior("cauchy(0, 2)", class = "sd", lb = 0))
  ml <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  for (f in list(map, ml)) {
    h <- hypothesis(f, c("x = 0", "x > 0"))
    hs <- hypothesis(f, "Intercept = 0", class = "sd", group = "g")
    expect_true(all(is.na(c(h$hypothesis$Evid.Ratio,
                            h$hypothesis$Post.Prob,
                            hs$hypothesis$Evid.Ratio))))
    out <- paste(utils::capture.output(print(hs)), collapse = "\n")
    expect_false(grepl("maximum-likelihood", out, fixed = TRUE))
    expect_true(grepl("brms computes them from posterior draws", out,
                      fixed = TRUE))
    expect_true(grepl("a fit has none", out, fixed = TRUE))
  }
})
