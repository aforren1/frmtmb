# Canary against accidental elementwise operations during taping (SPEC.md
# section 2a): taping a large GLMM must stay fast. A vectorized objective
# tapes this model in well under a second; an observation-length loop or
# elementwise sub-assignment would blow past the bound by orders of
# magnitude.
test_that("tape construction stays bounded for a large GLMM", {
  skip_on_cran()
  set.seed(7)
  n_g <- 500L
  n <- 100000L
  dd <- data.frame(
    x = rnorm(n),
    g = factor(sample.int(n_g, n, replace = TRUE))
  )
  dd$y <- rpois(n, exp(0.2 + 0.3 * dd$x + rnorm(n_g, 0, 0.4)[dd$g]))

  elapsed <- system.time({
    fr <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd,
                 dry_run = "frame")
    nll <- frmtmb:::build_objective(fr)
    obj <- RTMB::MakeADFun(nll, fr$par_template, random = "b",
                           silent = TRUE)
    obj$fn(obj$par)
  })[["elapsed"]]

  expect_lt(elapsed, 20)
})
