# A quadrature fit escapes a family's stationary point on the Laplace
# objective (quad_fit()) and again on the rule (fit_assembled()); the
# second escape's record replaced the first's, so a fit the Laplace escape
# moved and the rule's left alone reported a gain of 0
# (stationary_escape_keep(), release review of 0.82.0, minor 6).

test_that("the escape that moved a quadrature fit keeps its record", {
  lap <- c(starts = 2, gain = 24.27)
  # the rule's escape ran and gained nothing: the Laplace record stays
  opt <- list(objective = 1, stationary_escape = c(starts = 2, gain = 0))
  got <- stationary_escape_keep(opt, lap)
  expect_identical(got$stationary_escape, lap)
  expect_identical(got$stationary_escape_laplace, lap)
  # the rule's escape moved the fit: its record stays, the Laplace one
  # beside it
  rule <- c(starts = 2, gain = 0.5)
  opt$stationary_escape <- rule
  got <- stationary_escape_keep(opt, lap)
  expect_identical(got$stationary_escape, rule)
  expect_identical(got$stationary_escape_laplace, lap)
})

test_that("the rule's escape runs after the Laplace escape ran", {
  # The rule is optimized from the Laplace escape's optimum and can walk
  # back onto alpha = 0. A build that skipped the rule's escape whenever
  # the Laplace one had run stopped this fit (alpha ~ x, both coefficients
  # started at 0) with both at 1e-15, 6.67 below 0.82.0, code 0 and no
  # warning (review of lane quadseam2, F1;
  # dev/reviews/quadseam2-review-esc2.R, ALPHAX seed 11; the data of the
  # two designs before it are drawn first, as there)
  skip_on_cran()
  sn <- function(n, alpha) {
    dl <- alpha / sqrt(1 + alpha^2)
    dl * (abs(stats::rnorm(n)) - sqrt(2 / pi)) +
      sqrt(1 - dl^2) * stats::rnorm(n)
  }
  set.seed(1011)
  stats::rnorm(60 + 15 + 60)
  stats::rnorm(50 + 10)
  sn(50, 1.5)
  d <- data.frame(x = stats::rnorm(80), g = factor(rep(1:16, each = 5)))
  d$y <- 0.5 * d$x + stats::rnorm(16, 0, 1)[d$g] + sn(80, 2)
  fo <- bf(y ~ x + (1 | g), sigma ~ 1, alpha ~ x)
  # alpha runs toward its end (the half-normal limit), so the stop carries
  # the optimizer's false-convergence code; under OpenBLAS 0.3.26 and
  # 0.3.32 with LAPACK it also carries this standard-error warning, which
  # the reference BLAS does not give (dev/quadseam2-findings.md, short
  # check). Any other warning reaches the test
  fit <- allow_warnings(
    frm(fo, family = skew_normal(), data = d, quadrature = TRUE,
        start = list(betad = c(0, 0, 0))),
    allowed = c("Optimizer did not report convergence: false convergence",
                paste("Standard errors are not available for 4 of 6",
                      "parameters. (Intercept), x, alpha_(Intercept),",
                      "alpha_x: the likelihood is flat")))
  esc <- fit$opt[["stationary_escape"]]
  lap <- fit$opt[["stationary_escape_laplace"]]
  expect_false(is.null(lap))
  # the rule's own escape gained beyond the Laplace one's, far above the
  # optimizer's tolerance, and its record is the one kept
  expect_false(identical(esc, lap))
  expect_gt(esc[["gain"]], abs(as.numeric(logLik(fit))) * 1e-3)
  a <- fit$estimates$betad[grep("^alpha", names(fit$estimates$betad))]
  expect_true(all(abs(a) > 1))
})
