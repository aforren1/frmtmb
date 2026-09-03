#' @srrstats {G5.7} Performance is tested across a range of data sizes,
#'   not at one size only. The scaling test below fits the same model at
#'   two sample sizes a hundredfold apart and requires the run time to
#'   grow with the data and to stay within a linear envelope, which is
#'   the claimed complexity: the objective is vectorized over
#'   observations and the Laplace step over the sparse random-effect
#'   block. The assertion is deliberately coarse, because wall clock on a
#'   shared machine is not reproducible; a quadratic regression in the
#'   objective would still break it by orders of magnitude. The tape
#'   canary in this file pins the constant factor at one large size.
#' @noRd
NULL

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

test_that("fit time grows with n and stays within a linear envelope", {
  skip_on_cran()
  make <- function(n, seed) {
    set.seed(seed)
    d <- data.frame(x = stats::rnorm(n))
    d$y <- stats::rpois(n, exp(0.3 + 0.4 * d$x))
    d
  }
  fit_time <- function(d) {
    # these fits exist to be timed; a hairline gradient-threshold
    # warning on one BLAS build is irrelevant to the scaling assertion
    f <- function() suppressWarnings(frm(bf(y ~ x) + poisson(), data = d))
    f()                                    # warm the tape machinery
    stats::median(replicate(3, system.time(f())[["elapsed"]]))
  }
  small <- make(1000L, 71)
  large <- make(100000L, 72)               # a hundredfold more data
  t_small <- fit_time(small)
  t_large <- fit_time(large)

  expect_gt(t_large, t_small)
  # a hundredfold more data must not cost more than a hundredfold more
  # time; the floor keeps a sub-millisecond t_small from making the
  # bound meaningless
  expect_lt(t_large, 100 * max(t_small, 0.01))
})
