#' @srrstats {G5.7} Performance is tested across a range of data sizes,
#'   not at one size only. The scaling tests below build the same model
#'   at two sample sizes a hundredfold apart and require its cost to
#'   grow with the data and to stay within a linear envelope, which is
#'   the claimed complexity: the objective is vectorized over
#'   observations and the Laplace step over the sparse random-effect
#'   block. The cost is COUNTED, not timed: the operations on the tape
#'   and the bytes R allocates to build and evaluate it. Both are the
#'   same on a loaded machine and an idle one, where wall clock is not
#'   (dev/ciharden-findings.md).
#' @noRd
NULL

# Each count is bounded by TWICE the fold of n, not by the fold itself,
# because neither is exactly linear in n:
#
# - The tape's operations per observation depend on the response:
#   CppAD records nothing for a product by a constant 0 or 1, so a row
#   with y = 0 costs fewer operations, and the ratio follows the share
#   of zeros in the two samples. Over 31 GLM seed pairs it ran from
#   98.6 to 101.0 for a hundredfold n (dev/ciharden-rev-perfseeds.R).
# - Allocation: R's allocator and hash tables round some sizes up to a
#   power of two, and one-time allocations land in whichever call comes
#   first in the process, so the ratio depends on what ran before it
#   (58.8 in one order, 97.3 in another).
#
# A factor of 2 costs nothing against what these guard, which is
# superlinear: a quadratic is about fold^2. The allocation count is the
# canary for an elementwise sub-assignment loop in the objective (SPEC.md
# section 2a; dev/rtmb-pitfalls.md item 14), which copies the vector at
# every step; the tape count does NOT see that loop, since the loop
# records one operation per element either way (the mutation test below,
# and dev/ciharden-rev-perfmut.R). The tape count guards the tape's
# size. Neither count reads a floating-point result, so both are the
# same on every BLAS (dev/ciharden-perfcount2.R).
perf_census <- function(formula, d, random = NULL) {
  profmem <- isTRUE(capabilities("profmem"))
  tf <- tempfile()
  on.exit(unlink(tf))
  if (profmem) Rprofmem(tf, threshold = 0)
  fr <- frm(formula, data = d, dry_run = "frame")
  nll <- build_objective(fr)
  obj <- RTMB::MakeADFun(nll, fr$par_template, random = random,
                         silent = TRUE)
  obj$fn(obj$par)
  obj$gr(obj$par)
  if (profmem) Rprofmem(NULL)
  by <- if (profmem) {
    sum(suppressWarnings(as.numeric(sub(" *:.*", "", readLines(tf)))),
        na.rm = TRUE)
  } else {
    NA_real_
  }
  c(nodes = nrow(RTMB::GetTape(obj)$data.frame()), bytes = by)
}

expect_linear_cost <- function(small, large, fold) {
  expect_gt(large[["nodes"]], small[["nodes"]])
  expect_lte(large[["nodes"]] / small[["nodes"]], 2 * fold)
  skip_if_not(isTRUE(capabilities("profmem")),
              "R was built without memory profiling")
  expect_gt(large[["bytes"]], small[["bytes"]])
  expect_lte(large[["bytes"]] / small[["bytes"]], 2 * fold)
}

perf_glm_data <- function(n, seed) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n))
  d$y <- stats::rpois(n, exp(0.3 + 0.4 * d$x))
  d
}

test_that("a GLM's cost grows with n and stays linear in it", {
  # was a wall-clock bound, t_large < 100 * t_small, which failed under
  # load at 1.29 s and 1.31 s against 1.00
  skip_on_cran()
  f <- bf(y ~ x) + poisson()
  perf_census(f, perf_glm_data(1000L, 70))   # the first call loads caches
  small <- perf_census(f, perf_glm_data(1000L, 71))
  large <- perf_census(f, perf_glm_data(100000L, 72))
  expect_linear_cost(small, large, 100)
})

test_that("a GLMM's taping stays linear in n", {
  # The canary against an observation-length loop in the objective. It
  # was a 20 s wall-clock bound on taping n = 100000; the count catches
  # the same regression by the same orders of magnitude with no clock.
  skip_on_cran()
  make <- function(n, seed) {
    set.seed(seed)
    n_g <- 500L
    d <- data.frame(x = stats::rnorm(n),
                    g = factor(sample.int(n_g, n, replace = TRUE),
                               levels = seq_len(n_g)))
    d$y <- stats::rpois(n, exp(0.2 + 0.3 * d$x +
                                 stats::rnorm(n_g, 0, 0.4)[d$g]))
    d
  }
  f <- bf(y ~ x + (1 | g)) + poisson()
  perf_census(f, make(1000L, 6), random = "b")
  small <- perf_census(f, make(1000L, 7), random = "b")
  large <- perf_census(f, make(100000L, 8), random = "b")
  expect_linear_cost(small, large, 100)
})

test_that("the census catches an elementwise loop in the objective", {
  # The guard on the guard: item 14's loop put INSIDE the objective that
  # build_objective() returns must fail the allocation bound. A tenfold
  # n keeps it cheap; the loop's ratio there is about 97 against a
  # bound of 20 (dev/ciharden-rev-perfmut.R; 9676 for a hundredfold n).
  # The tape count passes it, which is why the allocation count is the
  # canary.
  skip_on_cran()
  skip_if_not(isTRUE(capabilities("profmem")),
              "R was built without memory profiling")
  orig <- build_objective
  n_now <- 0L
  local_mocked_bindings(build_objective = function(fr, ...) {
    nll <- orig(fr, ...)
    n <- n_now
    function(pars) {
      val <- nll(pars)
      "[<-" <- RTMB::ADoverload("[<-")
      b <- pars[[1L]][1L]
      m <- b * rep(1, n)
      for (i in seq_len(n)) m[i] <- b * 2
      val + 1e-12 * sum(m)
    }
  })
  f <- bf(y ~ x) + poisson()
  n_now <- 1000L
  perf_census(f, perf_glm_data(1000L, 70))
  small <- perf_census(f, perf_glm_data(1000L, 71))
  n_now <- 10000L
  large <- perf_census(f, perf_glm_data(10000L, 72))
  expect_gt(large[["bytes"]] / small[["bytes"]], 2 * 10)
  expect_lte(large[["nodes"]] / small[["nodes"]], 2 * 10)
})
