# Grammar fuzz tier. Off by default: it costs minutes and its job is to
# find new failures, not to guard old ones, so it must never be part of
# the signal a normal run gives. See helper-fuzz.R for the design.
#
#   Sys.setenv(FRMTMB_FUZZ = "true")
#   testthat::test_file("tests/testthat/test-fuzz.R")
#   FRMTMB_FUZZ_N = 20   subsample the plan for a quick smoke run
#
# The tier is green: every finding the full plan produces is on one of
# the three lists in helper-fuzz.R, each carrying the reason it stands -
# FUZZ_KNOWN_PENDING (fixed on an unmerged sibling branch),
# FUZZ_KNOWN_REFUSAL (refused on purpose), FUZZ_KNOWN_DIVERGENCE (a
# deliberate departure from brms). dev/fuzz-findings.md is the long
# form. A red line here is therefore a new defect, and the failure
# message prints it.

#' @srrstats {G5.6b} Parameter recovery is checked across many random
#'   seeds. G5.6b allows this to live in an extended rather than a
#'   regular test suite when it is long-running, which is the case here.
#'   `fuzz_plan()` derives a distinct seed per generated model
#'   (`seed + i * 977L`) from one plan seed, so a run is a sweep of
#'   several hundred seeds, and the `confint_coverage` invariant pools
#'   them: the 95% Wald intervals across the whole plan must cover the
#'   known simulation truth `FUZZ_BETA_X = 0.4` at least
#'   `qbinom(1e-4, n, 0.95)` times. That is a one-in-ten-thousand
#'   binomial tail, so only a broken interval, not an unlucky plan, can
#'   trip it.
#' @srrstats {G5.9b} The same mechanism is the noise-susceptibility test
#'   for seed variation: running under hundreds of different seeds and
#'   initial conditions must not change the conclusion, and the
#'   invariants (`predict_eq_fitted`, `loglik_identity`,
#'   `row_permutation`, `unit_weights`, `simulate_mean`, `vcov_dim`,
#'   `vcov_psd`, `summary_prints`, `confint_wald`) are
#'   required to hold on every one of them. `row_permutation` is a
#'   metamorphic version of the same idea: reordering the rows must not
#'   change the likelihood. Stability under different initial conditions
#'   is also checked outside the fuzzer, by `frm_allfit()` agreeing
#'   across optimizers and by the profiled fit reproducing the plain one.
#' @noRd
NULL

test_that("pairwise grammar fuzz finds no new invariant violations", {
  skip_if_not(identical(Sys.getenv("FRMTMB_FUZZ"), "true"),
              "set FRMTMB_FUZZ=true to run the grammar fuzz tier")

  size <- suppressWarnings(as.integer(Sys.getenv("FRMTMB_FUZZ_N")))
  # fuzz_plan() caps the plan itself, so a smoke run keeps the prefix of
  # the greedy cover (the rows that buy the most pairs) plus the refusal
  # probes, rather than a random slice of the full plan
  plan <- fuzz_plan(seed = 20260901L,
                    size = if (is.na(size)) 300L else max(size, 1L))
  res <- fuzz_run(plan)
  sm <- fuzz_summary(res)

  new <- Filter(function(f) identical(f$class, "real_new"), sm$triaged)
  gen <- Filter(function(f) identical(f$class, "generator"), sm$triaged)
  expect_equal(length(gen), 0L,
               info = paste(vapply(gen, fuzz_format, ""), collapse = "\n"))
  expect_equal(length(new), 0L,
               info = paste(vapply(new, fuzz_format, ""), collapse = "\n"))
})
