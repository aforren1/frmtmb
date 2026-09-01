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
