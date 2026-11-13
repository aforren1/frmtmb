# Review of lane arspat: how often frmtmb_autocor_gap fires across the core
# test files that fit ar(cov = TRUE). One file per call; counts only.
# Usage: Rscript arspat-review-noise.R <file> ...
source("dev/reviews/arspat-review-common.R")
suppressPackageStartupMessages(library(testthat))
Sys.setenv(NOT_CRAN = "true")
for (f in commandArgs(trailingOnly = TRUE)) {
  cnt <- 0L
  where <- character(0)
  withCallingHandlers(
    suppressWarnings(testthat::test_file(
      file.path("tests/testthat", f), package = "frmtmb",
      env = testthat::test_env("frmtmb"), reporter = "silent")),
    frmtmb_autocor_gap = function(m) {
      cnt <<- cnt + 1L
      where <<- c(where, substr(conditionMessage(m), 1, 50))
    })
  cat(sprintf("%-36s gap messages: %d %s\n", f, cnt,
              paste(unique(where), collapse = " | ")))
}
