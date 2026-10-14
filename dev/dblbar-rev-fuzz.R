# Reviewer, lane dblbar: the grammar fuzz tier (FRMTMB_FUZZ=true), whose
# helper the lane edited (dbar/dbar_f block counts, one known-pending
# entry removed) without running it. Usage: Rscript ... base|lane
# Output: stdout; the triaged findings are printed by the test.
arm <- commandArgs(TRUE)[1]
lib <- if (identical(arm, "base")) character(0) else
  "C:/Users/adf44/source/r/wt-dblbar-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_FUZZ = "true")
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
cat("lib:", find.package("frmtmb"), "\n")
f <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar/tests/testthat/test-fuzz.R"
r <- as.data.frame(testthat::test_file(f, package = "frmtmb",
                                       env = testthat::test_env("frmtmb"),
                                       reporter = "location"))
cat(sprintf("RESULT %s fuzz pass=%d fail=%d err=%d skip=%d\n", arm,
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped)))
