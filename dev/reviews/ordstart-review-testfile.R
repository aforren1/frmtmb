# Review of lane ordstart: run one test file of the worktree on one arm.
#   Rscript ordstart-review-testfile.R arm file
args <- commandArgs(trailingOnly = TRUE)
libs <- switch(args[1], base = "C:/Users/adf44/source/r/rellib-r20",
               lane = c("C:/Users/adf44/source/r/wt-ordstart-lib",
                        "C:/Users/adf44/source/r/rellib-r20"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
})
cat("lib:", find.package("frmtmb"), "\n")
r <- as.data.frame(testthat::test_file(args[2], package = "frmtmb",
                                       env = testthat::test_env("frmtmb"),
                                       reporter = "silent"))
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d\n",
            basename(args[2]), sum(r$passed), sum(r$failed),
            sum(r$error), sum(r$skipped), sum(r$warning)))
for (i in which(r$failed > 0 | r$error)) cat("  FAILED:", r$test[i], "\n")
