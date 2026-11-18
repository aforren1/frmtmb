# Review of lane squant: rerun the lane's new test file on the lane build.
# Run (worktree root), NOT_CRAN=true:
#   Rscript dev/reviews/squant-review-testrun.R > dev/reviews/squant-review-testrun.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
p <- "frmtmb"
library(p, character.only = TRUE)
cat("lib:", find.package(p), "\n")
r <- as.data.frame(testthat::test_file(
  "tests/testthat/test-families-quantile.R", package = p,
  env = testthat::test_env(p), reporter = "silent"))
cat(sprintf("pass=%d fail=%d err=%d skip=%d warn=%d secs=%.1f\n",
            sum(r$nb) - sum(r$failed) - sum(r$skipped), sum(r$failed),
            sum(r$error), sum(r$skipped), sum(r$warning), sum(r$real)))
