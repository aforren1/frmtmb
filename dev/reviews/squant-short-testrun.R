# Short check of lane squant fix round 1, item 6: the lane's test file on
# the lane build ("lane") or on rellib-r22 alone ("base"). NOT_CRAN=true.
# Run (worktree root):
#   Rscript dev/reviews/squant-short-testrun.R lane > dev/reviews/squant-short-testrun-lane.txt 2>&1
#   Rscript dev/reviews/squant-short-testrun.R base > dev/reviews/squant-short-testrun-base.txt 2>&1
arm <- commandArgs(TRUE)[1]
Sys.setenv(NOT_CRAN = "true")
libs <- c("C:/Users/adf44/source/r/rellib-r22",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-squant-lib", libs)
.libPaths(libs)
p <- "frmtmb"
library(p, character.only = TRUE)
cat("arm:", arm, "| lib:", find.package(p), format(packageVersion(p)),
    "| quantile_elf exported:", exists("quantile_elf"), "\n")
res <- testthat::test_file(
  "tests/testthat/test-families-quantile.R", package = p,
  env = testthat::test_env(p), reporter = "silent")
r <- as.data.frame(res)
cat(sprintf("pass=%d fail=%d err=%d skip=%d warn=%d secs=%.1f tests=%d\n",
            sum(r$nb) - sum(r$failed) - sum(r$skipped), sum(r$failed),
            sum(r$error), sum(r$skipped), sum(r$warning), sum(r$real),
            nrow(r)))
bad <- r[r$failed > 0 | r$error, c("test", "failed", "error")]
if (nrow(bad)) {
  cat("failing or erroring tests:\n")
  for (i in seq_len(nrow(bad))) {
    cat(sprintf("  [fail %d err %s] %s\n", bad$failed[i], bad$error[i],
                bad$test[i]))
  }
  # the first failure message of each failing test
  for (t in res) {
    for (x in t$results) {
      if (inherits(x, c("expectation_failure", "expectation_error"))) {
        cat(sprintf("  -- %s: %s\n", t$test,
                    substr(gsub("\n", " ", conditionMessage(x)), 1, 200)))
        break
      }
    }
  }
}
cat("DONE\n")
