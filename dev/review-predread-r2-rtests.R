# Short check of fix round 2: the tests that go with release fixes F1 to F4. Each block is
# self-contained so it can move into test-me-newdata.R as it is. Run on the lane's build,
# where each must fail, and on the build with the fixes, where each must pass.
#   Rscript dev/review-predread-r2-rtests.R <library>
a <- commandArgs(TRUE)
.libPaths(c(a[1L], "C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
cat("lib:", dirname(find.package("frmtmb")), "\n")
r <- as.data.frame(test_file("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-r2-rtests-file.R",
                             reporter = "silent", stop_on_failure = FALSE, env = test_env("frmtmb")))
print(r[, c("test", "failed", "error", "passed")])
cat("FAILED tests:", sum(r$failed > 0 | r$error), "of", length(unique(r$test)), "\n")
