# Run one test file in its own process against whichever library
# SKEWINIT_LIB names, and print counts that include errors.
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
f <- commandArgs(TRUE)[[1L]]
library(testthat)
library(frmtmb)
cat("frmtmb from", dirname(dirname(getNamespaceInfo("frmtmb", "path"))), "\n")
setwd("C:/Users/adf44/source/r/frmtmb-wt-skewinit/tests/testthat")
res <- test_file(f, reporter = "silent")
d <- as.data.frame(res)
cat(sprintf("FILE %s PASS=%d FAIL=%d ERROR=%d WARN=%d SKIP=%d\n", f,
            sum(d$passed), sum(d$failed), sum(d$error), sum(d$warning),
            sum(d$skipped)))
for (i in seq_len(nrow(d))) {
  if (d$failed[i] > 0 || d$error[i] > 0) {
    cat("--", d$test[i], "\n")
    for (r in res[[i]]$results) {
      if (inherits(r, "expectation_failure") ||
            inherits(r, "expectation_error")) {
        cat("   ", gsub("\n", " ", conditionMessage(r)), "\n")
      }
    }
  }
}
