# Item 3: the OLD test files (read-only copy in frmtmb-wt-r16base) against the
# lane's build. Which old tests fail, and on which expectation?
# Usage: Rscript dev/review-impsamp-oldtests.R <file.R> [pkgdir]
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
args <- commandArgs(TRUE)
root <- "C:/Users/adf44/source/r/frmtmb-wt-r16base/tests/testthat"
pkg <- "frmtmb"
if (length(args) >= 2) {
  root <- file.path("C:/Users/adf44/source/r/frmtmb-wt-r16base/extensions",
                    args[[2]], "tests/testthat")
  pkg <- args[[2]]
  library(pkg, character.only = TRUE)
}
cat("lib:", find.package("frmtmb"), "\n")
res <- test_file(file.path(root, args[[1]]), package = pkg,
                 env = test_env(pkg), reporter = "silent",
                 stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n", args[[1]],
            sum(df$passed), sum(df$failed), sum(df$error), sum(df$skipped),
            sum(df$warning)))
for (r in res) {
  bad <- Filter(function(x) inherits(x, c("expectation_failure", "expectation_error")),
                r$results)
  if (length(bad)) {
    cat("\nFAILED TEST:", r$test, "\n")
    for (b in bad) cat("  -", substr(conditionMessage(b), 1, 300), "\n")
  }
}
