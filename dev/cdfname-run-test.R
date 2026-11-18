# One test file in one R process.
#   Rscript dev/cdfname-run-test.R <pkg> <test file path> [base]
# "base" drops the lane library, for the before arm on rellib-r22.
args <- commandArgs(trailingOnly = TRUE)
pkg <- args[[1L]]
f <- args[[2L]]
base <- length(args) > 2L && args[[3L]] == "base"
lib <- "C:/Users/adf44/source/r/wt-cdfname-lib"
.libPaths(c(if (!base) lib, "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "|", pkg, find.package(pkg), "\n")
library(testthat)
res <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                           reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
fails <- df[df$failed > 0 | df$error, , drop = FALSE]
for (i in seq_len(nrow(fails))) cat("FAILED:", fails$test[i], "\n")
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), sum(df$passed), sum(df$failed), sum(df$error),
            sum(df$skipped), sum(df$warning)))
if (nrow(fails)) {
  for (r in res) for (x in r$results) {
    if (inherits(x, c("expectation_failure", "expectation_error"))) {
      cat("--", r$test, ":", conditionMessage(x), "\n")
    }
  }
}
