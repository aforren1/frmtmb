# Run one test file against a library, print a RESULT line and every
# failure. Usage:
#   Rscript dev/quadseam2-test.R <libpath|base> <package> <test file path>
a <- commandArgs(TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
pkg <- a[2]
f <- a[3]
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
cat("lib:", dirname(find.package("frmtmb")), "\n")
cat("pkglib:", dirname(find.package(pkg)), "\n")
suppressPackageStartupMessages(library(testthat))
res <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                           reporter = "silent")
df <- as.data.frame(res)
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d\n",
            basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
for (r in res) {
  for (x in r$results) {
    if (inherits(x, c("expectation_failure", "expectation_error",
                      "expectation_warning"))) {
      cat("FAIL [", r$test, "] ", gsub("\\s+", " ",
                                       conditionMessage(x)), "\n", sep = "")
    }
  }
}
