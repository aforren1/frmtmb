# Lane silent2: one test file in one process.
#   Rscript dev/silent2-runtest.R base|lane <package> <test file path>
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
pkg <- args[2]
f <- args[3]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-silent2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(testthat)
  library(pkg, character.only = TRUE)
})
cat("lib:", pkg, find.package(pkg), "frmtmb", find.package("frmtmb"), "\n")
r <- as.data.frame(testthat::test_file(f, package = pkg,
                                       env = testthat::test_env(pkg),
                                       reporter = "silent"))
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d\n",
            basename(f), sum(r$passed), sum(r$failed), sum(r$error),
            sum(r$skipped),
            sum(r$warning)))
bad <- r[r$failed > 0 | r$error | r$warning > 0,
         c("test", "failed", "error", "warning")]
if (nrow(bad)) print(bad)
if (sum(r$warning) > 0) {
  # the escaped warnings' texts, from a run of the same file that keeps
  # them
  res <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                             reporter = "silent")
  for (t in res) {
    for (e in t$results) {
      if (inherits(e, "expectation_warning")) {
        cat("WARNING in", t$test, ":", conditionMessage(e), "\n")
      }
    }
  }
}
