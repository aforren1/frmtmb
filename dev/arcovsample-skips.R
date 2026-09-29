# Lane wt-arcovsample: list the SKIP reasons of one test file, because a
# skip count is not evidence until its reasons are read (lane-rules.md).
#   Rscript dev/arcovsample-skips.R <package> <file>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-arcovsample-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(testthat))
suppressMessages(library(a[1L], character.only = TRUE))
r <- test_file(a[2L], package = a[1L],
               env = testthat::test_env(a[1L]), reporter = "silent")
for (tr in r) {
  for (x in tr$results) {
    if (inherits(x, "expectation_skip")) {
      cat("SKIP [", tr$test, "] ", conditionMessage(x), "\n", sep = "")
    }
  }
}
