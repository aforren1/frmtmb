# One test file, one process, on the lane's core: for the suites the lane did not
# rerun. Usage: Rscript dev/review-impsamp-fr1-suite.R <pkg> <file>
a <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-impsamp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(a[[1]], character.only = TRUE)
  library(testthat)
})
f <- file.path("C:/Users/adf44/source/r/frmtmb-wt-impsamp/extensions", a[[1]], "tests/testthat", a[[2]])
r <- as.data.frame(test_file(f, package = a[[1]], env = test_env(a[[1]]), reporter = "silent",
                             stop_on_failure = FALSE))
cat(sprintf("RESULT %s %s core=%s pass=%d fail=%d error=%d skip=%d warn=%d\n", a[[1]], a[[2]],
            dirname(find.package("frmtmb")), sum(r$passed), sum(r$failed), sum(r$error),
            sum(r$skipped), sum(r$warning)))
