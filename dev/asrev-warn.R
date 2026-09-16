## Identify the one WARN in the in-check suite. R CMD check leaves
## NOT_CRAN unset, so this run must too; the check reporter only counts
## warnings, so this one uses a reporter that names them.
.libPaths(c("C:/Users/adf44/source/r/asrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.unsetenv("NOT_CRAN")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
setwd("tests/testthat")
res <- testthat::test_dir(".", package = "frmtmb", reporter = "silent",
                          stop_on_failure = FALSE)
d <- as.data.frame(res)
cat(sprintf("\nTOTAL pass %d fail %d error %d skip %d warn %d\n",
            sum(d$passed), sum(d$failed), sum(d$error), sum(d$skipped),
            sum(d$warning)))
w <- d[d$warning > 0L, c("file", "test", "warning")]
cat("\nfiles with a warning:\n")
print(w, row.names = FALSE)
for (r in res) {
  for (x in r$results) {
    if (inherits(x, "expectation_warning")) {
      cat("\n---- WARNING ----\n")
      cat("file :", r$file, "\n")
      cat("test :", r$test, "\n")
      cat("msg  :", conditionMessage(x), "\n")
      cat("srcref:", paste(utils::capture.output(print(x$srcref)),
                           collapse = " "), "\n")
    }
  }
}
