# Run one test file and print a summary line that cannot be mistaken
# for a clean one when the file aborted: the file count, the skip count
# and the test count are all on it, and `error` is summed as well as
# `failed`, which a previous runner in this repository did not do.
a <- commandArgs(TRUE)
Sys.setenv(NOT_CRAN = "true")
suppressMessages({
  library(testthat)
  library(frmtmb)
  library(frmtmb.ode)
})
setwd(a[1])
r <- tryCatch(test_file(a[2], reporter = "summary",
                        package = "frmtmb.ode"),
              error = function(e) e)
if (inherits(r, "condition")) {
  cat("\nFILE ", a[2], " ABORTED: ", conditionMessage(r), "\n", sep = "")
} else {
  d <- as.data.frame(r)
  cat("\nFILE ", a[2], " pass=", sum(d$passed), " fail=", sum(d$failed),
      " err=", sum(d$error), " skip=", sum(d$skipped), " tests=",
      nrow(d), "\n", sep = "")
  bad <- d[d$failed > 0 | d$error > 0, "test"]
  if (length(bad)) cat("BROKE: ", paste(bad, collapse = " | "), "\n",
                       sep = "")
}
