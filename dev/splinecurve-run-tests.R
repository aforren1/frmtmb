# Lane splinecurve: one test file, one R process, against a named library.
# A copy of dev/release/run-tests.R with the library as an argument, so
# the same file can run against the lane's private library and against
# the read-only release library for the seen-to-fail record.
#
# Usage: Rscript splinecurve-run-tests.R <lib> <package> <test-file>

a <- commandArgs(trailingOnly = TRUE)
LIB <- a[1]
p <- a[2]
f <- a[3]
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
suppressMessages(library(testthat))
suppressMessages(library(p, character.only = TRUE))
cat("LIB ", LIB, " frmtmb ", format(packageVersion("frmtmb")), " from ",
    dirname(find.package("frmtmb")), "; ", p, " ",
    format(packageVersion(p)), " from ", dirname(find.package(p)), "\n",
    sep = "")

res <- tryCatch(
  test_file(f, package = p, env = testthat::test_env(p),
            reporter = "silent"),
  error = function(e) {
    cat("RESULT ", basename(f), " LOADERROR ", conditionMessage(e), "\n",
        sep = "")
    NULL
  })

if (!is.null(res)) {
  # the failing expectations themselves, since the silent reporter keeps
  # them and a count alone says nothing about which assertion broke
  for (t in res) {
    for (ex in t$results) {
      if (inherits(ex, c("expectation_failure", "expectation_error"))) {
        cat("FAILED [", t$test, "] ", conditionMessage(ex), "\n", sep = "")
      }
    }
  }
  r <- as.data.frame(res)
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      "\n", sep = "")
}
