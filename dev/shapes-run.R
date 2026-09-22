# One test file, one R process, against the lane's private library.
# Modeled on dev/release/run-tests.R; the library is the lane's, and the
# Makevars assertion is kept because a Stan program only compiles with
# the user flag.
#   Rscript dev/shapes-run.R <package> <path-to-test-file>

LIB <- "C:/Users/adf44/source/r/shapes-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}

suppressMessages(library(testthat))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
suppressMessages(library(p, character.only = TRUE))

r <- tryCatch(
  as.data.frame(test_file(f, package = p, env = testthat::test_env(p),
                          reporter = "silent")),
  error = function(e) {
    cat("RESULT ", basename(f), " LOADERROR ", conditionMessage(e), "\n",
        sep = "")
    NULL
  })

if (!is.null(r)) {
  for (i in seq_len(nrow(r))) {
    if (r$failed[i] > 0 || isTRUE(r$error[i])) {
      cat("  [", if (isTRUE(r$error[i])) "ERROR" else "FAIL", "] ",
          r$test[i], "\n", sep = "")
      for (res in r$result[[i]]) {
        if (inherits(res, c("expectation_failure", "expectation_error"))) {
          cat("      ", gsub("\n", "\n      ", conditionMessage(res)),
              "\n", sep = "")
        }
      }
    }
  }
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      "\n", sep = "")
}
