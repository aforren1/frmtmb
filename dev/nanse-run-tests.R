# One test file, one R process, against the lane library, recording
# every firing of the new standard-error warning (se_lost_message() is
# built exactly once per firing) even where a test muffles it.
#   Rscript dev/nanse-run-tests.R <package> <test file> <fire log>
LIB <- "C:/Users/adf44/source/r/wt-nanse-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
suppressMessages(library(testthat))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
fire <- a[3]
suppressMessages(library(p, character.only = TRUE))
cat("lib: ", find.package("frmtmb"), " | ", find.package(p), "\n", sep = "")
ns <- asNamespace("frmtmb")
if (exists("se_lost_message", ns)) {
  assign(".nanse_fire", fire, envir = globalenv())
  suppressMessages(trace(
    "se_lost_message", where = ns, print = FALSE,
    exit = quote(cat(gsub("[\r\n\t]+", " ", returnValue()), "\n",
                     file = get(".nanse_fire", envir = globalenv()),
                     append = TRUE))))
}
r <- tryCatch(
  as.data.frame(test_file(f, package = p, env = testthat::test_env(p),
                          reporter = "silent")),
  error = function(e) {
    cat("RESULT ", basename(f), " LOADERROR ", conditionMessage(e), "\n",
        sep = "")
    NULL
  })
if (!is.null(r)) {
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      " warn=", sum(r$warning), "\n", sep = "")
  for (i in seq_len(nrow(r))) {
    if (r$failed[i] > 0 || isTRUE(r$error[i]) || r$warning[i] > 0) {
      cat(sprintf("  TEST %s | fail=%d err=%s warn=%d\n", r$test[i],
                  r$failed[i], r$error[i], r$warning[i]))
    }
  }
}
