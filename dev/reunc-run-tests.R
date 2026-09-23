# Lane wt-reunc: dev/release/run-tests.R with the lane's library order.
# The release script is not edited; this is a copy whose ONLY change is
# the library path, which it takes from REUNC_LIB (default: the lane's
# private library, then the round's read-only base build rellib-r3).
# REUNC_REPORTER=check prints each failure in full, which is what a
# seen-failing log needs; the default is the release script's silent
# reporter and its one RESULT line.
#
# Usage: Rscript dev/reunc-run-tests.R <package> <path-to-test-file>

LIB <- Sys.getenv("REUNC_LIB", "C:/Users/adf44/source/r/reunc-lib")
.libPaths(unique(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))

mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
stopifnot(utils::packageVersion("tmbstan") >= "1.2.1",
          any(grepl("-std=gnu++17",
                    readLines(tools::makevars_user(), warn = FALSE),
                    fixed = TRUE)))

suppressMessages(library(testthat))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
suppressMessages(library(p, character.only = TRUE))
cat("LIB ", find.package(p), " frmtmb ", find.package("frmtmb"), "\n",
    sep = "")

rep <- Sys.getenv("REUNC_REPORTER", "silent")
r <- tryCatch(
  as.data.frame(test_file(f, package = p, env = testthat::test_env(p),
                          reporter = rep, stop_on_failure = FALSE)),
  error = function(e) {
    cat("RESULT ", basename(f), " LOADERROR ", conditionMessage(e), "\n",
        sep = "")
    NULL
  })

if (!is.null(r)) {
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      "\n", sep = "")
}
