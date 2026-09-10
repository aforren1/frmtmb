# One test file, one R process. Invoked by the drivers beside it.
#
# WHY THIS EXISTS AS A FILE IN THE REPOSITORY. The release harness used
# to live in the session's temp directory, which is cleaned between
# sessions. It evaporated mid-consolidation at 0.55.3 and
# `powershell -File <missing>` exits 0, so the driver reported success
# having run nothing. A harness is a guard, and every guard this
# project has built failed open on its first try; this one now lives
# where it can be read, reviewed and diffed.
#
# Usage: Rscript run-tests.R <package> <path-to-test-file>

LIB <- "C:/Users/adf44/source/r/rellib-r3"
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

# The pin has to win over the user library's StanHeaders 2.39.1, which
# rstan 2.32.7 cannot compile against. A populated Stan cache hides that
# everywhere except code that compiles something new, so a green suite
# is the weaker evidence and this assertion is the stronger one.
stopifnot(identical(as.character(utils::packageVersion("StanHeaders")),
                    "2.32.10"))

suppressMessages(library(testthat))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
suppressMessages(library(p, character.only = TRUE))

# Tests must run with the package NAMESPACE as parent, the way
# test_check() does under R CMD check. From the global environment an
# importFrom()'d symbol is invisible, so a file using one reports
# "could not find function" and reads as a regression. frmtmb.learn
# imports five functions from frmtmb.eam, which is how this was found.
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
      "\n", sep = "")
}
