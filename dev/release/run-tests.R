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
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

# The StanHeaders 2.32.10 pin is gone (2026-09-17). Its two reasons were a
# tmbstan build that sampled a standard normal, fixed in tmbstan 1.2.1,
# and rstan failing to compile against StanHeaders 2.39.1, fixed by
# `CXX17FLAGS += -std=gnu++17` in the user Makevars
# (dev/tmbstan121-findings.md). R only reads <HOME>/.R/Makevars.win and
# HOME depends on the launcher, so the file is named explicitly. A
# populated Stan cache hides a compile failure everywhere except code
# that compiles something new, so these assertions are the stronger
# evidence and a green suite the weaker.
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
