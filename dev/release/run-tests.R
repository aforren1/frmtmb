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

# Everything runs inside local(), so the runner leaves no variable in
# the global environment: a test whose model reads a variable the data
# lacks falls back, by R's rule, to the formula environment and then to
# the global one, where this runner's `f` (the test file's path) once
# answered for a missing `gp(x, by = f)` column (0.68.0 consolidation,
# test-gp-by.R). R CMD check's test_check() has no such globals.
local({
# FRMTMB_RELLIB puts a library ahead of the release one, for a "before"
# arm that holds core alone (dev/rel068-mutants.R); the release library
# stays behind it for the extensions
LIB <- "C:/Users/adf44/source/r/rellib-r7"
.libPaths(c(if (nzchar(Sys.getenv("FRMTMB_RELLIB"))) {
  Sys.getenv("FRMTMB_RELLIB")
}, LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

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
# which build ran, so a log can be checked against the library it claims
cat("lib: ", p, " ", as.character(utils::packageVersion(p)), " ",
    dirname(find.package(p)), "; frmtmb ",
    as.character(utils::packageVersion("frmtmb")), " ",
    dirname(find.package("frmtmb")), "\n", sep = "")

# Tests must run with the package NAMESPACE as parent, the way
# test_check() does under R CMD check. From the global environment an
# importFrom()'d symbol is invisible, so a file using one reports
# "could not find function" and reads as a regression. frmtmb.learn
# imports five functions from frmtmb.eam, which is how this was found.
res <- NULL
r <- tryCatch({
  res <- test_file(f, package = p, env = testthat::test_env(p),
                   reporter = "silent")
  as.data.frame(res)
}, error = function(e) {
  cat("RESULT ", basename(f), " LOADERROR ", conditionMessage(e), "\n",
      sep = "")
  NULL
})

# every failure, error, skip and escaped warning with its test, so a
# count in a summary can be read back to what it counts
if (!is.null(res)) {
  for (t in res) {
    for (x in t$results) {
      k <- if (inherits(x, "expectation_failure")) "FAIL" else
        if (inherits(x, "expectation_error")) "ERROR" else
          if (inherits(x, "expectation_skip")) "SKIP" else
            if (inherits(x, "expectation_warning")) "WARN" else NA
      if (!is.na(k)) {
        cat("DETAIL ", k, " [", t$test, "] ",
            gsub("\n", " | ", conditionMessage(x)), "\n", sep = "")
      }
    }
  }
}

if (!is.null(r)) {
  # warn= last, so a reader that parses the four counts before it is
  # unchanged; an escaped warning is a defect under dev/lane-rules.md
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      " warn=", sum(r$warning), "\n", sep = "")
}
})
