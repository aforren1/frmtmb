# Reviewer of lane defects: run one generated brms-suite file from a
# scratch copy of the test directory that carries the probe helper.
#   Rscript dev/defects-rev-probe-run.R <pkg> <dir> <file> <log>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-defects-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
Sys.setenv(R_MAKEVARS_USER = mk, NOT_CRAN = "true",
           FRMTMB_BRMS_FIT_TESTS = "true", FRMTMB_BRMSPORT_PKG = a[1],
           DEFREV_PROBE_LOG = normalizePath(a[4], winslash = "/",
                                            mustWork = FALSE))
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
suppressMessages(library(a[1], character.only = TRUE))
r <- as.data.frame(test_file(file.path(a[2], a[3]), package = a[1],
                             env = testthat::test_env(a[1]),
                             reporter = "silent", stop_on_failure = FALSE))
cat("RESULT", a[1], a[3], "pass=", sum(r$passed), "fail=", sum(r$failed),
    "err=", sum(r$error), "skip=", sum(r$skipped), "warn=",
    sum(r$warning), "\n")
