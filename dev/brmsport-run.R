# One generated brms-suite test file, one R process, against this lane's
# private library.
#
#   Rscript dev/brmsport-run.R <package> <test file> [record.tsv]
#
# With a third argument the harness RECORDS every assertion's outcome to
# that file and asserts nothing (dev/brmsport-ledger.R reads it); without
# one the file runs as the gated tier runs it, asserting the verdicts.
lib <- "C:/Users/adf44/source/r/brmsport-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
stopifnot(identical(as.character(utils::packageVersion("StanHeaders")),
                    "2.32.10"))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_BRMSPORT_PKG = p)
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE =
               "C:/Users/adf44/source/r/frmtmb-wt-brmsport/dev/stan-cache")
}
if (length(a) >= 3L) {
  # absolute, because test_file() runs from the test directory
  rec <- file.path(normalizePath(dirname(a[3]), winslash = "/"),
                   basename(a[3]))
  unlink(rec)
  Sys.setenv(FRMTMB_BRMSPORT_RECORD = rec)
}
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
suppressMessages(library(p, character.only = TRUE))
# the uncapped reporter: the summary reporter stops at ten failures
r <- as.data.frame(test_file(f, package = p, env = testthat::test_env(p),
                             reporter = "silent", stop_on_failure = FALSE))
cat("RESULT ", p, " ", basename(f), " pass=", sum(r$passed), " fail=",
    sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
    " blocks=", nrow(r), "\n", sep = "")
fails <- r[r$failed > 0 | r$error, "test"]
if (length(fails)) cat("FAILING BLOCK:", fails, sep = "\n  ")
