# One generated brms-suite test file, one R process, against a chosen
# library.
#
#   Rscript dev/brmsport-run.R <package> <test file> [record.tsv]
#
# With a third argument the harness RECORDS every assertion's outcome to
# that file and asserts nothing (dev/brmsport-ledger.R reads it); without
# one the file runs as the gated tier runs it, asserting the verdicts.
#
# THE LIBRARY AND THE TREE ARE ARGUMENTS, not constants. The first
# spelling of this file named `brmsport-lib` and the worktree
# `frmtmb-wt-brmsport`, neither of which outlived that lane: a later
# round could not run the tier at all without editing the runner. Set
# `FRMTMB_PORT_LIB` to the library holding the build under test; the
# round's reference build is the default.
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r3")
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

# The StanHeaders 2.32.10 pin is gone (2026-09-17); this file used to
# assert it. What replaces it is the user Makevars flag, because a
# populated Stan cache hides a compile failure everywhere except code
# that compiles something new. dev/release/run-tests.R carries the same
# two assertions, and the reason.
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
stopifnot(utils::packageVersion("tmbstan") >= "1.2.1",
          any(grepl("-std=gnu++17",
                    readLines(tools::makevars_user(), warn = FALSE),
                    fixed = TRUE)))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_BRMSPORT_PKG = p)
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                               winslash = "/",
                                               mustWork = FALSE))
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
