# One test file, in one R process, with the failure cap LIFTED and the
# package named so that internals resolve.
#
#   Rscript dev/rlddm-scripts/rlddm-testfile.R <lib> <pkg> <file> [env=...]
#
# The cap matters: testthat's summary reporter stops at ten failures and
# says so, and a lane in this repository has already reported a capped
# number as a count. `err` is printed beside `fail` for the same reason:
# a file that aborts halfway reports a clean failure count.

args <- commandArgs(trailingOnly = TRUE)
lib <- args[[1L]]
pkg <- args[[2L]]
file <- args[[3L]]
for (kv in args[-(1:3)]) {
  p <- strsplit(kv, "=", fixed = TRUE)[[1L]]
  do.call(Sys.setenv, stats::setNames(list(p[[2L]]), p[[1L]]))
}

.libPaths(c(lib,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", TESTTHAT_MAX_FAILS = "10000")

# ONLY the package under test is attached, which is what its own
# tests/testthat.R does. Attaching its siblings here would put their
# exports on the search path and hide exactly the defect R CMD check
# found in this lane: five assertions in test-rlddm-ndt.R called
# `ndt_time()` unqualified, which resolves under a runner that attaches
# frmtmb.eam and does not resolve under `R CMD check`.
suppressPackageStartupMessages({
  library(testthat)
  library(pkg, character.only = TRUE)
})
cat("lib :", lib, "\n")
cat("pkg :", pkg, format(packageVersion(pkg)), "at",
    dirname(system.file(package = pkg)), "\n")
cat("file:", file, "\n")

res <- testthat::test_file(
  file, package = pkg,
  reporter = testthat::SummaryReporter$new(max_reports = 10000L))
df <- as.data.frame(res)
cat("\nRESULT file=", basename(file),
    " pass=", sum(df$passed), " fail=", sum(df$failed),
    " err=", sum(df$error), " skip=", sum(df$skipped),
    " tests=", nrow(df), "\n", sep = "")
