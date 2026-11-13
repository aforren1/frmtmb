# Release review of 0.82.0: one test file on rellib-r21, one process.
# usage: rel082-review-test1.R <package> <test file path>
args <- commandArgs(TRUE)
p <- args[1]
f <- args[2]
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(p, character.only = TRUE))
cat("lib:", find.package(p), as.character(packageVersion(p)),
    "| frmtmb", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
t0 <- proc.time()[["elapsed"]]
r <- as.data.frame(testthat::test_file(f, package = p,
                                       env = testthat::test_env(p),
                                       reporter = "silent",
                                       stop_on_failure = FALSE))
el <- proc.time()[["elapsed"]] - t0
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d secs=%.0f\n",
            basename(f), sum(r$passed), sum(r$failed), sum(r$error),
            sum(r$skipped), sum(r$warning), el))
bad <- r[r$failed > 0 | r$error | r$warning > 0, c("test", "failed", "error",
                                                   "warning")]
if (nrow(bad)) print(bad)
