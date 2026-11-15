# Short check of lane silent2's fix round: one test file in one process,
# on the lane library.
#   Rscript dev/reviews/silent2-short-runtest.R <package> <test file path>
args <- commandArgs(trailingOnly = TRUE)
pkg <- args[1]
f <- args[2]
.libPaths(c("C:/Users/adf44/source/r/wt-silent2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(testthat)
  library(pkg, character.only = TRUE)
})
cat("lib:", pkg, find.package(pkg), "frmtmb", find.package("frmtmb"),
    "R", R.home(), "\n")
t0 <- proc.time()[["elapsed"]]
r <- as.data.frame(testthat::test_file(f, package = pkg,
                                       env = testthat::test_env(pkg),
                                       reporter = "silent"))
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d time=%.0fs\n",
            basename(f), sum(r$passed), sum(r$failed), sum(r$error),
            sum(r$skipped), sum(r$warning), proc.time()[["elapsed"]] - t0))
bad <- r[r$failed > 0 | r$error | r$warning > 0 | r$skipped,
         c("test", "failed", "error", "warning", "skipped")]
if (nrow(bad)) print(bad)
