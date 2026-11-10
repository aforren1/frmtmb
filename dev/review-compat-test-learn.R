# Reviewer scratch (round 2): one frmtmb.learn test file on the lane core.
#   Rscript dev/review-compat-test-learn.R <file>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-compat-lib",
            "C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
p <- "frmtmb.learn"
suppressPackageStartupMessages({library(testthat); library(p, character.only = TRUE)})
cat("lib:", find.package("frmtmb"), "|", find.package(p), "\n")
r <- as.data.frame(testthat::test_file(args[1], package = p,
                                       env = testthat::test_env(p),
                                       reporter = "silent"))
cat(sprintf("%s pass=%d fail=%d err=%d skip=%d warn=%d\n", basename(args[1]),
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped),
            sum(r$warning)))
