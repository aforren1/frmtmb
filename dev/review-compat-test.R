# Reviewer scratch: run the lane's test-compat.R against a library list.
#   Rscript dev/review-compat-test.R <lib;lib> <file>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(args[1], ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
cat("lib:", find.package("frmtmb"), "\n")
r <- as.data.frame(testthat::test_file(args[2], package = "frmtmb",
                                       env = testthat::test_env("frmtmb"),
                                       reporter = "silent"))
cat(sprintf("pass=%d fail=%d err=%d skip=%d warn=%d\n", sum(r$passed),
            sum(r$failed), sum(r$error), sum(r$skipped), sum(r$warning)))
bad <- r[r$failed > 0 | r$error, c("test", "failed")]
print(bad, row.names = FALSE)
