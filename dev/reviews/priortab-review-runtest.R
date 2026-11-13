# Reviewer of lane priortab: run one test file on a library arm.
# Rscript priortab-review-runtest.R <pkg> <file> <lane|base>
a <- commandArgs(TRUE)
pkg <- a[1]; f <- a[2]; arm <- a[3]
libs <- c(if (arm == "lane") "C:/Users/adf44/source/r/wt-priortab-lib",
          "C:/Users/adf44/source/r/rellib-r20",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
library(testthat)
cat("lib:", dirname(find.package("frmtmb")), "\n")
r <- as.data.frame(testthat::test_file(f, package = pkg,
                                       env = testthat::test_env(pkg),
                                       reporter = "silent"))
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d\n", basename(f),
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped),
            sum(r$warning)))
bad <- r[r$failed > 0 | r$error, "test"]
if (length(bad)) cat("FAILING:", paste(unique(bad), collapse = " | "), "\n")
