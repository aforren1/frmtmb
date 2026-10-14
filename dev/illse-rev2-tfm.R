# Reviewer of lane illse: one test file on the trial merge illse + optstart.
#   Rscript dev/illse-rev-tf.R <file>
.libPaths(c("C:/Users/adf44/source/r/illse-rev-lib/merge", "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(frmtmb)})
f <- commandArgs(TRUE)[1]
cat("lib:", find.package("frmtmb"), "| R", R.home(), "\n")
res <- as.data.frame(test_file(f, package = "frmtmb", env = test_env("frmtmb"), reporter = "silent"))
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n", basename(f), sum(res$passed),
            sum(res$failed), sum(res$error), sum(res$skipped), sum(res$warning)))
