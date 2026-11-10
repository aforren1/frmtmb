# Reviewer: run one lane test file against a build. args: file, lib
a <- commandArgs(TRUE)
.libPaths(c(strsplit(a[2], ";")[[1]], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(frmtmb)})
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
r <- as.data.frame(testthat::test_file(a[1], package = "frmtmb",
                   env = testthat::test_env("frmtmb"), reporter = "silent"))
bad <- r[r$failed > 0 | r$error, c("test", "failed", "error")]
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d\n", basename(a[1]),
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped), sum(r$warning)))
print(bad, row.names = FALSE)
