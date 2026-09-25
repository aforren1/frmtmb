.libPaths(unique(c( "C:/Users/adf44/source/r/rellib-r3",
  "C:/Users/adf44/source/r/pinlib",
  "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
Sys.setenv(NOT_CRAN = "true")
a <- commandArgs(TRUE)
suppressPackageStartupMessages({library(testthat); library(a[2], character.only = TRUE)})
cat("from", find.package(a[2]), find.package("frmtmb.eam"), "\n")
r <- as.data.frame(test_file(a[1], reporter = "progress", package = a[2],
                             load_package = "installed"))
print(r[, c("test", "nb", "failed", "skipped", "error", "passed")])
cat("RESULT", basename(a[1]), "blocks", nrow(r), "expectations", sum(r$nb),
    "pass", sum(r$passed), "fail", sum(r$failed), "err", sum(r$error),
    "skip", sum(r$skipped), "\n")
