# Review of lane prilabel: one test file in one process on the lane
# library. args: package, test file path.
a <- commandArgs(TRUE)
p <- a[1]
f <- a[2]
.libPaths(c("C:/Users/adf44/source/r/wt-prilabel-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(testthat)
  library(p, character.only = TRUE)
})
cat("lib:", dirname(find.package("frmtmb")), "\n")
r <- as.data.frame(testthat::test_file(f, package = p,
                                       env = testthat::test_env(p),
                                       reporter = "silent"))
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d\n",
            basename(f), sum(r$passed), sum(r$failed), sum(r$error),
            sum(r$skipped), sum(r$warning)))
bad <- r[r$failed > 0 | r$error, "test"]
if (length(bad)) cat("BAD:", bad, sep = "\n  ")
