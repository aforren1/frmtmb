# One test file per process. Usage:
#   Rscript frailty-runtest.R <file> [LIB]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) > 1L) args[2L] else
  "C:/Users/adf44/source/r/frailtylib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
setwd("C:/Users/adf44/source/r/frmtmb-wt-frailty/extensions/frmtmb.spline")
suppressMessages({
  library(testthat)
  library(frmtmb)
  library(frmtmb.spline)
})
cat("frmtmb.spline from", find.package("frmtmb.spline"), "\n")
# tests/testthat.R runs test_check(), which evaluates each file in an
# environment whose parent is the package NAMESPACE, so a test that
# reaches an internal function works there. A bare test_file() does
# not, and three files in this suite reach one; running them without
# this env reports an ERROR that is the runner's and not the code's.
res <- testthat::test_file(
  file.path("tests/testthat", args[1L]), reporter = "summary",
  env = testthat::test_env("frmtmb.spline"))
df <- as.data.frame(res)
cat("\nFILE", args[1L], "PASS", sum(df$passed), "FAIL", sum(df$failed),
    "ERROR", sum(df$error), "SKIP", sum(df$skipped), "WARN",
    sum(df$warning), "\n")
if (sum(df$failed) + sum(df$error) > 0) quit(status = 1L)
