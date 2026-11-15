# Review of lane slopescale: run one core test file against the lane
# library, one file per process, and print a RESULT line.
# Usage: Rscript slopescale-review-testfile.R <test file name>
f <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-slopescale-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
})
cat("lib:", find.package("frmtmb"), "\nBLAS:", extSoftVersion()[["BLAS"]],
    "\nLAPACK:", La_library(), "\n")
p <- "frmtmb"
path <- file.path("C:/Users/adf44/source/r/frmtmb-wt-slopescale/tests/testthat",
                  f)
r <- as.data.frame(testthat::test_file(path, package = p,
                                       env = testthat::test_env(p),
                                       reporter = "silent"))
fails <- r[r$failed > 0 | r$error, c("test", "failed", "error")]
if (nrow(fails)) print(fails)
cat("RESULT", f, "pass=", sum(r$passed), "fail=", sum(r$failed),
    "err=", sum(r$error), "skip=", sum(r$skipped), "warn=", sum(r$warning),
    "\n")
