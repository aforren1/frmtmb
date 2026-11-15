# Review of lane slopescale, attack items 4 and 9: the lane merged with
# lane selimit2 (owner of R/se-check.R), loaded from source with
# pkgload::load_all() (nothing installed). Runs one test file of the
# merged tree and prints a RESULT line.
# Usage: Rscript slopescale-review-merged.R <merged tree> <test file>
a <- commandArgs(TRUE)
tree <- a[1]
f <- a[2]
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(testthat)
  pkgload::load_all(tree, quiet = TRUE, export_all = FALSE,
                    helpers = FALSE, attach_testthat = FALSE)
})
cat("frmtmb from", getNamespaceInfo("frmtmb", "path"), "\n")
r <- as.data.frame(testthat::test_file(
  file.path(tree, "tests/testthat", f), package = "frmtmb",
  env = testthat::test_env("frmtmb"), reporter = "silent",
  load_package = "none"))
fails <- r[r$failed > 0 | r$error, c("test", "failed", "error")]
if (nrow(fails)) print(fails)
cat("RESULT merged", f, "pass=", sum(r$passed), "fail=", sum(r$failed),
    "err=", sum(r$error), "skip=", sum(r$skipped), "warn=", sum(r$warning),
    "\n")
