# One test file against the lane's frmtmb.sample as it was before fix
# round 1 (the reviewer's read-only copy in wt-priordraws-rv-lib/floor),
# on the lane's current core: the fix-round tests seen to fail.
#   Rscript dev/priordraws-runtest-prev.R <test file>
f <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-priordraws-rv-lib/floor",
            "C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(testthat); library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), "|", find.package("frmtmb.sample"), "\n")
r <- as.data.frame(testthat::test_file(f, package = "frmtmb.sample",
                                       env = testthat::test_env("frmtmb.sample"),
                                       reporter = "summary"))
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d\n", basename(f),
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped)))
