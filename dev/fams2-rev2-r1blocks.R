# Reviewer, item 5: the lane's new test file against the reviewer's
# round-1 build (wt-fams2-rev-lib), per block.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-rev-lib", "C:/Users/adf44/source/r/rellib-r3", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages(library(testthat))
cat("frmtmb from", find.package("frmtmb"), "; lbeta_ad exists:", exists("lbeta_ad", asNamespace("frmtmb")), "\n")
r <- as.data.frame(test_file("C:/Users/adf44/source/r/frmtmb-wt-fams2/tests/testthat/test-xbeta-zibb-hurdle-cum.R",
  package = "frmtmb", env = test_env("frmtmb"), reporter = "silent", load_package = "installed"))
print(r[, c("test", "nb", "passed", "failed", "error", "warning")], right = FALSE)
