# Short check: test-reparam.R on the lane's library, one process.
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib","C:/Users/adf44/source/r/rellib-r16","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(testthat); library(frmtmb.sample)})
res <- as.data.frame(test_file("extensions/frmtmb.sample/tests/testthat/test-reparam.R", package="frmtmb.sample", env=test_env("frmtmb.sample"), reporter="silent", stop_on_failure=FALSE))
cat(sprintf("RERUN test-reparam.R pass=%d fail=%d error=%d skip=%d warn=%d\n", sum(res$passed),sum(res$failed),sum(as.integer(res$error)),sum(as.integer(res$skipped)),sum(res$warning)))
