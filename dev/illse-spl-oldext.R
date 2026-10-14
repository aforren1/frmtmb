# Lane illse: frmtmb.spline's test-illcond.R with the lane core and the
# 0.70.0 frmtmb.spline (rellib-r8), to see it fail where the extension
# formed C V C' from the whole matrix. The lane core is a copy in
# dev/illse-tmplib so that the lane library's own frmtmb.spline is not
# on the path.
.libPaths(c("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-tmplib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(frmtmb.spline)})
cat("frmtmb from", dirname(find.package("frmtmb")), "| frmtmb.spline from",
    dirname(find.package("frmtmb.spline")), "\n")
res <- as.data.frame(test_file(
  "C:/Users/adf44/source/r/frmtmb-wt-illse/extensions/frmtmb.spline/tests/testthat/test-illcond.R",
  package = "frmtmb.spline", env = test_env("frmtmb.spline"),
  reporter = "summary"))
cat(sprintf("RESULT pass=%d fail=%d error=%d\n", sum(res$passed),
            sum(res$failed), sum(res$error)))
