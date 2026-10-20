# The frmtmb.sample half of lane catmulti seen to fail alone: the lane's
# core with rellib-r9's frmtmb.sample 0.19.0, then the lane's own
# frmtmb.sample, on the two frmtmb.sample tests of the lane.
#   Rscript dev/catmulti-mixarm.R <sample: base|lane>
# Output: dev/catmulti-out/mixarm-<arm>.txt
arm <- commandArgs(TRUE)[1]
lib <- "C:/Users/adf44/source/r/wt-catmulti-lib"
base <- "C:/Users/adf44/source/r/rellib-r9"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(lib, base, user))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb, lib.loc = lib)
  library(frmtmb.sample, lib.loc = if (arm == "base") base else lib)
  library(testthat)
})
cat("frmtmb", dirname(find.package("frmtmb")), "| frmtmb.sample",
    dirname(find.package("frmtmb.sample")), "\n")
for (f in c("test-linpred-stack.R", "test-simplex-mvbf-draws.R")) {
  res <- test_file(file.path("extensions/frmtmb.sample/tests/testthat", f),
                   package = "frmtmb.sample",
                   env = test_env("frmtmb.sample"), reporter = "summary",
                   stop_on_failure = FALSE)
  df <- as.data.frame(res)
  cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d\n", f,
              sum(df$passed), sum(df$failed), sum(df$error),
              sum(df$skipped)))
}
