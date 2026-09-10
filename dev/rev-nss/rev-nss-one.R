# One test file, one R process. A whole-suite run in one process has
# repeatedly hidden state leakage here, and a runner that sums `failed`
# and not `error` prints a clean line for a file that aborted halfway,
# so both are summed and printed.
#
# Usage: Rscript rev-nss-one.R <test-file> [lane|ref]
args <- commandArgs(trailingOnly = TRUE)
Sys.setenv(REV_ARM = if (length(args) > 1L) args[[2L]] else "lane")
source("C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss/prelude.R")
suppressMessages({
  library(testthat); library(frmtmb); library(frmtmb.ode)
})
setwd(paste0("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/",
             "frmtmb.ode/tests/testthat"))
# package = : without it local_mocked_bindings() cannot find the
# namespace and reports "No packages loaded with pkgload", which shows
# up as an ERROR that belongs to the runner and not to the code
r <- test_file(args[[1L]], reporter = "silent",
               package = "frmtmb.ode")
d <- as.data.frame(r)
cat(sprintf("REV %s pass=%d fail=%d err=%d warn=%d skip=%d\n",
            args[[1L]], sum(d$passed), sum(d$failed), sum(d$error),
            sum(d$warning), sum(d$skipped)))
