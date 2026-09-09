# lane tmbstan: does the grad_log_prob block SKIP on a broken build?
# Run in isolation, because the fixture blocks in the same file reset
# the detector's memo on exit and would un-poison it.
.libPaths(c("C:/Users/adf44/source/r/lanelib-tmbstan",
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
library(testthat); library(frmtmb); library(frmtmb.sample)
setwd("C:/Users/adf44/source/r/frmtmb-wt-tmbstan/extensions/frmtmb.sample")
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("../../dev/stan-cache"))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) && identical(a[1L], "broken")) {
  assign("cached", TRUE,
         envir = environment(frmtmb.sample:::tmbstan_build_broken))
}
cat("### detector:", frmtmb.sample:::tmbstan_build_broken(), "\n")
r <- as.data.frame(test_file(
  "tests/testthat/test-tmbstan-build-guard.R", package = "frmtmb.sample",
  desc = "the reverse-mode gradient is the model's, not a standard normal",
  reporter = "silent"))
for (k in c("failed", "error", "skipped", "passed")) {
  cat("###  ", k, ": ", sum(r[[k]], na.rm = TRUE), "\n", sep = "")
}
