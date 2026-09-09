# lane tmbstan: the one frmtmb.eam block that reaches frm_sample().
LIB <- "C:/Users/adf44/source/r/lanelib-tmbstan"
.libPaths(c(LIB, "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
library(testthat); library(frmtmb); library(frmtmb.eam)
setwd("C:/Users/adf44/source/r/frmtmb-wt-tmbstan/extensions/frmtmb.eam")
a <- commandArgs(trailingOnly = TRUE)
if (length(a) && identical(a[1L], "broken")) {
  loadNamespace("frmtmb.sample")
  assign("cached", TRUE,
         envir = environment(frmtmb.sample:::tmbstan_build_broken))
  cat("### arm: BROKEN\n")
} else cat("### arm: clean\n")
r <- as.data.frame(test_file("tests/testthat/test-sampling.R",
                             package = "frmtmb.eam", reporter = "silent"))
cat("### blocks:", nrow(r), "\n")
for (k in c("failed", "error", "skipped", "passed")) {
  cat("###  ", k, ": ", sum(r[[k]], na.rm = TRUE), "\n", sep = "")
}
