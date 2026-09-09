# Reviewer instrument for lane tmbstan.
#
# Runs ONE test file in this process and emits one TSV row per
# test_that block, so the two arms can be compared block by block
# rather than by a summary line. The summary reporter caps failures at
# ten and a runner that sums `failed` alone misses a file that aborted,
# which is why every field is emitted per block.
#
# usage: Rscript dev/rev-tmbstan-sweep.R <pkgdir> <testfile> <arm>
#   arm: "clean" or "broken". "broken" pokes the detector's memo.
LIB <- "C:/Users/adf44/source/r/lanelib-tmbstan"
.libPaths(c(LIB, "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
args <- commandArgs(trailingOnly = TRUE)
pkgdir <- args[[1L]]
f <- args[[2L]]
arm <- args[[3L]]
pkg <- basename(pkgdir)

library(testthat)
library(frmtmb)
suppressMessages(library(frmtmb.sample))
root <- "C:/Users/adf44/source/r/frmtmb-wt-tmbstan"
Sys.setenv(FRMTMB_STAN_CACHE = file.path(root, "dev", "stan-cache"))
setwd(file.path(root, pkgdir))

env <- environment(frmtmb.sample:::tmbstan_build_broken)
if (identical(arm, "broken")) assign("cached", TRUE, envir = env)
cat("## detector:", frmtmb.sample:::tmbstan_build_broken(), "arm:", arm,
    "\n")

r <- as.data.frame(test_file(file.path("tests/testthat", f),
                             package = pkg, reporter = "silent"))
for (i in seq_len(nrow(r))) {
  cat(paste("ROW", arm, f, r$test[[i]], r$passed[[i]], r$failed[[i]],
            as.integer(r$error[[i]]), r$skipped[[i]], sep = "\t"), "\n",
      sep = "")
}
cat("## TOTAL", arm, f, nrow(r), sum(r$passed), sum(r$failed),
    sum(r$error), sum(r$skipped), "\n")
