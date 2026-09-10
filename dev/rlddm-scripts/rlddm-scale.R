# One scale row, in a fresh process.
#
# Usage (from the worktree root):
#   Rscript dev/rlddm-scripts/rlddm-scale.R <row> <lib> <out.tsv>
#
# `lib` selects the arm: the lane library holds this change, and
# `C:/Users/adf44/source/r/rellib-0552` is the round's shared reference
# build of the base commit, which supplies the BEFORE numbers.
#
# The tier's own seed is 20260908, stated in
# extensions/frmtmb.learn/tests/testthat/test-scale.R.

args <- commandArgs(trailingOnly = TRUE)
row <- args[[1L]]
lib <- args[[2L]]
out <- args[[3L]]

.libPaths(c(lib,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

Sys.setenv(FRMTMB_SCALE_TESTS = "true", NOT_CRAN = "true",
           FRMTMB_SCALE_ROW = row, FRMTMB_SCALE_OUT = out)

cat("row  :", row, "\n")
cat("lib  :", lib, "\n")
cat("learn:", format(packageVersion("frmtmb.learn")), "at",
    dirname(system.file(package = "frmtmb.learn")), "\n")
cat("eam  :", format(packageVersion("frmtmb.eam")), "at",
    dirname(system.file(package = "frmtmb.eam")), "\n")

suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
  library(frmtmb.eam)
  library(frmtmb.learn)
})

t0 <- Sys.time()
res <- testthat::test_file(
  "extensions/frmtmb.learn/tests/testthat/test-scale.R",
  package = "frmtmb.learn",
  reporter = testthat::MultiReporter$new(list(
    testthat::SummaryReporter$new(max_reports = 1000L),
    testthat::CheckReporter$new())))
cat("\nwall:", format(as.numeric(difftime(Sys.time(), t0, units = "secs")),
                      digits = 6), "s\n")
df <- as.data.frame(res)
cat("FILE test-scale.R  pass=", sum(df$passed), " fail=", sum(df$failed),
    " err=", sum(df$error), " skip=", sum(df$skipped), "\n", sep = "")
