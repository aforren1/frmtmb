# Reviewer, fix round 1: one frmtmb.eam test file on a given core.
#   Rscript cdfs-rev-eamtf.R <arm> <file>; frmtmb.eam from rellib-r10
a <- commandArgs(TRUE)
lp <- c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (a[1] == "lane") lp <- c("C:/Users/adf44/source/r/wt-cdfs-rev-lib", lp)
.libPaths(lp)
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(frmtmb); library(frmtmb.eam)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.eam"), "\n")
td <- "C:/Users/adf44/source/r/frmtmb-wt-cdfs/extensions/frmtmb.eam/tests/testthat"
nw <- 0L
r <- withCallingHandlers(
  as.data.frame(test_file(file.path(td, a[2]), package = "frmtmb.eam",
                          env = testthat::test_env("frmtmb.eam"), reporter = "silent",
                          stop_on_failure = FALSE)),
  warning = function(w) { nw <<- nw + 1L })
cat(sprintf("RESULT %s %s pass=%d fail=%d err=%d skip=%d warn=%d\n", a[1], a[2],
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped), nw))
bad <- r[r$failed > 0 | r$error, "test"]
if (length(bad)) cat("FAILING:", bad, sep = "\n  ")
cat("\n")
