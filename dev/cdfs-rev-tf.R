# Reviewer: one test file, one process. Rscript cdfs-rev-tf.R <arm> <file>
a <- commandArgs(TRUE)
lp <- c("C:/Users/adf44/source/r/rellib-r10",
        "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (a[1] == "lane") lp <- c("C:/Users/adf44/source/r/wt-cdfs-rev-lib", lp)
.libPaths(lp)
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(frmtmb)})
cat("lib:", find.package("frmtmb"), "\n")
f <- file.path("C:/Users/adf44/source/r/frmtmb-wt-cdfs/tests/testthat", a[2])
nw <- 0L
r <- withCallingHandlers(
  as.data.frame(test_file(f, package = "frmtmb", env = testthat::test_env("frmtmb"),
                          reporter = "silent", stop_on_failure = FALSE)),
  warning = function(w) { nw <<- nw + 1L })
cat(sprintf("RESULT %s %s pass=%d fail=%d err=%d skip=%d warn=%d\n", a[1], a[2],
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped), nw))
bad <- r[r$failed > 0 | r$error, "test"]
if (length(bad)) cat("FAILING:", bad, sep = "\n  ")
cat("\n")
