# Short review of lane alapse fix round 1: run one test file on the lane
# library and print a RESULT line.
#   Rscript dev/reviews/alapse-short-tf.R <test file name>
f <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-alapse-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
library(testthat)
library(frmtmb)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
t0 <- proc.time()
r <- testthat::test_file(file.path("tests/testthat", f), package = "frmtmb",
                         env = testthat::test_env("frmtmb"),
                         reporter = "silent")
d <- as.data.frame(r)
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d time=%.0fs\n",
            f, sum(d$passed), sum(d$failed), sum(d$error),
            sum(d$skipped), sum(d$warning), (proc.time() - t0)[[3]]))
bad <- d[d$failed > 0 | d$error | d$warning > 0, "test"]
if (length(bad)) cat("BAD:", bad, sep = "\n  ")
