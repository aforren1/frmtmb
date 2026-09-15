## One test file per R process, which is the standing rule.
##
##   COH_TEST=test-coherence.R Rscript dev/coh-runtest.R
##
## `COH_DESC` runs one test_that() block of it while iterating.
.libPaths(c("C:/Users/adf44/source/r/coh-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.coupling))
setwd("extensions/frmtmb.coupling/tests/testthat")
f <- Sys.getenv("COH_TEST", "test-coherence.R")
desc <- Sys.getenv("COH_DESC", "")
t0 <- Sys.time()
res <- if (nzchar(desc)) {
  testthat::test_file(f, desc = desc, reporter = "summary")
} else {
  testthat::test_file(f, reporter = "summary")
}
df <- as.data.frame(res)
## Sum ERROR as well as failed: a runner that adds only `failed` prints
## a clean line for a file that aborted halfway, and one did.
cat(sprintf("\nFILE %s  PASS=%d FAIL=%d ERROR=%d SKIP=%d WARN=%d  %.1fs\n",
            f, sum(df$passed), sum(df$failed), sum(df$error),
            sum(df$skipped), sum(df$warning),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
