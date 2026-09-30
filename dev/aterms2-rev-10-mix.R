# Reviewer, claim 8: the sample test with the lane's core and the BASE
# frmtmb.sample (a copy of rellib-r3's installed package in the session
# scratchpad, first on the path). Log: dev/aterms2-rev-log-10-mix.txt
.libPaths(c("C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad/aterms2-rev-mixlib",
            "C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(testthat); library(frmtmb.sample)})
cat("core", find.package("frmtmb"), " sample", find.package("frmtmb.sample"), "\n")
f <- "C:/Users/adf44/source/r/frmtmb-wt-aterms2/extensions/frmtmb.sample/tests/testthat/test-subset-rate-draws.R"
setwd(dirname(f))
options(testthat.summary.max_reports = 10000L)
res <- testthat::test_file(f, package = "frmtmb.sample",
                           env = testthat::test_env("frmtmb.sample"),
                           reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("\nRESULT mix pass=%d fail=%d error=%d skip=%d warn=%d\n",
            sum(df$passed), sum(df$failed), sum(df$error), sum(df$skipped),
            sum(df$warning)))
