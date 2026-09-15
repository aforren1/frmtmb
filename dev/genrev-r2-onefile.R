# genrev round 2: reproduce the single-process regression WITHOUT
# R CMD check.  Run ONE test file in one process, with brms's namespace
# loaded first (which is what an earlier file did under R CMD check),
# or not.  Counts failures AND errors.
a <- commandArgs(trailingOnly = TRUE)
LIB <- a[1]; FILE <- a[2]; pre <- identical(a[3], "brms")
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
if (pre) suppressMessages(loadNamespace("brms"))
library(testthat)
suppressMessages(library(frmtmb))
setup <- file.path(dirname(FILE), "setup.R")
if (file.exists(setup)) sys.source(setup, envir = globalenv())
res <- as.data.frame(test_file(FILE, reporter = "silent", package = "frmtmb"))
cat(sprintf("LIB %s brms-preloaded %s\n", basename(LIB), pre))
cat(sprintf("BLOCKS %d ASSERT %d PASS %d FAIL %d ERROR %d SKIP %d\n",
            nrow(res), sum(res$passed) + sum(res$failed) + sum(res$error),
            sum(res$passed), sum(res$failed), sum(res$error), sum(res$skipped)))
cat("GENREVDONE\n")
