# Release reviewer of 0.71.0: one test file in one process, with the
# given library path (";"-separated, the user library appended), the
# package attached; prints the lib lines, each failure and the counts.
#   Rscript dev/relrev071-runtest.R <libs> <pkg> <file>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(a[1], ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
pkg <- a[2]
suppressMessages({library(testthat); library(pkg, character.only = TRUE)})
cat("lib:", pkg, as.character(packageVersion(pkg)), find.package(pkg),
    "| frmtmb", as.character(packageVersion("frmtmb")),
    find.package("frmtmb"), "\n")
res <- as.data.frame(test_file(a[3], package = pkg, env = test_env(pkg),
                               reporter = "summary",
                               stop_on_failure = FALSE))
cat(sprintf("RESULT %s pass %d fail %d error %d skip %d warn %d\n",
            basename(a[3]), sum(res$nb) - sum(res$failed),
            sum(res$failed), sum(res$error), sum(res$skipped),
            sum(res$warning)))
