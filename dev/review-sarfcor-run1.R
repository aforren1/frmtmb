# Reviewer: run one test file against a given library stack.
# Rscript dev/review-sarfcor-run1.R <lib> <pkg> <file>
args <- commandArgs(TRUE)
lib <- args[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
p <- args[2]
f <- args[3]
suppressMessages(library(p, character.only = TRUE))
cat("lib:", p, find.package(p), format(packageVersion(p)), "| frmtmb",
    find.package("frmtmb"), "\n")
r <- as.data.frame(testthat::test_file(f, package = p,
                                       env = testthat::test_env(p),
                                       reporter = "summary"))
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d\n", basename(f),
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped),
            sum(r$warning)))
