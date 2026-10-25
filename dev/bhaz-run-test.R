# One test file in one process, against a library stack.
#   Rscript dev/bhaz-run-test.R <pkg> <file> <lane|base>
# lane: the bhaz lane library, then rellib-r12, then the user library.
# base: rellib-r12, then the user library (the "before" arm).
args <- commandArgs(trailingOnly = TRUE)
pkg <- args[[1L]]
f <- args[[2L]]
arm <- if (length(args) >= 3L) args[[3L]] else "lane"
libs <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-bhaz-lib", libs)
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
cat("lib:", find.package("frmtmb"), " ", pkg, "at", find.package(pkg),
    "\n")
res <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                           reporter = testthat::ProgressReporter$new(
                             show_praise = FALSE))
df <- as.data.frame(res)
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
