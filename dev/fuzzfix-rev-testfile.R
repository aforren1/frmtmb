# Reviewer: one test file against a library, printing the BLAS and
# LAPACK in use so the log shows which R ran it.
#   Rscript dev/fuzzfix-rev-testfile.R <lib|base> <test file> <label>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(frmtmb)})
cat("R:", R.home(), "| BLAS:", extSoftVersion()[["BLAS"]], "| LAPACK:",
    La_library(), La_version(), "\n")
cat("lib: frmtmb", format(packageVersion("frmtmb")), find.package("frmtmb"),
    "\n")
res <- testthat::test_file(a[2], package = "frmtmb",
                           env = testthat::test_env("frmtmb"),
                           reporter = "check", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("RESULT %s %s pass=%d fail=%d err=%d skip=%d warn=%d\n", a[3],
            basename(a[2]), sum(df$passed), sum(df$failed), sum(df$error),
            sum(df$skipped), sum(df$warning)))
