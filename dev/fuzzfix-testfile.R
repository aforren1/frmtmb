# Run one test file against a library and print one RESULT line.
#   Rscript dev/fuzzfix-testfile.R <lib|base> <pkg> <test file> [label]
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
Sys.setenv(NOT_CRAN = "true")
p <- a[2]; f <- a[3]; lbl <- if (length(a) > 3) a[4] else a[1]
suppressMessages({library(testthat); library(p, character.only = TRUE)})
cat("R:", R.home(), "| LAPACK:", La_version(), "\n")
cat("lib:", p, format(packageVersion(p)), find.package(p),
    "| frmtmb", format(packageVersion("frmtmb")), find.package("frmtmb"),
    "\n")
res <- testthat::test_file(f, package = p, env = testthat::test_env(p),
                           reporter = "check", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("RESULT %s %s pass=%d fail=%d err=%d skip=%d warn=%d\n", lbl,
            basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
