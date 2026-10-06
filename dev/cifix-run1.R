# Run one test file of one package against a chosen library stack.
# Usage: Rscript dev/cifix-run1.R <lib or "base"> <package> <test file>
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
libs <- if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user)
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true")
p <- args[2]
f <- args[3]
suppressPackageStartupMessages(library(p, character.only = TRUE))
cat("lib:", find.package("frmtmb"), " ", find.package(p), "\n")
cat("BLAS speed probe (s):",
    system.time({m <- matrix(1, 600, 600); m %*% m})[["elapsed"]], "\n")
res <- testthat::test_file(f, package = p, env = testthat::test_env(p),
                           reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
