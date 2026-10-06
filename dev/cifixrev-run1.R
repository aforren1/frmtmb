# Reviewer: run one test file of one package against a library stack.
# Usage: Rscript dev/cifixrev-run1.R <libs, ';'-separated> <package> <file>
# The user library is always appended last.
args <- commandArgs(TRUE)
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(strsplit(args[1], ";", fixed = TRUE)[[1]], user))
Sys.setenv(NOT_CRAN = "true")
p <- args[2]
f <- args[3]
Sys.setenv(FRMTMB_REV_TAG = paste0(p, "/", basename(f)))
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
