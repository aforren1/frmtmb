# Run one test file against one arm. Usage:
# Rscript dev/reslot-test1.R before|after <path> [package]
args <- commandArgs(trailingOnly = TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
f <- args[2]
p <- if (length(args) >= 3) args[3] else "frmtmb"
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages(library(p, character.only = TRUE))
library(testthat)
r <- as.data.frame(testthat::test_file(f, package = p,
                                       env = testthat::test_env(p),
                                       reporter = "summary"))
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), sum(r$nb) - sum(r$failed) - sum(r$skipped),
            sum(r$failed), sum(r$error), sum(r$skipped), sum(r$warning)))
