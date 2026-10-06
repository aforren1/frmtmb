# Reviewer: one test file with a verbose reporter, to read a failure.
#   Rscript dev/nanse-rev-onefile.R <merge|lane|release> <package> <file>
a <- commandArgs(trailingOnly = TRUE)
libs <- switch(a[1],
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"),
  release = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  base = "C:/Users/adf44/source/r/rellib-r5")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressMessages(library(testthat))
p <- a[2]
suppressMessages(library(p, character.only = TRUE))
cat("lib:", find.package("frmtmb"), "\n")
setwd(dirname(a[3]))
options(warn = 1)
r <- as.data.frame(test_file(a[3], package = p, env = test_env(p),
                             reporter = ProgressReporter$new(
                               show_praise = FALSE)))
cat("RESULT pass=", sum(r$passed), " fail=", sum(r$failed), " err=",
    sum(r$error), " skip=", sum(r$skipped), " warn=", sum(r$warning),
    "\n", sep = "")
