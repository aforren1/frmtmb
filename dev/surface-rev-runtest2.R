# Reviewer copy of dev/surface-runtest.R (reviewer library). Run one test file of one package against one build, in
# its own process, and print a RESULT line from every expectation.
#
#   Rscript dev/surface-runtest.R <pkg> <test file path> lane|base
#
# Every gate on: NOT_CRAN, FRMTMB_BRMS_FIT_TESTS, FRMTMB_DRMTMB_FIT_TESTS,
# FRMTMB_FUZZ.
a <- commandArgs(TRUE)
pkg <- a[1]
f <- a[2]
arm <- a[3]
source("C:/Users/adf44/source/r/frmtmb-wt-surface/dev/surface-rev-env.R")
rev_env(arm)
Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", FRMTMB_DRMTMB_FIT_TESTS = "true",
           FRMTMB_FUZZ = "true")
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
cat("lib:", find.package(pkg), "| frmtmb", find.package("frmtmb"), "\n")
res <- testthat::test_file(f, package = pkg,
                           env = testthat::test_env(pkg),
                           reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
fails <- 0L
for (i in seq_len(nrow(df))) {
  r <- res[[i]]$results
  for (x in r) {
    cl <- class(x)[1]
    if (cl %in% c("expectation_failure", "expectation_error")) {
      cat("FAIL [", df$test[i], "]:", substr(conditionMessage(x), 1, 400),
          "\n")
    }
    if (cl == "expectation_warning") {
      cat("WARN [", df$test[i], "]:", substr(conditionMessage(x), 1, 200),
          "\n")
    }
    if (cl == "expectation_skip") {
      cat("SKIP [", df$test[i], "]:", substr(conditionMessage(x), 1, 200),
          "\n")
    }
  }
}
cat(sprintf("RESULT %s %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            pkg, basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
