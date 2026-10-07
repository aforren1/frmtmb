# Reviewer of lane surface: one test file, one process.
#   Rscript dev/surface-rev-runtest.R <arm> <pkg> <file>
a <- commandArgs(TRUE)
source("dev/surface-rev-env.R")
rev_env(a[1])
Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")
pkg <- a[2]
suppressPackageStartupMessages({
  library(testthat); library(pkg, character.only = TRUE)
})
tdir <- if (pkg == "frmtmb") "tests/testthat" else
  file.path("extensions", pkg, "tests/testthat")
res <- as.data.frame(test_file(file.path(tdir, a[3]), package = pkg,
                               env = test_env(pkg), reporter = "silent",
                               stop_on_failure = FALSE))
cat(sprintf("RESULT %s %s %s pass=%d fail=%d error=%d skip=%d warn=%d | lib: %s\n",
            a[1], pkg, a[3], sum(res$passed), sum(res$failed),
            sum(res$error), sum(res$skipped), sum(res$warning),
            find.package("frmtmb")))
