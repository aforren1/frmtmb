source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("base")
suppressPackageStartupMessages(library(testthat))
f <- file.path(REV_WT, "tests/testthat/test-gp-noncentered.R")
r <- test_file(f, package = "frmtmb", env = test_env("frmtmb"), reporter = "silent", stop_on_failure = FALSE)
for (tt in r) { res <- tt$results; nf <- sum(vapply(res, inherits, NA, "expectation_failure")); ne <- sum(vapply(res, inherits, NA, "expectation_error")); ok <- sum(vapply(res, inherits, NA, "expectation_success"))
  msg <- ""; if (ne) msg <- substr(conditionMessage(Filter(function(x) inherits(x, "expectation_error"), res)[[1]]), 1, 90)
  cat(sprintf("RSF %-58s ok %d fail %d err %d %s\n", substr(tt$test, 1, 58), ok, nf, ne, gsub("\n", " ", msg))) }
