source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
Sys.setenv(FRMTMB_STAN_CACHE = file.path(REV_WT, "dev/review-ncgp-stan-cache"))
suppressPackageStartupMessages(library(testthat))
f <- "C:/Users/adf44/source/r/rev-ncgp-src/tests/test-brms-lik-gp-exact.R"
r <- testthat::test_file(f, package = "frmtmb", env = testthat::test_env("frmtmb"), reporter = "silent", stop_on_failure = FALSE)
for (t in r) if (grepl("gr = FALSE", t$test)) for (x in t$results) if (!inherits(x, "expectation_success")) cat(substr(conditionMessage(x), 1, 400), "\n---\n")
