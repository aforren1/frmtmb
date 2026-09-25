# Punch round 1: the new and changed test-records.R run against the
# ROUND-0 frm_ode_records() (before the punch), recovered from the check
# tarball built for the round-0 R CMD check
# (C:/Users/adf44/source/r/phase3a-round1/frmtmb.ode/R/records.R), put
# into the lane namespace in place of the current one. The released
# build has no frm_ode_records() at all, so this is the run that shows
# each new assertion failing on the code the review found wanting.
source("C:/Users/adf44/source/r/frmtmb-wt-phase3a/dev/phase3a-lib.R")
Sys.setenv(NOT_CRAN = "true")
suppressMessages({
  library(testthat); library(frmtmb); library(frmtmb.ode)
})
setwd("C:/Users/adf44/source/r/frmtmb-wt-phase3a/extensions/frmtmb.ode")
ns <- asNamespace("frmtmb.ode")
old <- new.env(parent = ns)
sys.source("C:/Users/adf44/source/r/phase3a-round1/frmtmb.ode/R/records.R",
           envir = old)
for (nm in c("frm_ode_records", "ode_rec_column")) {
  unlockBinding(nm, ns)
  assign(nm, get(nm, envir = old), envir = ns)
  lockBinding(nm, ns)
}
res <- as.data.frame(testthat::test_file(
  "tests/testthat/test-records.R",
  reporter = testthat::SilentReporter$new(),
  env = testthat::test_env("frmtmb.ode")))
for (i in seq_len(nrow(res))) {
  cat(sprintf("%-62s pass %3d fail %3d error %s\n", substr(res$test[i], 1, 62),
              res$passed[i], res$failed[i], res$error[i]))
}
cat(sprintf("TOTAL blocks %d pass %d fail %d error %d\n", nrow(res),
            sum(res$passed), sum(res$failed), sum(res$error)))
