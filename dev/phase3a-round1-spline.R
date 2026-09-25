# Punch round 2 of item 3.6: the new and changed spline tests run
# against the ROUND-1 rp-check.R and royston-parmar.R (copied before the
# round-2 edit to C:/Users/adf44/source/r/phase3a-round1/frmtmb.spline/R/),
# put into the lane namespace in place of the current functions. Round 1
# refused on censored rows; the tests now expect a warning, so each must
# fail there.
source("C:/Users/adf44/source/r/frmtmb-wt-phase3a/dev/phase3a-lib.R")
Sys.setenv(NOT_CRAN = "true")
suppressMessages({
  library(testthat); library(frmtmb); library(frmtmb.spline)
})
setwd("C:/Users/adf44/source/r/frmtmb-wt-phase3a/extensions/frmtmb.spline")
ns <- asNamespace("frmtmb.spline")
old <- new.env(parent = ns)
src <- "C:/Users/adf44/source/r/phase3a-round1/frmtmb.spline/R/"
for (f in c("rp-check.R", "royston-parmar.R")) {
  sys.source(paste0(src, f), envir = old)
}
for (nm in ls(old)) {
  if (!exists(nm, envir = ns, inherits = FALSE)) next
  unlockBinding(nm, ns)
  assign(nm, get(nm, envir = old), envir = ns)
  lockBinding(nm, ns)
}
# the exported copies too, so a test calling rp_floored() by name reaches
# the round-1 one
exp_env <- as.environment("package:frmtmb.spline")
for (nm in intersect(ls(old), ls(exp_env))) {
  unlockBinding(nm, exp_env)
  assign(nm, get(nm, envir = old), envir = exp_env)
  lockBinding(nm, exp_env)
}
for (tf in c("test-rp-floored.R", "test-frailty.R")) {
  res <- as.data.frame(testthat::test_file(
    file.path("tests/testthat", tf),
    reporter = testthat::SilentReporter$new(),
    env = testthat::test_env("frmtmb.spline")))
  cat("\n==", tf, "\n")
  for (i in seq_len(nrow(res))) {
    cat(sprintf("%-62s pass %3d fail %3d error %s\n",
                substr(res$test[i], 1, 62), res$passed[i], res$failed[i],
                res$error[i]))
  }
  cat(sprintf("TOTAL blocks %d pass %d fail %d error %d\n", nrow(res),
              sum(res$passed), sum(res$failed), sum(res$error)))
}
