# Mutation check for item 3.6: with the interval upper end NOT read, the
# interval test must fail; with the censored rows NOT tested at all, the
# no-events test and the flipped frailty assertions must fail. Runs the
# lane build with one internal function replaced in the namespace.
source("C:/Users/adf44/source/r/frmtmb-wt-phase3a/dev/phase3a-lib.R")
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(frmtmb); library(frmtmb.spline)})
setwd("C:/Users/adf44/source/r/frmtmb-wt-phase3a/extensions/frmtmb.spline")
ns <- asNamespace("frmtmb.spline")
orig <- get("sp_rp_fitted", ns)
run <- function(label, fun) {
  unlockBinding("sp_rp_fitted", ns)
  assign("sp_rp_fitted", fun, envir = ns)
  lockBinding("sp_rp_fitted", ns)
  res <- as.data.frame(testthat::test_file(
    "tests/testthat/test-rp-floored.R",
    reporter = testthat::SilentReporter$new(),
    env = testthat::test_env("frmtmb.spline")))
  cat(sprintf("%-28s %s\n", label, paste(sprintf("[%s] fail %d err %d",
    substr(res$test, 1, 40), res$failed, as.integer(res$error))[
      res$failed > 0 | res$error], collapse = "; ")))
}
no_hi <- function(object, fam) {
  f <- orig(object, fam)
  f$detadx_hi[] <- NA_real_
  f
}
# every censored row reported as monotone: the 0.7.0 behavior
no_cens <- function(object, fam) {
  f <- orig(object, fam)
  f$detadx[f$cens != 0] <- 1
  f$detadx_hi[] <- NA_real_
  f
}
run("unmutated", orig)
run("upper end not read", no_hi)
run("censored rows not tested", no_cens)
