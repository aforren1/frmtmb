# Fix round 2: the rlddm false alarm, before and after type = "continuous" on a copy of
# frmtmb.learn. Arg: the library that holds the patched copy, or "lib2" for the shipped one.
a <- commandArgs(trailingOnly = TRUE)
libs <- c(if (a[1] != "lib2") a[1], "C:/Users/adf44/source/r/wt-nlfix-lib2", "C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(libs)
suppressPackageStartupMessages({ library(frmtmb); library(frmtmb.learn); library(testthat) })
cat("frmtmb", dirname(find.package("frmtmb")), " learn", dirname(find.package("frmtmb.learn")), " type", rlddm(subject = id, trial = trial)[["type"]], " bandit", bandit2arm_delta(subject = id, trial = trial)[["type"]], "\n")
set.seed(13003)
dl <- data.frame(id = rep(1:6, each = 40), trial = rep(1:40, times = 6), rt = 0.4 + rnorm(240, 0, 0.004), choice = rbinom(240, 1, 0.5), pay1 = rbinom(240, 1, 0.7), pay2 = rbinom(240, 1, 0.3))
w <- character(0)
r <- withCallingHandlers(frm(bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift ~ 1, bs ~ 1, ndt ~ 1, bias ~ 1), family = rlddm(subject = id, trial = trial), data = dl),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
cat("rlddm logLik", as.numeric(logLik(r)), " check fires:", any(grepl("every response is discrete", w, fixed = TRUE)), "\n")
if (length(a) > 1) {
  for (f in list.files(a[2], "^test-.*[.]R$", full.names = TRUE)) {
    x <- as.data.frame(test_file(f, package = "frmtmb.learn", env = test_env("frmtmb.learn"), reporter = "silent", stop_on_failure = FALSE))
    cat(sprintf("FILE %-34s pass=%d fail=%d err=%d skip=%d\n", basename(f), sum(x$passed), sum(x$failed), sum(x$error), sum(x$skipped)))
  }
}
