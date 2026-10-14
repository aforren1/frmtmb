# Does test-se-conv.R's "an sd held near zero short of its maximum is not
# at the maximum" pin lane illse's R2-1 fix? Runs that one test of the
# release tree's file on a library, as installed or with
# se_sd_gain_up() stubbed to 0 (the headroom without R2-1, as the
# release review's dev/relrev071-seconv.R removes it). The test must
# fail on 0.70.0, fail with the stub, and pass on the merged build.
#   Rscript dev/rel071-seconv3.R <lib> <installed|stub>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(frmtmb)})
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), format(packageVersion("frmtmb")),
    "| arm", a[2], "| BLAS", R.home(), "\n")
if (identical(a[2], "stub")) {
  stopifnot(exists("se_sd_gain_up", ns))
  unlockBinding("se_sd_gain_up", ns)
  assign("se_sd_gain_up", function(fit, k, ...) 0, envir = ns)
}
f <- "C:/Users/adf44/source/r/frmtmb-wt-release/tests/testthat/test-se-conv.R"
r <- test_file(f, package = "frmtmb", env = test_env("frmtmb"),
               reporter = "silent",
               desc = "an sd held near zero short of its maximum is not at the maximum")
for (t in r) for (x in t$results) {
  k <- if (inherits(x, "expectation_failure")) "FAIL" else
    if (inherits(x, "expectation_error")) "ERROR" else
      if (inherits(x, "expectation_success")) "ok" else class(x)[1]
  cat(k, "|", substr(gsub("\n", " ", conditionMessage(x)), 1, 160), "\n")
}
d <- as.data.frame(r)
cat(sprintf("RESULT pass=%d fail=%d err=%d\n", sum(d$nb) - sum(d$failed),
            sum(d$failed), sum(d$error)))
