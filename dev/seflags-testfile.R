# Lane seflags: run one test file of one package against the lane build
# ("lane": wt-seflags-lib over rellib-r7), the 0.69.0 base ("base":
# rellib-r7) or the trial merge with lane optrap ("merge",
# dev/seflags-merge.sh, rellib-r7 behind it), one file per process;
# prints every failure, error, skip and escaped warning, then the
# counts.
#   Rscript dev/seflags-testfile.R <lane|base|merge> <test file> [package]
a <- commandArgs(trailingOnly = TRUE)
libs <- switch(a[1],
  lane = c("C:/Users/adf44/source/r/wt-seflags-lib",
           "C:/Users/adf44/source/r/rellib-r7"),
  base = "C:/Users/adf44/source/r/rellib-r7",
  merge = c("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-merge-out/lib",
            "C:/Users/adf44/source/r/rellib-r7"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
pkg <- if (length(a) >= 3) a[3] else "frmtmb"
suppressMessages({
  library(testthat)
  library(pkg, character.only = TRUE)
})
cat("lib:", find.package(pkg), "| frmtmb:", find.package("frmtmb"),
    "| R:", R.home(), "| BLAS probe",
    system.time({m <- matrix(1, 600, 600); m %*% m})[["elapsed"]], "\n")
res <- test_file(a[2], package = pkg, env = test_env(pkg),
                 reporter = "silent")
for (t in res) for (x in t$results) {
  k <- if (inherits(x, "expectation_failure")) "FAIL" else
    if (inherits(x, "expectation_error")) "ERROR" else
      if (inherits(x, "expectation_skip")) "SKIP" else
        if (inherits(x, "expectation_warning")) "WARN" else NA
  if (!is.na(k)) {
    cat("DETAIL", k, "[", t$test, "]",
        substr(gsub("\n", " | ", conditionMessage(x)), 1, 400), "\n")
  }
}
df <- as.data.frame(res)
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d\n",
            basename(a[2]), sum(df$passed), sum(df$failed), sum(df$error),
            sum(df$skipped), sum(df$warning)))
