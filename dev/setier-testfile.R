# Lane setier: one test file against a library ("base" = rellib-r6).
#   Rscript dev/setier-testfile.R <lib or "base"> <test file> [package]
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(a[1], "base")) c(base, user) else c(a[1], base, user))
p <- if (length(a) >= 3) a[3] else "frmtmb"
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(p, character.only = TRUE)})
cat("lib:", find.package("frmtmb"), " BLAS probe",
    system.time({m <- matrix(1, 600, 600); m %*% m})[["elapsed"]], "\n")
res <- test_file(a[2], package = p, env = test_env(p), reporter = "silent")
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
