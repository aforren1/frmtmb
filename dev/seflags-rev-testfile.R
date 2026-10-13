# Reviewer of lane seflags: one test file against a trial merge build
# (seflags + moparity: "mo"; seflags + optrap: "op"), one file per
# process, counts and every non-pass.
#   Rscript dev/seflags-rev-testfile.R <mo|op> <test file> [package]
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(paste0("C:/Users/adf44/source/r/seflags-rev-lib/", a[1]),
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
pkg <- if (length(a) >= 3) a[3] else "frmtmb"
suppressMessages({
  library(testthat)
  library(pkg, character.only = TRUE)
})
cat("lib:", find.package(pkg), "| frmtmb:", find.package("frmtmb"), "\n")
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
