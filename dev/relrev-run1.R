# Reviewer: run one test file on rellib-r6, the package attached, one
# process per file. A third argument "globalf" defines a global `f` (the
# test file's path, as the 0.67.0 runner leaked it) before the run;
# otherwise the global environment is left empty.
#   Rscript dev/relrev-run1.R <package> <test file> [globalf]
a <- commandArgs(trailingOnly = TRUE)
if (length(a) >= 3 && a[3] == "globalf") assign("f", a[2], envir = globalenv())
local({
  .libPaths(c("C:/Users/adf44/source/r/rellib-r6",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  p <- a[1]; tf <- a[2]
  suppressMessages({library(testthat); library(p, character.only = TRUE)})
  cat("lib:", p, format(packageVersion(p)), dirname(find.package(p)),
      "; frmtmb", format(packageVersion("frmtmb")),
      dirname(find.package("frmtmb")), "\n")
  cat("globals:", paste(ls(globalenv()), collapse = ","), "\n")
  cat("start:", format(Sys.time()), "\n")
  res <- test_file(tf, package = p, env = testthat::test_env(p),
                   reporter = "silent")
  for (t in res) for (x in t$results) {
    k <- if (inherits(x, "expectation_failure")) "FAIL" else
      if (inherits(x, "expectation_error")) "ERROR" else
        if (inherits(x, "expectation_skip")) "SKIP" else
          if (inherits(x, "expectation_warning")) "WARN" else NA
    if (!is.na(k)) cat("DETAIL", k, "[", t$test, "]",
                       gsub("\n", " | ", conditionMessage(x)), "\n")
  }
  r <- as.data.frame(res)
  cat("RESULT", basename(tf), "pass=", sum(r$passed), "fail=", sum(r$failed), "err=", sum(r$error),
      "skip=", sum(r$skipped), "warn=", sum(r$warning), "\n")
})
