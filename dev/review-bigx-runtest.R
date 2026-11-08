# usage: Rscript review-bigx-runtest.R <frmtmb lib or "before"> <test file>
a <- commandArgs(TRUE)
lib <- a[1]; tf <- a[2]
.libPaths(c(if (lib != "before") lib, "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressPackageStartupMessages(library(frmtmb))
suppressPackageStartupMessages(library(testthat))
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
cat("file:", tf, "\n")
p <- "frmtmb"
res <- testthat::test_file(tf, package = p, env = testthat::test_env(p),
                           reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("pass=%d fail=%d err=%d skip=%d warn=%d\n", sum(df$passed),
            sum(df$failed), sum(df$error), sum(df$skipped), sum(df$warning)))
bad <- df[df$failed > 0 | df$error | df$warning > 0, "test"]
if (length(bad)) cat("failing tests:\n", paste0("  ", bad, "\n"))
for (r in res) for (e in r$results) {
  if (inherits(e, c("expectation_failure", "expectation_error",
                    "expectation_warning"))) {
    cat("--", r$test, ":", substr(conditionMessage(e), 1, 300), "\n")
  }
}
