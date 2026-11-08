# Review of lane nlstart: one test file, one process.
# Usage: REV_ARM=after|before|<lib> Rscript dev/review-nlstart-test.R <file>
# <file> is a name in tests/testthat or a path to a copy; a copy gets the
# helpers of tests/testthat sourced into its environment.
f <- commandArgs(trailingOnly = TRUE)[1]
Sys.setenv(NOT_CRAN = "true")
source("dev/review-nlstart-common.R")
suppressPackageStartupMessages(library(testthat))
path <- if (file.exists(f)) f else file.path("tests/testthat", f)
env <- testthat::test_env("frmtmb")
if (normalizePath(dirname(path)) != normalizePath("tests/testthat")) {
  for (h in list.files("tests/testthat", "^helper.*[.]R$", full.names = TRUE)) {
    sys.source(h, envir = env)
  }
}
t0 <- proc.time()[["elapsed"]]
res <- testthat::test_file(path, package = "frmtmb", env = env,
                           reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("FILE %s pass=%d fail=%d err=%d skip=%d warn=%d secs=%.0f\n", f,
            sum(df$passed), sum(df$failed), sum(df$error), sum(df$skipped),
            sum(df$warning), proc.time()[["elapsed"]] - t0))
for (i in which(df$failed > 0 | df$error | df$warning > 0)) {
  cat("  BAD:", df$test[i], "\n")
  r <- res[[i]]$results
  for (x in r) {
    if (inherits(x, c("expectation_failure", "expectation_error",
                      "expectation_warning"))) {
      cat("    -", gsub("\n", " ", substr(conditionMessage(x), 1, 300)), "\n")
    }
  }
}
