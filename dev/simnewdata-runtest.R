# One test file, one process, against the lane build or the base build.
#   Rscript dev/simnewdata-runtest.R <test file> [package]
#   SIMNEWDATA_LIB=base Rscript dev/simnewdata-runtest.R <test file>
# Prints every failing or erroring expectation and a RESULT line that
# counts failures AND errors, so a block that aborts is not read as a
# pass.
source("dev/simnewdata-prelude.R")
a <- commandArgs(trailingOnly = TRUE)
f <- a[1]
p <- if (length(a) >= 2L) a[2] else "frmtmb"
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
if (p != "frmtmb") suppressMessages(library(p, character.only = TRUE))
cat("library:", find.package(p), "\n")
r <- test_file(f, package = p, env = testthat::test_env(p),
               reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(r)
for (blk in r) {
  for (res in blk$results) {
    if (!inherits(res, c("expectation_success", "expectation_skip"))) {
      cat("--", class(res)[1], "in [", blk$test, "]\n  ",
          substr(gsub("\n", " ", conditionMessage(res)), 1, 400), "\n")
    }
  }
}
cat("RESULT ", basename(f), " pass=", sum(df$passed), " fail=",
    sum(df$failed), " err=", sum(df$error), " skip=", sum(df$skipped),
    " blocks=", nrow(df), "\n", sep = "")
