# Run one test file against a build, one process.
#   Rscript dev/nanse-testfile.R <lib> <pkg> <file>
a <- commandArgs(TRUE)
lib <- a[1]; pkg <- a[2]; f <- a[3]
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(pkg, character.only = TRUE))
cat("lib:", find.package(pkg), "\n")
r <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                         reporter = "summary")
d <- as.data.frame(r)
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d\n", basename(f),
            sum(d$passed), sum(d$failed), sum(d$error), sum(d$skipped),
            sum(d$warning)))
for (i in seq_len(nrow(d))) {
  if (d$failed[i] > 0 || d$error[i]) {
    cat(sprintf("  FAILED %-70s fail=%d err=%s\n", substr(d$test[i], 1, 70),
                d$failed[i], d$error[i]))
  }
}
