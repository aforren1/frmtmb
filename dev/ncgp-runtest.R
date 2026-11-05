# One test file in one process, against the library stack NCGP_TLIB names
# (`;`-separated, the user library is always last). Prints each failure's
# message, cut short, under its test.
#   Rscript dev/ncgp-runtest.R <package> <file>
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
libs <- strsplit(Sys.getenv(
  "NCGP_TLIB",
  "C:/Users/adf44/source/r/wt-ncgp-lib;C:/Users/adf44/source/r/rellib-r16"),
  ";")[[1]]
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(testthat)
  library(p, character.only = TRUE)
})
cat("lib:", dirname(find.package(p)), "version",
    as.character(packageVersion(p)), "\n")
cat("blas:", extSoftVersion()[["BLAS"]], "lapack:", La_library(), "\n")
if ("frmtmb" != p) {
  cat("core lib:", dirname(find.package("frmtmb")), "\n")
}
r <- test_file(f, package = p, env = test_env(p), reporter = "silent",
               stop_on_failure = FALSE)
res <- as.data.frame(r)
for (tt in r) {
  for (e in tt$results) {
    if (inherits(e, c("expectation_failure", "expectation_error"))) {
      cat("FAILMSG [", tt$test, "]", substr(conditionMessage(e), 1, 300), "\n")
    } else if (inherits(e, "expectation_warning")) {
      cat("WARNMSG [", tt$test, "]", substr(conditionMessage(e), 1, 200), "\n")
    }
  }
}
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), sum(res$passed), sum(res$failed),
            sum(as.integer(res$error)), sum(as.integer(res$skipped)),
            sum(res$warning)))
bad <- res[res$failed > 0 | res$error | res$warning > 0, "test"]
if (length(bad)) cat("BAD:", paste(unique(bad), collapse = " | "), "\n")
