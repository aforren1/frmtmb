# One test file in one process, against the library stack GPGAP_TLIB names
# (`;`-separated, the user library is always last).
#   Rscript dev/gpgap-runtest.R <package> <file> [<helper to source over>]
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
libs <- strsplit(Sys.getenv(
  "GPGAP_TLIB",
  "C:/Users/adf44/source/r/wt-gpgap-lib;C:/Users/adf44/source/r/rellib-r16"),
  ";")[[1]]
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(testthat)
  library(p, character.only = TRUE)
})
cat("lib:", dirname(find.package(p)), "version",
    as.character(packageVersion(p)), "\n")
# extSoftVersion() and La_library() are empty on Windows: the R copy names it
cat("blas: R home", R.home(), "| threads",
    Sys.getenv("OPENBLAS_NUM_THREADS"), "\n")
if ("frmtmb" != p) {
  cat("core lib:", dirname(find.package("frmtmb")), "\n")
}
res <- as.data.frame(test_file(f, package = p, env = test_env(p),
                               reporter = "silent", stop_on_failure = FALSE))
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), sum(res$passed), sum(res$failed),
            sum(as.integer(res$error)), sum(as.integer(res$skipped)),
            sum(res$warning)))
bad <- res[res$failed > 0 | res$error | res$warning > 0, "test"]
if (length(bad)) cat("BAD:", paste(unique(bad), collapse = " | "), "\n")
