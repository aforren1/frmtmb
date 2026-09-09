# One test file, one R process, with a count that cannot be capped or
# silently truncated: every expectation is classified, and errors are
# counted separately from failures.
#
# Usage: Rscript rev-lincmt-testfile.R <lib> <pkgdir> <file>
args <- commandArgs(trailingOnly = TRUE)
LIB <- args[[1L]]
PKG <- args[[2L]]
FILE <- args[[3L]]
.libPaths(c(LIB, "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(testthat); library(frmtmb); library(frmtmb.ode)
})
cat("file:", FILE, "\n")
cat("frmtmb.ode from:", dirname(getNamespaceInfo("frmtmb.ode", "path")),
    "\n")
p <- file.path(PKG, "tests", "testthat", FILE)
res <- test_file(p, reporter = "silent", package = NULL)
df <- as.data.frame(res)
cnt <- colSums(df[, c("failed", "skipped", "error", "warning",
                      "passed"), drop = FALSE] *
                 rep(1, nrow(df)))
cat(sprintf("PASS %d  FAIL %d  ERROR %d  SKIP %d  WARN %d\n",
            sum(df$passed), sum(df$failed), sum(df$error),
            sum(df$skipped), sum(df$warning)))
for (i in seq_len(nrow(df))) {
  if (df$failed[[i]] > 0 || df$error[[i]] > 0) {
    cat("  -- ", df$test[[i]], " failed", df$failed[[i]], " error",
        df$error[[i]], "\n")
    for (r in res[[i]]$results) {
      if (inherits(r, "expectation_failure") ||
            inherits(r, "expectation_error")) {
        cat("     ", substr(paste(conditionMessage(r), collapse = " "),
                            1, 220), "\n")
      }
    }
  }
}
