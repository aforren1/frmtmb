# Run ONE test file in this process against a named library, counting
# failures AND errors (a runner that sums only `failed` prints a clean
# line for a file that aborted).
a <- commandArgs(trailingOnly = TRUE)
LIB <- a[1]; FILE <- a[2]
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
library(testthat)
library(frmtmb)
if (file.exists(file.path(dirname(FILE), "setup.R"))) {
  sys.source(file.path(dirname(FILE), "setup.R"), envir = globalenv())
}
res <- test_file(FILE, reporter = "silent", package = "frmtmb")
d <- as.data.frame(res)
np <- sum(d$passed); nf <- sum(d$failed); ne <- sum(d$error)
ns <- sum(d$skipped)
cat(sprintf("FILE %s LIB %s\n", basename(FILE), LIB))
cat(sprintf("BLOCKS %d ASSERT %d PASS %d FAIL %d ERROR %d SKIP %d\n",
            nrow(d), np + nf + ne, np, nf, ne, ns))
for (i in seq_len(nrow(d))) {
  if (d$failed[i] > 0 || d$error[i] > 0) {
    cat(sprintf("  FAILBLOCK: %s (failed=%d error=%d)\n", d$test[i],
                d$failed[i], d$error[i]))
  }
}
cat("GENREVDONE\n")
