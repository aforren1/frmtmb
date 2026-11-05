# Review of lane gpgap, item 7: one test file against one mutant.
#   REV_MUT=<mutant> Rscript dev/review-gpgap-testmut.R <file> <lib> [testdir]
a <- commandArgs(trailingOnly = TRUE)
f <- a[1]
lib <- a[2]
dir <- if (length(a) >= 3) a[3] else
  "C:/Users/adf44/source/r/rev-gpgap-src/tests/testthat"
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
res <- as.data.frame(test_file(file.path(dir, f), package = "frmtmb",
                               env = test_env("frmtmb"), reporter = "silent",
                               stop_on_failure = FALSE))
cat(sprintf("TMUT %s %s lib=%s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            Sys.getenv("REV_MUT", "none"), f, basename(dirname(find.package("frmtmb"))),
            sum(res$passed), sum(res$failed), sum(as.integer(res$error)),
            sum(as.integer(res$skipped)), sum(res$warning)))
bad <- res[res$failed > 0 | res$error, c("test", "failed", "error")]
for (i in seq_len(nrow(bad))) {
  cat(sprintf("  BAD fail=%d err=%s: %s\n", bad$failed[i], bad$error[i],
              substr(bad$test[i], 1, 70)))
}
