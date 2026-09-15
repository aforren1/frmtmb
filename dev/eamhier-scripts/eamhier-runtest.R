# Lane eamhier: one test file per process.
#
# Counts errors as well as failures, and prints the skip count, because
# a runner that sums `failed` alone prints a clean line for a file that
# aborted halfway.
#
# Run:  Rscript --vanilla dev/eamhier-scripts/eamhier-runtest.R \
#         <package dir> <test file>
source("dev/eamhier-scripts/eamhier-common.R")
eamhier_libs()
suppressMessages({
  library(testthat)
  library(frmtmb)
})
args <- commandArgs(trailingOnly = TRUE)
pkgdir <- args[[1L]]
f <- args[[2L]]
pkg <- basename(pkgdir)
suppressMessages(library(pkg, character.only = TRUE))
setwd(file.path(pkgdir, "tests", "testthat"))
r <- as.data.frame(test_file(f, reporter = "silent", package = pkg))
cat(sprintf("FILE %s tests %d pass %d fail %d error %d skip %d\n",
            f, nrow(r), sum(r$passed), sum(r$failed), sum(r$error),
            sum(r$skipped)))
bad <- r[r$failed > 0 | r$error, , drop = FALSE]
if (nrow(bad)) {
  for (i in seq_len(nrow(bad))) {
    cat("  BAD:", bad$test[i], "| failed", bad$failed[i], "| error",
        bad$error[i], "\n")
  }
}
