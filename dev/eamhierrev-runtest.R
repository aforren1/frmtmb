# REVIEW of lane eamhier: one test file per process, in a directory
# given on the command line, so that a MUTATED copy of a test file can
# be run without touching the lane's tree.
#
# Run: Rscript --vanilla dev/eamhierrev-runtest.R <dir> <file> [<pkg>]

.libPaths(c("C:/Users/adf44/source/r/eamhierrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(testthat)
  library(frmtmb)
})
args <- commandArgs(trailingOnly = TRUE)
dir <- args[[1L]]
f <- args[[2L]]
pkg <- if (length(args) >= 3L) args[[3L]] else "frmtmb.eam"
suppressMessages(library(pkg, character.only = TRUE))
setwd(dir)
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
