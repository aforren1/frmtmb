# Same runner as dev/generics-runtests.R, for the sibling package this
# lane had to edit.
#
#   Rscript dev/generics-runtests-sample.R <one test file>
av <- commandArgs(trailingOnly = TRUE)
LIB <- "C:/Users/adf44/source/r/generics-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
# a second argument names a library to put FIRST, which is how the same
# test file is run against the base commit's frmtmb and frmtmb.sample
if (length(av) >= 2L) .libPaths(c(av[[2L]], .libPaths()))
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
  library(frmtmb.sample)
})
Sys.setenv(NOT_CRAN = "true")
message("frmtmb.sample from: ", find.package("frmtmb.sample")[[1L]],
        "  frmtmb from: ", find.package("frmtmb")[[1L]])

f <- av[[1L]]
res <- testthat::test_file(f, reporter = "silent",
                           package = "frmtmb.sample")
df <- as.data.frame(res)
tot <- function(nm) sum(df[[nm]], na.rm = TRUE)
cat("\nFILE   ", f, "\n")
cat("BLOCKS ", nrow(df), "\n")
cat("ASSERT ", tot("nb"), "\n")
cat("PASS   ", tot("nb") - tot("failed"), "\n")
cat("FAIL   ", tot("failed"), "\n")
cat("ERROR  ", tot("error"), "\n")
cat("SKIP   ", tot("skipped"), "\n")
cat("WARN   ", tot("warning"), "\n")
for (i in seq_len(nrow(df))) {
  if (isTRUE(df$failed[i] > 0) || isTRUE(df$error[i])) {
    cat("BAD    ", df$test[i], "\n")
    for (r in res[[i]]$results) {
      if (inherits(r, c("expectation_failure", "expectation_error"))) {
        cat("        ", gsub("\n", "\n         ",
                             conditionMessage(r)), "\n")
      }
    }
  }
}
if (tot("failed") > 0 || tot("error") > 0) quit(status = 1L)
