# Lane `generics`: run ONE test file in THIS process and report counts
# that cannot be mistaken for a clean run.
#
# testthat's summary reporter caps its failure list at ten, and a runner
# that sums `failed` and not `error` prints PASS n FAIL 0 for a file
# that aborted halfway. So this reads the whole results frame and
# prints blocks, assertions, errors and skips beside the failures, and
# exits non-zero on any error or failure.
#
#   Rscript dev/generics-runtests.R <one test file> [library to put FIRST]
#
# The second argument is how the SAME test file is run against the
# unfixed build: point it at a library holding the base commit's
# frmtmb, and the test must fail.
av <- commandArgs(trailingOnly = TRUE)
LIB <- "C:/Users/adf44/source/r/generics-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (length(av) >= 2L) .libPaths(c(av[[2L]], .libPaths()))
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
})
Sys.setenv(NOT_CRAN = "true")
message("frmtmb from: ", find.package("frmtmb")[[1L]],
        "  version ", as.character(packageVersion("frmtmb")))

f <- av[[1L]]
t0 <- Sys.time()
res <- testthat::test_file(f, reporter = "silent", package = "frmtmb")
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
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
cat("SECS   ", round(secs, 1), "\n")
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
