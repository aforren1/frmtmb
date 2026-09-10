# Lane `latent`: run ONE test file in THIS process and report counts
# that cannot be mistaken for a clean run.
#
# The summary reporter caps its failure list at ten, and a runner that
# sums `failed` and not `error` prints PASS n FAIL 0 for a file that
# aborted halfway. So this reads the whole results frame, prints the
# block, assertion and SKIP counts before the failure count, and exits
# non-zero on any error or failure.
#
#   Rscript dev/latent-runtests.R <path to one test file>

source("dev/latent-env.R")
av <- commandArgs(trailingOnly = TRUE)
# a SECOND argument names a library to put first, which is how a
# deliberately broken build is run without touching this lane's own
if (length(av) >= 2L) .libPaths(c(av[[2L]], .libPaths()))
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
  library(frmtmb.latent)
})
Sys.setenv(NOT_CRAN = "true")
message("frmtmb.latent from: ", find.package("frmtmb.latent")[[1L]])

f <- av[[1L]]
t0 <- Sys.time()
res <- testthat::test_file(f, reporter = "silent",
                           package = "frmtmb.latent")
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
cat("SECONDS", format(secs, digits = 4), "\n")

for (i in seq_len(nrow(df))) {
  flag <- if (isTRUE(df$error[i])) "ERROR" else
    if (isTRUE(df$failed[i] > 0)) "FAIL" else
      if (isTRUE(df$skipped[i])) "SKIP" else
        if (isTRUE(df$warning[i] > 0)) "WARN" else NA_character_
  if (is.na(flag)) next
  x <- res[[i]]
  cat("  [", flag, "] ", df$test[i], "\n", sep = "")
  for (r in x$results) {
    if (inherits(r, c("expectation_failure", "expectation_error",
                      "expectation_skip", "expectation_warning"))) {
      cat("      ", gsub("\n", "\n      ", conditionMessage(r)),
          "\n", sep = "")
    }
  }
}
bad <- tot("failed") + tot("error")
quit(save = "no", status = if (bad > 0) 1L else 0L)
