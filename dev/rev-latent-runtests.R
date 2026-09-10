# Reviewer: run ONE test file in THIS process against the REVIEWER's
# library, and report counts that cannot be mistaken for a clean run.
#
# Same shape as dev/latent-runtests.R, pointed at rev-latent-lib so the
# counts are taken against a package built from the worktree by this
# review rather than against whatever sits in the lane's library.
#
# `library` may be overridden so the same harness can run a MUTANT build:
#   Rscript dev/rev-latent-runtests.R <test file> [library] [test dir]

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env.R")
args <- commandArgs(trailingOnly = TRUE)
f <- args[[1L]]
if (length(args) >= 2L && nzchar(args[[2L]])) {
  .libPaths(c(args[[2L]], .libPaths()))
}
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
  library(frmtmb.latent)
})
Sys.setenv(NOT_CRAN = "true")
cat("frmtmb.latent from:", find.package("frmtmb.latent")[[1L]], "\n")

t0 <- Sys.time()
res <- testthat::test_file(f, reporter = "silent", package = "frmtmb.latent")
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
      cat("      ", gsub("\n", "\n      ", conditionMessage(r)), "\n",
          sep = "")
    }
  }
}
bad <- tot("failed") + tot("error")
quit(save = "no", status = if (bad > 0) 1L else 0L)
