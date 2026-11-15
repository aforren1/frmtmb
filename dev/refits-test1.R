# One test file in one R process. Usage:
#   Rscript dev/refits-test1.R <test file> [package]
# ARM=before reads rellib-r21 alone; anything else puts the lane lib first.
args <- commandArgs(TRUE)
f <- normalizePath(args[1L], winslash = "/")
p <- if (length(args) > 1L) args[2L] else "frmtmb"
arm <- Sys.getenv("ARM", "after")
lanelib <- "C:/Users/adf44/source/r/wt-refits-lib"
base <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "before")) base else c(lanelib, base))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(p, character.only = TRUE))
cat("lib:", find.package("frmtmb"), " pkg:", find.package(p), "\n")
t0 <- proc.time()[["elapsed"]]
res <- testthat::test_file(f, package = p, env = testthat::test_env(p),
                           reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
for (i in seq_len(nrow(df))) {
  r <- res[[i]]
  for (e in r$results) {
    if (inherits(e, c("expectation_failure", "expectation_error",
                      "expectation_warning", "expectation_skip"))) {
      cat("[", class(e)[1L], "] ", r$test, ": ",
          substr(gsub("[\r\n]+", " ", conditionMessage(e)), 1L, 600L),
          "\n", sep = "")
    }
  }
}
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d tests=%d time=%.1f\n",
            basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning), nrow(df),
            proc.time()[["elapsed"]] - t0))
