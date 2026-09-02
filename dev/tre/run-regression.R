# Regression runner for the wt-tre lane. Writes to dev/tre/regression.txt.
# Run: Rscript dev/tre/run-regression.R [file ...]
args <- commandArgs(trailingOnly = TRUE)
fs <- if (length(args)) args else c(
  "test-tre", "test-covstruct", "test-parse", "test-frame",
  "test-methods", "test-methods-audit", "test-simulate-ergonomics",
  "test-priors-bounds-grcov", "test-id-kron", "test-data2",
  "test-compat", "test-review-v29")
sink("dev/tre/regression.txt")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
library(testthat)
tot <- data.frame()
for (f in fs) {
  p <- file.path("tests/testthat", paste0(f, ".R"))
  if (!file.exists(p)) { cat("MISSING ", f, "\n"); next }
  cat("\n#### ", f, "\n")
  r <- as.data.frame(test_file(p, reporter = "silent"))
  cat(sprintf("  pass %d  fail %d  error %d  warn %d  skip %d\n",
              sum(r$passed), sum(r$failed), sum(r$error),
              sum(r$warning), sum(r$skipped)))
  if (sum(r$failed) || sum(r$error)) {
    bad <- r[r$failed > 0 | r$error, ]
    for (i in seq_len(nrow(bad))) {
      cat("  X ", bad$test[i], "\n")
      for (res in bad$result[[i]]) {
        if (inherits(res, c("expectation_failure", "expectation_error"))) {
          cat("      ", gsub("\n", "\n      ",
                             conditionMessage(res)), "\n")
        }
      }
    }
  }
  tot <- rbind(tot, data.frame(file = f, pass = sum(r$passed),
                               fail = sum(r$failed), err = sum(r$error),
                               skip = sum(r$skipped)))
}
cat("\n==== totals ====\n")
print(tot, row.names = FALSE)
cat(sprintf("\nTOTAL pass %d fail %d error %d skip %d\n", sum(tot$pass),
            sum(tot$fail), sum(tot$err), sum(tot$skip)))
sink()
cat("done\n")
