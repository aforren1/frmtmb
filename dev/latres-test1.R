# One test file against the lane's library, with counts that no reporter caps.
# Usage: Rscript dev/latres-test1.R <file> [package]
args <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
f <- args[[1]]
p <- if (length(args) > 1L) args[[2]] else "frmtmb"
library(testthat)
library(p, character.only = TRUE)
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
res <- test_file(f, package = p, env = test_env(p), reporter = "silent",
                 stop_on_failure = FALSE)
df <- as.data.frame(res)
for (r in res) {
  for (e in r$results) {
    if (inherits(e, "expectation_success")) next
    cl <- if (inherits(e, "expectation_failure")) "FAIL" else
      if (inherits(e, "expectation_error")) "ERROR" else
        if (inherits(e, "expectation_skip")) "SKIP" else
          if (inherits(e, "expectation_warning")) "WARN" else class(e)[1L]
    cat(sprintf("%s [%s] %s\n  %s\n", cl, r$file, r$test,
                substr(gsub("\n", " | ", conditionMessage(e)), 1L, 600L)))
  }
}
cat(sprintf("RESULT %s %s: tests=%d pass=%d fail=%d err=%d skip=%d warn=%d time=%.1f\n",
            p, basename(f), sum(df$nb), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning), sum(df$real)))
