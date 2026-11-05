# One test file, one process, on the build a command-line library picks:
#   Rscript dev/impsamp-testfile.R <lib or ""> <package> <path to test file>
args <- commandArgs(TRUE)
lib <- args[[1]]
pkg <- args[[2]]
f <- args[[3]]
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "|", find.package(pkg), "\n")
Sys.setenv(NOT_CRAN = "true")
t0 <- proc.time()[["elapsed"]]
res <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                           reporter = "silent")
el <- proc.time()[["elapsed"]] - t0
df <- as.data.frame(res)
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d secs=%.1f\n",
            basename(f), sum(df$nb) - sum(df$failed) - sum(df$skipped),
            sum(df$failed), sum(df$error), sum(df$skipped), sum(df$warning),
            el))
bad <- df[df$failed > 0 | df$error | df$warning > 0, ]
for (i in seq_len(nrow(bad))) {
  cat("BAD:", bad$test[i], "| fail", bad$failed[i], "err", bad$error[i],
      "warn", bad$warning[i], "\n")
  for (r in res[[which(df$test == bad$test[i])[1]]]$results) {
    if (!inherits(r, "expectation_success")) {
      cat("   ", substr(conditionMessage(r), 1, 600), "\n")
    }
  }
}
print(df[order(-df$real), c("test", "real")][1:min(8, nrow(df)), ],
      row.names = FALSE)
if (length(args) >= 4L) {
  utils::write.table(df[, c("test", "nb", "failed", "skipped", "error", "real")],
                     args[[4L]], sep = "\t", row.names = FALSE)
}
