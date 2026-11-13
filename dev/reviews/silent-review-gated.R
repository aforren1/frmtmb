# One test file in one R process, against the base build or the lane's.
# Usage: Rscript dev/silent-run-test.R base|lane <package> <test file>
args <- commandArgs(TRUE)
arm <- args[1]
pkg <- args[2]
f <- args[3]
libs <- c("C:/Users/adf44/source/r/rellib-r20",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) {
  libs <- c("C:/Users/adf44/source/r/wt-silent-lib", libs)
}
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true", FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-silent/dev/stan-cache", R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(testthat)
  library(pkg, character.only = TRUE)
})
cat("arm:", arm, " file:", f, "\n")
cat("lib:", find.package(pkg), " frmtmb:", find.package("frmtmb"), "\n")
t0 <- proc.time()[["elapsed"]]
res <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                           reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
for (k in seq_len(nrow(df))) {
  cat(sprintf("TEST %-70s nb=%d failed=%d skipped=%s error=%s warning=%d\n",
              substr(df$test[k], 1, 70), df$nb[k], df$failed[k],
              df$skipped[k], df$error[k], df$warning[k]))
}
cnt <- c(success = 0L, failure = 0L, error = 0L, skip = 0L, warning = 0L)
for (t in res) {
  for (e in t$results) {
    k <- sub("^expectation_", "", class(e)[1])
    if (k %in% names(cnt)) cnt[[k]] <- cnt[[k]] + 1L
    if (!inherits(e, "expectation_success")) {
      cat("---", class(e)[1], "in", t$test, "\n")
      cat(substr(conditionMessage(e), 1, 1500), "\n")
    }
  }
}
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d tests=%d secs=%.1f\n",
            basename(f), cnt[["success"]], cnt[["failure"]], cnt[["error"]],
            cnt[["skip"]], cnt[["warning"]], nrow(df),
            proc.time()[["elapsed"]] - t0))
