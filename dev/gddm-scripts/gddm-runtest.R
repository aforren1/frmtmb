# lane gddm: run ONE test file in ONE process, with the failure cap
# lifted, and print counts that include errors.
#
# testthat's summary reporter stops at ten failures and a runner that
# sums `failed` and not `error` prints a clean line for a file that
# aborted halfway. Both are counted here.
#
# GDDM_LIB picks the arm, GDDM_TEST the file.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/gddm-lib")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
options(testthat.progress.max_fails = 1e6)
suppressMessages({
  library(testthat); library(frmtmb); library(frmtmb.eam)
})
f <- Sys.getenv("GDDM_TEST")
root <- "C:/Users/adf44/source/r/frmtmb-wt-gddm/extensions/frmtmb.eam"
path <- file.path(root, "tests", "testthat", f)
stopifnot(file.exists(path))
res <- test_file(path, package = "frmtmb.eam", reporter = "silent")
df <- as.data.frame(res)
cat(sprintf("FILE %s  PASS %d  FAIL %d  ERROR %d  WARN %d  SKIP %d\n",
            f, sum(df$passed), sum(df$failed), sum(df$error),
            sum(df$warning), sum(df$skipped)))
for (i in seq_len(nrow(df))) {
  if (df$failed[[i]] > 0L || df$error[[i]]) {
    cat("  BAD: ", df$test[[i]], "\n", sep = "")
    for (r in res[[i]][["results"]]) {
      if (inherits(r, "expectation_failure") ||
          inherits(r, "expectation_error")) {
        cat("    ", gsub("\n", " ", conditionMessage(r)), "\n", sep = "")
      }
    }
  }
}
