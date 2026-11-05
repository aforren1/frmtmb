# One test file in one R process on this lane's library.
#   Rscript dev/refuse-runtest.R <package> <test file>
args <- commandArgs(trailingOnly = TRUE)
p <- args[1L]
f <- args[2L]
lib <- Sys.getenv("REFUSE_TEST_LIB", "C:/Users/adf44/source/r/wt-refuse-lib")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(testthat)
  library(p, character.only = TRUE)
})
cat("lib:", dirname(find.package("frmtmb")), "\n")
res <- testthat::test_file(f, package = p, env = testthat::test_env(p),
                           reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("RESULT %s %s pass=%d fail=%d error=%d skip=%d warn=%d\n", p,
            basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
for (i in which(df$failed > 0L | df$error)) {
  cat("DETAIL", df$test[i], "\n")
  for (r in res[[i]]$results) {
    if (inherits(r, c("expectation_failure", "expectation_error"))) {
      cat("  ", substr(conditionMessage(r), 1, 700), "\n")
    }
  }
}
for (i in which(df$warning > 0L)) {
  for (r in res[[i]]$results) {
    if (inherits(r, "expectation_warning")) {
      cat("WARN", df$test[i], ":", substr(conditionMessage(r), 1, 200), "\n")
    }
  }
}
for (i in which(df$skipped)) {
  for (r in res[[i]]$results) {
    if (inherits(r, "expectation_skip")) {
      cat("SKIP", df$test[i], ":", substr(conditionMessage(r), 1, 120), "\n")
    }
  }
}
