# Reviewer (round 2) copy of dev/gpby-runtest.R (own Stan cache copy). Run ONE test file against a library, one process per file.
# Usage: Rscript dev/gpby-runtest.R <pkg> <file> [base]
# Prints a RESULT line with pass/fail/error/skip/warn counts.
args <- commandArgs(TRUE)
pkg <- args[1]
file <- args[2]
base <- length(args) > 2 && args[3] == "base"
LIB <- "C:/Users/adf44/source/r/wt-gpby-lib"
.libPaths(c(if (!base) LIB, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-gpby/dev/gpby-rev2-stan-cache")
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
cat("lib:", find.package("frmtmb"), "|", find.package(pkg), "\n")
res <- testthat::test_file(file, package = pkg,
                           env = testthat::test_env(pkg),
                           reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
cnt <- c(pass = sum(df$passed), fail = sum(df$failed),
         error = sum(df$error), skip = sum(df$skipped),
         warn = sum(df$warning))
for (i in seq_len(nrow(df))) {
  if (df$failed[i] > 0 || df$error[i]) {
    cat("BAD:", df$test[i], "\n")
  }
}
for (r in res) {
  for (e in r$results) {
    if (inherits(e, c("expectation_failure", "expectation_error"))) {
      cat("----", r$test, "\n", conditionMessage(e), "\n")
    }
    if (inherits(e, "expectation_warning")) {
      cat("WARN", r$test, ":", conditionMessage(e), "\n")
    }
  }
}
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(file), cnt[["pass"]], cnt[["fail"]], cnt[["error"]],
            cnt[["skip"]], cnt[["warn"]]))
