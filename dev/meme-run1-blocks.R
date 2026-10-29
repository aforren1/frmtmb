# Lane meme: like dev/meme-run1.R, and one line per test_that() block
# with its pass/fail/error counts (summary reporters stop at ten).
a <- commandArgs(TRUE)
pkg <- a[1]; f <- a[2]; arm <- if (length(a) >= 3) a[3] else "lane"
libs <- c("C:/Users/adf44/source/r/rellib-r15",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-meme-lib", libs)
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-meme/dev/stan-cache")
suppressPackageStartupMessages({library(testthat); library(pkg, character.only = TRUE)})
cat("lib:", find.package("frmtmb"), "|", find.package(pkg), "\n")
res <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                           reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
for (i in seq_len(nrow(df))) {
  cat(sprintf("BLOCK pass=%d fail=%d error=%s skip=%s | %s\n",
              df$passed[i], df$failed[i], df$error[i], df$skipped[i],
              df$test[i]))
}
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), sum(df$passed), sum(df$failed), sum(df$error),
            sum(df$skipped), sum(df$warning)))
