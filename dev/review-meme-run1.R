# Review of lane meme: one test file in one process.
#   Rscript dev/review-meme-run1.R <pkg> <test file> <arm: lane|base>
a <- commandArgs(TRUE)
pkg <- a[1]
f <- a[2]
arm <- a[3]
libs <- c("C:/Users/adf44/source/r/rellib-r15",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-meme-lib", libs)
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-meme/dev/stan-cache")
suppressPackageStartupMessages({
  library(testthat)
  library(pkg, character.only = TRUE)
})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")),
    "|", pkg, find.package(pkg), "\n")
t0 <- proc.time()[["elapsed"]]
res <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                           reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("RESULT %s %s pass=%d fail=%d error=%d skip=%d warn=%d secs=%.0f\n",
            arm, basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning),
            proc.time()[["elapsed"]] - t0))
bad <- df[df$failed > 0 | df$error, c("test", "failed", "error")]
if (nrow(bad)) { cat("FAILING BLOCKS:\n"); print(bad, row.names = FALSE) }
