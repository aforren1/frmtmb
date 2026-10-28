## Reviewer: one test file in one R process.
##   Rscript dev/review-mointerop-run-test.R <lane|base> <package> <file>
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1L]; pkg <- args[2L]; f <- args[3L]
base <- c("C:/Users/adf44/source/r/rellib-r14",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mointerop-lib", base)
} else base)
Sys.setenv(NOT_CRAN = "true")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE =
               "C:/Users/adf44/source/r/frmtmb-wt-mointerop/dev/stan-cache")
}
suppressMessages(library(testthat))
suppressMessages(library(pkg, character.only = TRUE))
cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "|", pkg,
    dirname(system.file(package = pkg)), "\n")
res <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                           reporter = testthat::ProgressReporter$new(
                             show_praise = FALSE))
df <- as.data.frame(res)
cat(sprintf("RESULT %s %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            arm, basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
