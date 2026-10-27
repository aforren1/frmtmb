## One test file in one R process, against the lane library or the base.
##   Rscript dev/mocode-run-test.R <lane|base> <package> <test file path>
## Prints a RESULT line with pass, fail, error, skip and warning counts.
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1L]; pkg <- args[2L]; f <- args[3L]
base <- c("C:/Users/adf44/source/r/rellib-r13",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mocode-lib", base)
} else base)
Sys.setenv(NOT_CRAN = "true")
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE =
               "C:/Users/adf44/source/r/frmtmb-wt-mocode/dev/stan-cache")
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
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
