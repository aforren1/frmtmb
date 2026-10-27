## One test file in one R process against the reviewer's prototype build
## (wt-mocode-rv-lib, in front of the lane library).
##   Rscript dev/review-mocode-run-test.R <package> <test file path>
args <- commandArgs(trailingOnly = TRUE)
pkg <- args[1L]; f <- args[2L]
.libPaths(c("C:/Users/adf44/source/r/wt-mocode-rv-lib",
            "C:/Users/adf44/source/r/wt-mocode-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
Sys.setenv(FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-mocode/dev/stan-cache")
suppressMessages(library(testthat))
suppressMessages(library(pkg, character.only = TRUE))
cat("lib:", dirname(system.file(package = "frmtmb")), "|", pkg,
    dirname(system.file(package = pkg)), "\n")
res <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                           reporter = testthat::ProgressReporter$new(
                             show_praise = FALSE))
df <- as.data.frame(res)
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
