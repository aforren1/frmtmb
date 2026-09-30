# Reviewer, lane ordinal: run ONE test file of one package against the
# lane build (lane library, then rellib-r4, then the user library), with
# the package attached, and print a RESULT line.
# Usage: Rscript dev/ordinal-rev-runtest.R <pkg> <file> [lane|base]
args <- commandArgs(TRUE)
pkg <- args[1]
file <- args[2]
arm <- if (length(args) >= 3) args[3] else "lane"
libs <- c("C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) {
  libs <- c("C:/Users/adf44/source/r/wt-ordinal-lib", libs)
}
.libPaths(libs)
sp <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad"
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(sp, "ordrev-stan-cache"))
options(testthat.progress.max_fails = Inf)
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
cat("lib:", find.package(pkg), " frmtmb:", find.package("frmtmb"), "\n")
path <- if (file.exists(file)) file else {
  dir <- if (identical(pkg, "frmtmb")) file.path(wt, "tests", "testthat") else
    file.path(wt, "extensions", pkg, "tests", "testthat")
  file.path(dir, file)
}
res <- testthat::test_file(path, package = pkg,
                           env = testthat::test_env(pkg),
                           reporter = testthat::ProgressReporter$new(
                             show_praise = FALSE))
df <- as.data.frame(res)
cat(sprintf("RESULT %s %s %s pass=%d fail=%d err=%d skip=%d warn=%d\n",
            arm, pkg, basename(path), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
