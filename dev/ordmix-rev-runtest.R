# Reviewer of lane ordmix: run ONE test file against the lane build
# ("lane") or the base build ("base"), package attached, and print a
# RESULT line. Usage: Rscript dev/ordmix-rev-runtest.R <pkg> <file> <arm>
args <- commandArgs(TRUE)
pkg <- args[1]
file <- basename(args[2])
arm <- args[3]
libs <- c("C:/Users/adf44/source/r/rellib-r5",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) {
  libs <- c("C:/Users/adf44/source/r/wt-ordmix-lib", libs)
}
.libPaths(libs)
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordmix"
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_DRMTMB_FIT_TESTS = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(wt, "dev/ordmix-rev-stan-cache"))
options(testthat.progress.max_fails = Inf)
dir <- if (identical(pkg, "frmtmb")) file.path(wt, "tests", "testthat") else
  file.path(wt, "extensions", pkg, "tests", "testthat")
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
cat("lib:", find.package(pkg), " frmtmb:", find.package("frmtmb"), "\n")
res <- testthat::test_file(file.path(dir, file), package = pkg,
                           env = testthat::test_env(pkg),
                           reporter = testthat::ProgressReporter$new(
                             show_praise = FALSE))
df <- as.data.frame(res)
cat(sprintf("RESULT %s %s %s pass=%d fail=%d err=%d skip=%d warn=%d\n",
            arm, pkg, file, sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
