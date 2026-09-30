# Run ONE test file of one package against the lane build (default) or
# the base build (arm "base"), with the package attached, and print a
# RESULT line with pass, fail, error, skip and warning counts.
# Usage: Rscript dev/ordinal-runtest.R <pkg> <file> [lane|base|mut:<lib>]
#        [gated]
# "mut:<lib>" puts a mutant core build ahead of the lane library.
args <- commandArgs(TRUE)
pkg <- args[1]
file <- args[2]
arm <- if (length(args) >= 3) args[3] else "lane"
libs <- c("C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) {
  libs <- c("C:/Users/adf44/source/r/wt-ordinal-lib", libs)
}
if (startsWith(arm, "mut:")) {
  libs <- c(sub("^mut:", "", arm), "C:/Users/adf44/source/r/wt-ordinal-lib",
            libs)
}
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-ordinal/dev/stan-cache")
if (length(args) >= 4 && identical(args[4], "gated")) {
  Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")
}
# the progress reporter stops listing after ten failures unless told not
# to, and a capped count is not a count
options(testthat.progress.max_fails = Inf)
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
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
