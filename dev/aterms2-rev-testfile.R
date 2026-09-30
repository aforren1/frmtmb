# Reviewer's runner: one test file per process, the package attached.
#   Rscript dev/aterms2-rev-testfile.R <pkg> <abs test file> <lane|base>
# The test file's own directory supplies the helpers.
args <- commandArgs(TRUE)
pkg <- args[[1]]
f <- args[[2]]
arm <- args[[3]]
libs <- c("C:/Users/adf44/source/r/wt-aterms2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressPackageStartupMessages(library(testthat))
options(testthat.summary.max_reports = 10000L)
suppressMessages(library(pkg, character.only = TRUE))
cat("ARM", arm, "PKG", find.package(pkg), "\n")
res <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                           reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("\nRESULT %s %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            arm, basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
