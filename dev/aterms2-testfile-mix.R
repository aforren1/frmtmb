# Run one test file against this lane's library (or the base build with
# ARM=base), one file per process. Usage:
#   Rscript dev/aterms2-testfile.R <pkg> <file> [base]
args <- commandArgs(TRUE)
pkg <- args[[1]]
f <- args[[2]]
base <- length(args) > 2 && identical(args[[3]], "base")
libs <- c("C:/Users/adf44/source/r/wt-aterms2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (base) libs <- libs[-1]
if (identical(Sys.getenv("ATERMS2_MIX"), "1")) libs <- c("C:/Users/adf44/source/r/frmtmb-wt-aterms2/dev/aterms2-mixlib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(testthat))
# the summary reporter stops listing at ten failures otherwise
options(testthat.summary.max_reports = 10000L)
# attached, as run-tests.R and test_check() have it
suppressMessages(library(pkg, character.only = TRUE))
res <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                           reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("\nRESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
