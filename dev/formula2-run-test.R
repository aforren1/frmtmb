# One test file in one process. Args: <pkg> <file> [before]
# "before" drops the lane library, so the base build rellib-r3 is read.
args <- commandArgs(TRUE)
pkg <- args[1]
file <- args[2]
before <- length(args) > 2 && args[3] == "before"
lib <- "C:/Users/adf44/source/r/wt-formula2-lib"
.libPaths(c(if (!before) lib, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-formula2/dev/stan-cache")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-formula2"
dir <- if (pkg == "frmtmb") file.path(wt, "tests/testthat") else
  file.path(wt, "extensions", pkg, "tests/testthat")
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
cat("LIB", find.package(pkg), "\n")
res <- testthat::test_file(file.path(dir, file), package = pkg,
                           env = testthat::test_env(pkg),
                           reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("RESULT %s %s pass=%d fail=%d err=%d skip=%d warn=%d\n", pkg,
            file, sum(df$passed), sum(df$failed), sum(df$error),
            sum(df$skipped), sum(df$warning)))
