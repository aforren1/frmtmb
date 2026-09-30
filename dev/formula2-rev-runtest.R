# Reviewer: one test file in one process.
# Args: <arm before|after> <pkg> <file> <gated TRUE|FALSE>
args <- commandArgs(TRUE)
arm <- args[1]; pkg <- args[2]; file <- args[3]
gated <- identical(args[4], "TRUE")
base <- c("C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "after") c("C:/Users/adf44/source/r/wt-formula2-lib",
                                base) else base)
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = paste0(
             "C:/Users/adf44/AppData/Local/Temp/1/claude/",
             "c--Users-adf44-source-r-frmtmb/",
             "66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad/rev-stan-cache"))
if (gated) Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-formula2"
dir <- if (pkg == "frmtmb") file.path(wt, "tests/testthat") else
  file.path(wt, "extensions", pkg, "tests/testthat")
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
cat("LIB", find.package(pkg), find.package("frmtmb"), "\n")
res <- testthat::test_file(file.path(dir, file), package = pkg,
                           env = testthat::test_env(pkg),
                           reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("RESULT %s %s %s gated=%s pass=%d fail=%d err=%d skip=%d warn=%d\n",
            arm, pkg, file, gated, sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
