# Reviewer's own runner: one test file, one package, a given library
# stack, counts from the result object (not from the lane's runner).
#
#   Rscript dev/review-ci-runfile.R <test file> <pkg> <lib1;lib2;...>
args <- commandArgs(trailingOnly = TRUE)
f <- normalizePath(args[[1L]], winslash = "/")
p <- args[[2L]]
.libPaths(c(strsplit(args[[3L]], ";", fixed = TRUE)[[1L]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages(library(p, character.only = TRUE))
cat("lib:", p, as.character(packageVersion(p)), find.package(p), "\n")
cat("file:", f, "\n")
t0 <- proc.time()[["elapsed"]]
res <- testthat::test_file(f, package = p, env = testthat::test_env(p),
                           reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warning=%d tests=%d secs=%.0f\n",
            basename(f), sum(df$passed), sum(df$failed),
            sum(as.integer(df$error)), sum(as.integer(df$skipped)),
            sum(df$warning), nrow(df), proc.time()[["elapsed"]] - t0))
bad <- df[df$failed > 0 | df$error | df$warning > 0, c("test", "failed",
                                                       "error", "warning")]
if (nrow(bad)) print(bad)
