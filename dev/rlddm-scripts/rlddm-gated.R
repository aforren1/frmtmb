# One gated tier, in a fresh process.
#
#   Rscript dev/rlddm-scripts/rlddm-gated.R <pkg> <file> [VAR=value ...]
#
# The Stan identity tier is the one that would see a change to rlddm()'s
# density: it checks the recursion against an independent Stan program
# of the same model at the same estimates. The lane rules require
# checking the LIBRARY before believing a Stan-backed result, because
# rstan 2.32.7 cannot compile against the user library's StanHeaders
# 2.39.1 and a populated cache hides that everywhere except code that
# compiles something new. Both versions are printed below, and the
# cache is counted, so a run served entirely from cache says so.

args <- commandArgs(trailingOnly = TRUE)
pkg <- args[[1L]]
file <- args[[2L]]
for (kv in args[-(1:2)]) {
  p <- strsplit(kv, "=", fixed = TRUE)[[1L]]
  do.call(Sys.setenv, stats::setNames(list(p[[2L]]), p[[1L]]))
}
.libPaths(c("C:/Users/adf44/source/r/rlddm-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", TESTTHAT_MAX_FAILS = "10000")

cat("pkg         :", pkg, format(packageVersion(pkg)), "at",
    dirname(system.file(package = pkg)), "\n")
sh <- tryCatch(packageVersion("StanHeaders"), error = function(e) NA)
cat("StanHeaders :", format(sh), "at",
    tryCatch(dirname(system.file(package = "StanHeaders")),
             error = function(e) "-"), "\n")
cat("rstan       :",
    format(tryCatch(packageVersion("rstan"), error = function(e) NA)),
    "\n")
if (!identical(format(sh), "2.32.10")) {
  stop("StanHeaders is ", format(sh),
       " and must be 2.32.10; nothing Stan-backed below is believable",
       call. = FALSE)
}
cache <- file.path("dev", "stan-cache")
cat("stan cache  :", cache, length(list.files(cache)), "files\n")

suppressPackageStartupMessages({
  library(testthat)
  library(pkg, character.only = TRUE)
})
res <- testthat::test_file(
  file, package = pkg,
  reporter = testthat::SummaryReporter$new(max_reports = 10000L))
df <- as.data.frame(res)
cat("\nRESULT file=", basename(file), " pass=", sum(df$passed),
    " fail=", sum(df$failed), " err=", sum(df$error),
    " skip=", sum(df$skipped), " tests=", nrow(df), "\n", sep = "")
