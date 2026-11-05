# Rscript dev/nlfix-testrun.R <pkg> <test-file> <log>: one test file, one R
# process, on the lane's library (NLFIX_LIB overrides, for a "before" arm).
args <- commandArgs(trailingOnly = TRUE)
p <- args[1]
f <- args[2]
log <- args[3]
lib <- Sys.getenv("NLFIX_LIB",
                  "C:/Users/adf44/source/r/wt-nlfix-lib2;C:/Users/adf44/source/r/rellib-r16")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                               winslash = "/"))
}
library(p, character.only = TRUE)
suppressPackageStartupMessages(library(testthat))
det <- testthat::test_file(f, package = p, env = testthat::test_env(p),
                           reporter = testthat::ListReporter$new(),
                           stop_on_failure = FALSE)
res <- as.data.frame(det)
out <- c(sprintf("lib: %s", find.package(p)),
         sprintf("core lib: %s", find.package("frmtmb")),
         sprintf("FILE %s PASS=%d FAIL=%d ERROR=%d SKIP=%d WARN=%d", basename(f),
                 sum(res$passed), sum(res$failed), sum(res$error), sum(res$skipped),
                 sum(res$warning)))
sk <- res[res$skipped, "test"]
if (length(sk)) out <- c(out, paste("skipped:", sk))
msgs <- unlist(lapply(det, function(r) {
  lapply(r$results, function(x) {
    if (inherits(x, c("expectation_failure", "expectation_error",
                      "expectation_warning"))) {
      paste0("[", r$test, "] ", class(x)[1], ": ",
             substr(conditionMessage(x), 1, 500))
    }
  })
}))
if (length(msgs)) out <- c(out, "details:", msgs)
writeLines(out, log)
cat(out[3], "\n")
