# Count the warnings one test file emits, for the stray-warning item.
#   Rscript dev/priorform-warncount.R <testdir> <filter> ref|lane
# <testdir> may be the base checkout's tests/testthat (read, never
# written): testthat runs in a scratch COPY so nothing lands there.
args <- commandArgs(trailingOnly = TRUE)
src <- args[1]
filt <- args[2]
lib <- switch(args[3],
  ref = "C:/Users/adf44/source/r/rellib-r3",
  lane = "C:/Users/adf44/source/r/priorform-lib")
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
tmp <- file.path(tempdir(), "priorform-warncount")
dir.create(tmp, showWarnings = FALSE)
keep <- list.files(src, pattern = paste0("^(helper.*|setup.*|test-", filt,
                                         ")[.]R$"), full.names = TRUE)
invisible(file.copy(keep, tmp, overwrite = TRUE))
res <- testthat::test_dir(tmp, filter = filt, package = "frmtmb",
                          reporter = "silent", stop_on_failure = FALSE)
d <- as.data.frame(res)
cat(sprintf("PRIORFORM-WARN %s [%s, %s] pass %d fail %d warning %d\n",
            filt, src, args[3], sum(d$passed), sum(d$failed),
            sum(d$warning)))
w <- unlist(lapply(res, function(r) {
  vapply(Filter(function(x) inherits(x, "expectation_warning"), r$results),
         conditionMessage, "")
}))
if (length(w)) writeLines(paste("  warning:", unique(w)))
