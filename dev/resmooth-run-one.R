# Lane wt-resmooth. One test file in one R process, counts to stdout.
#   Rscript dev/resmooth-run-one.R <file> [base]
# `base` loads the 0.64.0 reference build instead of the lane library,
# which is how a new test is SEEN TO FAIL before the fix.
args <- commandArgs(trailingOnly = TRUE)
f <- args[[1L]]
base <- length(args) > 1L && identical(args[[2L]], "base")
.libPaths(c(if (!base) "C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressMessages(library(frmtmb))
pkg <- if (grepl("extensions/", f, fixed = TRUE)) {
  sub(".*extensions/([^/]+)/.*", "\\1", f)
} else "frmtmb"
if (pkg != "frmtmb") suppressMessages(library(pkg, character.only = TRUE))
res <- testthat::test_file(f, package = pkg,
                           env = testthat::test_env(pkg),
                           reporter = "silent")
df <- as.data.frame(res)
tot <- function(nm) sum(df[[nm]])
cat(sprintf("RESULT %s lib=%s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), if (base) "base" else "lane", tot("passed"),
            tot("failed"), tot("error"), tot("skipped"), tot("warning")))
bad <- df[df$failed > 0 | df$error > 0, , drop = FALSE]
if (nrow(bad)) {
  for (i in seq_len(nrow(bad))) {
    cat("  BAD:", bad$test[[i]], "\n")
    for (r in res[[which(vapply(res, function(x) identical(x$test,
                                                           bad$test[[i]]),
                                NA))[1L]]]$results) {
      if (inherits(r, c("expectation_failure", "expectation_error"))) {
        cat("    ", gsub("\n", " | ", substr(conditionMessage(r), 1, 400)),
            "\n")
      }
    }
  }
}
