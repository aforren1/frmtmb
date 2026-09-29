# Reviewer's own one-file-one-process runner.
#   Rscript dev/resmooth-rev-run-one.R <file> [base]
# Prints pass/fail/error/skip/warn AND the number of test blocks the file
# reported against the number of `test_that(` calls its source holds, so
# a file that aborted part way cannot print a clean line.
args <- commandArgs(trailingOnly = TRUE)
f <- args[[1L]]
base <- length(args) > 1L && identical(args[[2L]], "base")
.libPaths(c(if (!base) "C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
# REV_GATE=off leaves the brms fit tier gated, for the two files that
# compile a fresh Stan program per model and cost hours
if (!identical(Sys.getenv("REV_GATE"), "off")) {
  Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")
}
suppressMessages(library(frmtmb))
pkg <- if (grepl("extensions/", f, fixed = TRUE)) {
  sub(".*extensions/([^/]+)/.*", "\\1", f)
} else "frmtmb"
if (pkg != "frmtmb") suppressMessages(library(pkg, character.only = TRUE))
src <- readLines(f, warn = FALSE)
declared <- sum(grepl("^\\s*test_that\\(", src))
res <- tryCatch(testthat::test_file(f, package = pkg,
                                    env = testthat::test_env(pkg),
                                    reporter = "silent"),
                error = function(e) e)
if (inherits(res, "error")) {
  cat(sprintf("RESULT %s lib=%s ABORTED: %s\n", basename(f),
              if (base) "base" else "lane",
              gsub("\n", " | ", conditionMessage(res))))
  quit(save = "no", status = 0)
}
df <- as.data.frame(res)
tot <- function(nm) sum(as.numeric(df[[nm]]))
cat(sprintf(paste0("RESULT %s lib=%s pass=%d fail=%d error=%d skip=%d ",
                   "warn=%d blocks=%d/%d\n"),
            basename(f), if (base) "base" else "lane", tot("passed"),
            tot("failed"), tot("error"), tot("skipped"), tot("warning"),
            nrow(df), declared))
bad <- df[df$failed > 0 | as.numeric(df$error) > 0, , drop = FALSE]
if (nrow(bad)) {
  for (i in seq_len(nrow(bad))) cat("  BAD:", bad$test[[i]], "\n")
}
sk <- df[as.numeric(df$skipped) > 0, , drop = FALSE]
if (nrow(sk)) for (i in seq_len(nrow(sk))) cat("  SKIP:", sk$test[[i]], "\n")
