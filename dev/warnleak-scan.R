# List every warning a test file lets escape, with its location.
#
# testthat counts a warning no expectation caught as WARN in its
# summary. That count says nothing about where the warnings come from,
# and a leaked warning that was always there hides a new one next to it.
# This runs one file per process, as the release harness does, and
# writes one line per escaped warning.
#
# Usage: Rscript dev/warnleak-scan.R <package> <test file> <out file>
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]; f <- a[2]; out <- a[3]
suppressMessages({library(testthat); library(p, character.only = TRUE)})
r <- tryCatch(test_file(f, package = p, env = test_env(p),
                        reporter = "list"),
              error = function(e) NULL)
lines <- character()
for (x in r) for (e in x$results) {
  if (!inherits(e, "expectation_warning")) next
  sr <- e$srcref
  loc <- if (is.null(sr)) "?" else as.character(sr[1])
  msg <- gsub("[\r\n\t]+", " ", conditionMessage(e))
  lines <- c(lines, paste(p, basename(f), loc, x$test, substr(msg, 1, 160),
                          sep = "\t"))
}
writeLines(lines, out)
