# Lane wt-conditions: run the generated census test of every package
# source tree under a root, printing one RESULT line per package.
#   Rscript dev/conditions-census-run.R [root]
# `root` defaults to the worktree; the absent-case demonstration points
# it at a scratch copy with a bare stop() planted in it. The census
# reads sources only, so the installed build it loads (for the test
# helpers) does not matter.
av <- commandArgs(trailingOnly = TRUE)
root <- if (length(av)) av[1L] else "."
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(testthat))
dirs <- c(root, sort(list.dirs(file.path(root, "extensions"),
                               recursive = FALSE)))
for (d in dirs) {
  f <- file.path(d, "tests", "testthat", "test-conditions-census.R")
  if (!file.exists(f)) next
  pkg <- basename(normalizePath(d))
  if (!startsWith(pkg, "frmtmb.")) pkg <- "frmtmb"
  suppressMessages(library(pkg, character.only = TRUE))
  res <- test_file(f, reporter = "silent", package = pkg)
  r <- as.data.frame(res)
  cat(sprintf("RESULT %-16s pass %d fail %d error %d skip %d\n", pkg,
              sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped)))
  for (blk in res) for (x in blk$results) {
    if (inherits(x, c("expectation_failure", "expectation_error"))) {
      cat("   ", substr(gsub("\n", " / ", conditionMessage(x)), 1, 400),
          "\n")
    }
  }
}
