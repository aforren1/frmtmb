# Reviewer scratch: one test file in one process, printing a per-test
# row, so the base and lane arms of a data() change compare test by test.
#   Rscript dev/testiso-rev-onefile.R <lane|base> <file>
# base: rellib-r9 and the main checkout's tests (47cf659d);
# lane: wt-testiso-lib ahead of rellib-r9, and the worktree's tests.
local({
a <- commandArgs(TRUE)
arm <- a[1]
f <- a[2]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-testiso-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
dir <- if (arm == "lane") {
  "C:/Users/adf44/source/r/frmtmb-wt-testiso/tests/testthat"
} else "C:/Users/adf44/source/r/frmtmb/tests/testthat"
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
r <- testthat::test_file(file.path(dir, f), package = "frmtmb",
                         env = testthat::test_env("frmtmb"),
                         reporter = "silent", stop_on_failure = FALSE)
d <- as.data.frame(r)
for (i in seq_len(nrow(d))) {
  cat("ROW\t", f, "\t", d$test[i], "\t", d$passed[i], "\t", d$failed[i],
      "\t", d$skipped[i], "\t", d$error[i], "\n", sep = "")
}
for (t in r) for (x in t$results) {
  if (inherits(x, c("expectation_failure", "expectation_error",
                    "expectation_skip"))) {
    cat("NOTE\t", f, "\t", t$test, "\t", class(x)[1], "\t",
        gsub("\n", " | ", substr(conditionMessage(x), 1, 200)), "\n", sep = "")
  }
}
cat("TOTAL", arm, f, "pass", sum(d$passed), "fail", sum(d$failed), "skip",
    sum(d$skipped), "err", sum(d$error), "globals:",
    paste(ls(globalenv()), collapse = ","), "\n")
})
