# Lane wt-defects: one test file, one process.
#   DEFECTS_ARM=before|after Rscript dev/defects-runtest.R <pkg> <file>
# "before" reads rellib-r3 (base 0.65.0); "after" puts the lane library
# first. Every failing or erroring expectation's message is printed, so
# a log shows HOW a test failed on the base build, not only that it did.
a <- commandArgs(trailingOnly = TRUE)
arm <- Sys.getenv("DEFECTS_ARM", "after")
libs <- c("C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (!identical(arm, "before")) {
  libs <- c("C:/Users/adf44/source/r/wt-defects-lib", libs)
}
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true")
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (file.exists(mk)) Sys.setenv(R_MAKEVARS_USER = mk)
suppressMessages(library(testthat))
p <- a[1]
suppressMessages(library(p, character.only = TRUE))
res <- test_file(a[2], package = p, env = test_env(p), reporter = "silent",
                 stop_on_failure = FALSE)
r <- as.data.frame(res)
cat("RESULT", arm, p, format(packageVersion(p)), basename(a[2]),
    "pass=", sum(r$passed), "fail=", sum(r$failed), "err=",
    sum(r$error), "skip=", sum(r$skipped), "warn=", sum(r$warning), "\n")
for (t in res) {
  for (e in t$results) {
    if (inherits(e, c("expectation_failure", "expectation_error",
                      "expectation_warning"))) {
      cat("--", class(e)[1], "in:", t$test, "\n  ",
          gsub("\n", "\n   ", substr(conditionMessage(e), 1, 600)), "\n")
    }
  }
}
