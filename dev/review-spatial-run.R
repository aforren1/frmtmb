# Reviewer's runner: one test file in one process.
# Usage: Rscript dev/review-spatial-run.R <libs ;-separated> <pkg> <file>
local({
a <- commandArgs(trailingOnly = TRUE)
libs <- strsplit(a[1], ";", fixed = TRUE)[[1]]
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(testthat))
p <- a[2]
f <- a[3]
suppressMessages(library(p, character.only = TRUE))
cat("lib: ", p, " ", as.character(utils::packageVersion(p)), " ",
    dirname(find.package(p)), "; frmtmb ",
    as.character(utils::packageVersion("frmtmb")), " ",
    dirname(find.package("frmtmb")), "\n", sep = "")
res <- tryCatch(test_file(f, package = p, env = testthat::test_env(p),
                          reporter = "silent"),
                error = function(e) {
                  cat("LOADERROR", conditionMessage(e), "\n")
                  NULL
                })
r <- if (is.null(res)) NULL else as.data.frame(res)
for (t in res) for (x in t$results) {
  if (inherits(x, c("expectation_failure", "expectation_error",
                    "expectation_warning"))) {
    cat("DETAIL [", t$test, "] ", class(x)[1], ": ",
        substr(gsub("\n", " | ", conditionMessage(x)), 1, 600), "\n",
        sep = "")
  }
}
if (!is.null(r)) {
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      " warn=", sum(r$warning), "\n", sep = "")
}
})
