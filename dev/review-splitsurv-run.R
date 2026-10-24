# Reviewer's runner: one test file, one R process, lane library first.
# Usage: Rscript review-splitsurv-run.R <package> <test file>
# REV_LIB overrides the first library (default: the lane's).
local({
first <- Sys.getenv("REV_LIB", "C:/Users/adf44/source/r/wt-splitsurv-lib")
.libPaths(c(if (nzchar(first) && first != "none") first,
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(testthat))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]; f <- a[2]
suppressMessages(library(p, character.only = TRUE))
cat("lib: ", p, " ", as.character(utils::packageVersion(p)), " ",
    dirname(find.package(p)), "; frmtmb ",
    as.character(utils::packageVersion("frmtmb")), " ",
    dirname(find.package("frmtmb")), "\n", sep = "")
t0 <- proc.time()[["elapsed"]]
res <- test_file(f, package = p, env = testthat::test_env(p),
                 reporter = "silent")
r <- as.data.frame(res)
for (t in res) for (x in t$results) {
  k <- if (inherits(x, "expectation_failure")) "FAIL" else
    if (inherits(x, "expectation_error")) "ERROR" else
      if (inherits(x, "expectation_skip")) "SKIP" else
        if (inherits(x, "expectation_warning")) "WARN" else NA
  if (!is.na(k)) cat("DETAIL ", k, " [", t$test, "] ",
                     gsub("\n", " | ", conditionMessage(x)), "\n", sep = "")
}
cat("RESULT ", basename(f), " tests=", nrow(r), " pass=", sum(r$passed),
    " fail=", sum(r$failed), " err=", sum(r$error), " skip=",
    sum(r$skipped), " warn=", sum(r$warning), " elapsed=",
    sprintf("%.1f", proc.time()[["elapsed"]] - t0), "\n", sep = "")
})
