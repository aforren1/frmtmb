# Reviewer of lane loodraws: one test file, one R process.
#   Rscript dev/loodraws-rev-runtest.R <package> <test file> <rev|lane|base>
# rev: the reviewer's build of the lane (wt-loodraws-rev-lib); lane: the
# lane's own library; base: rellib-r10. Prints the lib: line, a DETAIL
# line per failure, error, skip and escaped warning, and a RESULT line.
# Uses the silent reporter, which writes nothing into the test directory.
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
arm <- a[3]
lib <- switch(arm, rev = "C:/Users/adf44/source/r/wt-loodraws-rev-lib",
              lane = "C:/Users/adf44/source/r/wt-loodraws-lib", base = NULL)
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (!nzchar(Sys.getenv("NOT_CRAN"))) Sys.setenv(NOT_CRAN = "true")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-loodraws/dev/stan-cache")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
suppressMessages(library(p, character.only = TRUE))
cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
    format(utils::packageVersion("frmtmb")), "|",
    dirname(system.file(package = p)), p,
    format(utils::packageVersion(p)), "\n")
res <- test_file(f, package = p, env = testthat::test_env(p),
                 reporter = "silent", stop_on_failure = FALSE)
for (t in res) {
  for (x in t$results) {
    k <- if (inherits(x, "expectation_failure")) "FAIL" else
      if (inherits(x, "expectation_error")) "ERROR" else
        if (inherits(x, "expectation_skip")) "SKIP" else
          if (inherits(x, "expectation_warning")) "WARN" else NA
    if (!is.na(k)) {
      cat("DETAIL ", k, " [", t$test, "] ",
          substr(gsub("\n", " | ", conditionMessage(x)), 1, 400), "\n",
          sep = "")
    }
  }
}
r <- as.data.frame(res)
cat("RESULT ", p, " ", basename(f), " pass=", sum(r$passed), " fail=",
    sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
    " warn=", sum(r$warning), " blocks=", nrow(r), "\n", sep = "")
