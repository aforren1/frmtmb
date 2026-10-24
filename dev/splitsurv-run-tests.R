# One test file, one R process, as dev/release/run-tests.R, plus the
# per-test elapsed time, which guardrail 3 of dev/core-boundary.md
# needs and the release runner does not record.
#
# Usage: Rscript splitsurv-run-tests.R <package> <test file> <tsv out>
#
# SPLITSURV_LIB, when set, goes ahead of everything (the lane's build);
# unset, the run reads rellib-r11, the base build.
local({
LIB <- "C:/Users/adf44/source/r/rellib-r11"
.libPaths(c(if (nzchar(Sys.getenv("SPLITSURV_LIB"))) {
  Sys.getenv("SPLITSURV_LIB")
}, LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
suppressMessages(library(testthat))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
out <- a[3]
suppressMessages(library(p, character.only = TRUE))
cat("lib: ", p, " ", as.character(utils::packageVersion(p)), " ",
    dirname(find.package(p)), "; frmtmb ",
    as.character(utils::packageVersion("frmtmb")), " ",
    dirname(find.package("frmtmb")), "\n", sep = "")
t0 <- proc.time()[["elapsed"]]
res <- NULL
r <- tryCatch({
  res <- test_file(f, package = p, env = testthat::test_env(p),
                   reporter = "silent")
  as.data.frame(res)
}, error = function(e) {
  cat("RESULT ", basename(f), " LOADERROR ", conditionMessage(e), "\n",
      sep = "")
  NULL
})
el <- proc.time()[["elapsed"]] - t0
if (!is.null(res)) {
  for (t in res) {
    for (x in t$results) {
      k <- if (inherits(x, "expectation_failure")) "FAIL" else
        if (inherits(x, "expectation_error")) "ERROR" else
          if (inherits(x, "expectation_skip")) "SKIP" else
            if (inherits(x, "expectation_warning")) "WARN" else NA
      if (!is.na(k)) {
        cat("DETAIL ", k, " [", t$test, "] ",
            gsub("\n", " | ", conditionMessage(x)), "\n", sep = "")
      }
    }
  }
}
if (!is.null(r)) {
  tsv <- data.frame(package = p, file = basename(f), test = r$test,
                    nb = r$nb, passed = r$passed, failed = r$failed,
                    skipped = r$skipped, error = r$error,
                    warning = r$warning, real = r$real)
  utils::write.table(tsv, out, sep = "\t", quote = TRUE,
                     row.names = FALSE)
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      " warn=", sum(r$warning), " elapsed=", sprintf("%.2f", el), "\n",
      sep = "")
}
})
