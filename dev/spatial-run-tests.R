# One test file, one R process, with its elapsed time: the lane's copy
# of dev/release/run-tests.R. The release runner reads rellib-r11 behind
# FRMTMB_RELLIB and prints no time; guardrail 3 of dev/core-boundary.md
# asks for each moved feature's test time and its share of the suite,
# so this one prints `secs=` on the RESULT line. Everything else is the
# release runner's logic, so the counts compare with its logs.
#
# Usage: Rscript dev/spatial-run-tests.R <package> <test file>
# Library: SPATIAL_LIBS, a ';'-separated list put ahead of the user
# library (the "before" arm passes rellib-r11 alone).
local({
libs <- strsplit(Sys.getenv("SPATIAL_LIBS"), ";", fixed = TRUE)[[1]]
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
suppressMessages(library(testthat))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
suppressMessages(library(p, character.only = TRUE))
cat("lib: ", p, " ", as.character(utils::packageVersion(p)), " ",
    dirname(find.package(p)), "; frmtmb ",
    as.character(utils::packageVersion("frmtmb")), " ",
    dirname(find.package("frmtmb")), "\n", sep = "")
res <- NULL
t0 <- proc.time()[["elapsed"]]
r <- tryCatch({
  res <- test_file(f, package = p, env = testthat::test_env(p),
                   reporter = "silent")
  as.data.frame(res)
}, error = function(e) {
  cat("RESULT ", basename(f), " LOADERROR ", conditionMessage(e), "\n",
      sep = "")
  NULL
})
secs <- proc.time()[["elapsed"]] - t0
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
  # per test_that() block, so the car/spde share of a mixed file can be
  # summed from the blocks that use them
  for (i in seq_len(nrow(r))) {
    cat("BLOCK ", basename(f), " :: ", r$test[i], " :: real=",
        format(r$real[i], digits = 4), " n=", r$nb[i], "\n", sep = "")
  }
}
if (!is.null(r)) {
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      " warn=", sum(r$warning), " secs=", format(secs, digits = 5), "\n",
      sep = "")
}
})
