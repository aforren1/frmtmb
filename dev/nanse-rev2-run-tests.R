# Reviewer, punch round 1: one test file, one R process, against the
# trial merge (nanse-rev-lib, then rellib-r6). Unlike
# dev/nanse-rev-run-tests.R it defines no global variable (a global `f`
# answered a test's lookup there) and does not trace frm(); it records
# each SE warning (se_lost_message()) to a .fire file.
#   Rscript dev/nanse-rev2-run-tests.R <package> <test file> <fire log>
local({
  lib <- "C:/Users/adf44/source/r/nanse-rev-lib"
  .libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r6",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
  if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
    Sys.setenv(R_MAKEVARS_USER = mk)
  }
})
suppressMessages(library(testthat))
.rev2_args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(.rev2_args[1], character.only = TRUE))
cat("lib: ", find.package("frmtmb"), " | ", find.package(.rev2_args[1]),
    "\n", sep = "")
Sys.setenv(NANSE_REV2_FIRE = .rev2_args[3])
suppressMessages(trace(
  "se_lost_message", where = asNamespace("frmtmb"), print = FALSE,
  exit = quote(cat(gsub("[\r\n\t]+", " ", returnValue()), "\n",
                   file = Sys.getenv("NANSE_REV2_FIRE"), append = TRUE))))
.rev2_r <- tryCatch(
  as.data.frame(test_file(.rev2_args[2], package = .rev2_args[1],
                          env = testthat::test_env(.rev2_args[1]),
                          reporter = "silent")),
  error = function(e) {
    cat("RESULT ", basename(.rev2_args[2]), " LOADERROR ",
        conditionMessage(e), "\n", sep = "")
    NULL
  })
if (!is.null(.rev2_r)) {
  .rev2_t <- .rev2_r
  cat("RESULT ", basename(.rev2_args[2]), " pass=", sum(.rev2_t$passed), " fail=",
      sum(.rev2_t$failed), " err=", sum(.rev2_t$error), " skip=", sum(.rev2_t$skipped),
      " warn=", sum(.rev2_t$warning), "\n", sep = "")
  for (.rev2_i in seq_len(nrow(.rev2_t))) {
    if (.rev2_t$failed[.rev2_i] > 0 || isTRUE(.rev2_t$error[.rev2_i]) || .rev2_t$warning[.rev2_i] > 0 ||
          .rev2_t$skipped[.rev2_i]) {
      cat(sprintf("  TEST %s | fail=%d err=%s warn=%d skip=%s\n", .rev2_t$test[.rev2_i],
                  .rev2_t$failed[.rev2_i], .rev2_t$error[.rev2_i], .rev2_t$warning[.rev2_i], .rev2_t$skipped[.rev2_i]))
    }
  }
}
