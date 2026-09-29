# REVIEW runner: ONE test file, ONE process.
#
#   Rscript dev/arcovsample-rev-run.R <lane|ref> <pkg> <test file> [gated]
#
# `lane` puts the worker's build first, `ref` puts rellib-r3 first and
# does NOT put the lane library on the path at all. Which build was
# actually used is PRINTED from a symbol only the lane build has, because
# both builds carry the same version number and a wrong path order would
# otherwise be invisible.

a <- commandArgs(trailingOnly = TRUE)
arm <- a[1L]; pkg <- a[2L]; f <- a[3L]
gated <- length(a) >= 4L && identical(a[4L], "gated")

LANE <- "C:/Users/adf44/source/r/wt-arcovsample-lib"
REF <- "C:/Users/adf44/source/r/rellib-r3"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(LANE, REF, USER) else c(REF, USER))

mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
stopifnot(utils::packageVersion("tmbstan") >= "1.2.1",
          any(grepl("-std=gnu++17",
                    readLines(tools::makevars_user(), warn = FALSE),
                    fixed = TRUE)))
Sys.setenv(NOT_CRAN = "true")
if (gated) Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")

suppressMessages(library(testthat))
suppressMessages(library(pkg, character.only = TRUE))
cat("ARM ", arm, " gated=", gated, " frmtmb=",
    format(packageVersion("frmtmb")), " at ",
    dirname(system.file("DESCRIPTION", package = "frmtmb")),
    " newexports=",
    all(c("arma_cond_resp", "arma_cond_dpars") %in%
          getNamespaceExports("frmtmb")), "\n", sep = "")

r <- tryCatch(
  as.data.frame(test_file(f, package = pkg, env = testthat::test_env(pkg),
                          reporter = "silent")),
  error = function(e) {
    cat("RESULT ", basename(f), " LOADERROR ", conditionMessage(e), "\n",
        sep = "")
    NULL
  })

if (!is.null(r)) {
  cat("RESULT ", arm, if (gated) " gated" else "", " ", basename(f),
      " pass=", sum(r$passed), " fail=", sum(r$failed), " err=",
      sum(r$error), " skip=", sum(r$skipped), "\n", sep = "")
  bad <- r[r$failed > 0 | r$error > 0, , drop = FALSE]
  for (i in seq_len(nrow(bad))) {
    cat("--- BAD: ", bad$test[i], "\n", sep = "")
    for (rs in bad$result[[i]]) {
      if (inherits(rs, c("expectation_failure", "expectation_error"))) {
        cat("    ", gsub("\n", "\n    ", conditionMessage(rs)), "\n",
            sep = "")
      }
    }
  }
  sk <- r[r$skipped > 0, , drop = FALSE]
  for (i in seq_len(nrow(sk))) {
    cat("--- SKIP: ", sk$test[i], "\n", sep = "")
    for (rs in sk$result[[i]]) {
      if (inherits(rs, "expectation_skip")) {
        cat("    ", conditionMessage(rs), "\n", sep = "")
      }
    }
  }
}
