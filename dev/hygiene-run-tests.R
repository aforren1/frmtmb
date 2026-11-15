# One test file, one R process, for lane hygiene. A copy of
# dev/release/run-tests.R with the library chosen by HYG_ARM:
#   mine  the lane's library, then rellib-r21, then the user library
#   base  rellib-r21, then the user library (the "before" arm)
# HYG_BARE=true attaches the package but runs test_file() without
# package= and env=, the way a reader runs one file by hand against the
# installed package (backlog items 764 and 769).
#
# Usage: Rscript dev/hygiene-run-tests.R <package> <path-to-test-file>
local({
MINE <- "C:/Users/adf44/source/r/wt-hygiene-lib"
BASE <- "C:/Users/adf44/source/r/rellib-r21"
USR <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
arm <- Sys.getenv("HYG_ARM", "mine")
.libPaths(c(if (arm == "mine") MINE, BASE, USR))

mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}

suppressMessages(library(testthat))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
bare <- identical(Sys.getenv("HYG_BARE"), "true")
suppressMessages(library(p, character.only = TRUE))
# HYG_PRELOAD="a,b" loads these first, as an earlier file of one session
# would have
for (q in strsplit(Sys.getenv("HYG_PRELOAD"), ",", fixed = TRUE)[[1]]) {
  if (nzchar(q)) suppressMessages(loadNamespace(q))
}
cat("lib: arm=", arm, " bare=", bare, " ", p, " ",
    as.character(utils::packageVersion(p)), " ",
    dirname(find.package(p)), "; frmtmb ",
    as.character(utils::packageVersion("frmtmb")), " ",
    dirname(find.package("frmtmb")), "\n", sep = "")

res <- NULL
r <- tryCatch({
  res <- if (bare) {
    test_file(f, reporter = "silent")
  } else {
    test_file(f, package = p, env = testthat::test_env(p),
              reporter = "silent")
  }
  as.data.frame(res)
}, error = function(e) {
  cat("RESULT ", basename(f), " LOADERROR ", conditionMessage(e), "\n",
      sep = "")
  NULL
})

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
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      " warn=", sum(r$warning), "\n", sep = "")
}
})
