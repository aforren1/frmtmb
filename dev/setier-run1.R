# Lane setier: one test file, one R process, against a library stack,
# recording every condition frmtmb raises (frm_warning(), frm_message())
# even where a test muffles it, so that a message's firing count is read
# from the run and not from the escaped warnings.
#   Rscript dev/setier-run1.R <libs, ';'-separated> <package> <file> <cond log>
# The user library is always appended last. Based on dev/release/run-tests.R.
local({
a <- commandArgs(trailingOnly = TRUE)
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(strsplit(a[1], ";", fixed = TRUE)[[1]], user))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
suppressMessages(library(testthat))
p <- a[2]
f <- a[3]
cond <- a[4]
suppressMessages(library(p, character.only = TRUE))
cat("lib: ", p, " ", as.character(utils::packageVersion(p)), " ",
    dirname(find.package(p)), "; frmtmb ",
    as.character(utils::packageVersion("frmtmb")), " ",
    dirname(find.package("frmtmb")), "; BLAS probe ",
    system.time({m <- matrix(1, 600, 600); m %*% m})[["elapsed"]], "\n",
    sep = "")
ns <- asNamespace("frmtmb")
logger <- function(kind) {
  substitute({
    .txt <- tryCatch({
      .a <- list(...)
      paste(vapply(.a, function(x) {
        if (inherits(x, "condition")) conditionMessage(x) else
          if (is.atomic(x)) paste(x, collapse = "") else ""
      }, ""), collapse = "")
    }, error = function(e) "?")
    cat(KIND, "\t", gsub("[\r\n\t]+", " ", .txt), "\n", sep = "",
        file = COND, append = TRUE)
  }, list(KIND = kind, COND = cond))
}
suppressMessages({
  trace("frm_warning", where = ns, print = FALSE, tracer = logger("W"))
  trace("frm_message", where = ns, print = FALSE, tracer = logger("M"))
})
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
