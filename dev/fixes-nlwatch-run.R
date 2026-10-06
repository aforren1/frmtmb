# Lane fixes, punch round 2 (after the reviewer's dev/fixes-rev-nlwatch-run.R):
# one test file with the fitted-point flat check nl_flat_message()
# wrapped, logging every call and whether it warned, and every offer of
# a recovery start.
#   Rscript dev/fixes-nlwatch-run.R <lib> <package> <test file> <log>
a <- commandArgs(trailingOnly = TRUE)
LIB <- a[1]
.libPaths(unique(c(LIB, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(testthat))
p <- a[2]; f <- a[3]; LOG <- a[4]
suppressMessages(library(p, character.only = TRUE))
ns <- asNamespace("frmtmb")
orig <- get("nl_flat_message", ns)
wrap <- function(obj, opt, frame) {
  bodies <- Filter(function(lp) !is.null(lp[["nl_body"]]),
                   frame[["linpreds"]] %||% list())
  res <- orig(obj, opt, frame)
  if (length(bodies)) {
    body_txt <- paste(vapply(bodies, function(lp) deparse1(lp[["nl_body"]]),
                             ""), collapse = " ; ")
    cat(sprintf("%s\t%s\t%s\n", basename(f), body_txt,
                if (is.null(res)) "pass" else
                  paste("WARNS:", substr(res, 1, 120))),
        file = LOG, append = TRUE)
  }
  res
}
environment(wrap) <- ns
assignInNamespace("nl_flat_message", wrap, ns = "frmtmb")
orig_rs <- get("fit_recovery_starts", ns)
wrap_rs <- function(...) {
  out <- orig_rs(...)
  cat(sprintf("%s\tRECOVERY\t%s\n", basename(f),
              paste(names(out), collapse = " | ")),
      file = LOG, append = TRUE)
  out
}
environment(wrap_rs) <- ns
assignInNamespace("fit_recovery_starts", wrap_rs, ns = "frmtmb")
res <- tryCatch(test_file(f, package = p, env = testthat::test_env(p),
                          reporter = "silent"),
                error = function(e) {
                  cat("RESULT ", basename(f), " LOADERROR ",
                      conditionMessage(e), "\n", sep = ""); NULL
                })
if (!is.null(res)) {
  r <- as.data.frame(res)
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      " warn=", sum(r$warning), "\n", sep = "")
}
