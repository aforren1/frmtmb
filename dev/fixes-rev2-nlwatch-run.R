# Reviewer of lane fixes, re-check: run one test file with
# check_nl_identified() wrapped, logging every call that sees a
# nonlinear body, whether a pair of its nonlinear parameters has
# identical stats::D() derivatives, and whether the guard refused.
#   Rscript dev/fixes-rev-nlwatch-run.R <lib> <package> <test file> <log>
a <- commandArgs(trailingOnly = TRUE)
LIB <- a[1]
.libPaths(unique(c(LIB, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(testthat))
p <- a[2]; f <- a[3]; LOG <- a[4]
suppressMessages(library(p, character.only = TRUE))
ns <- asNamespace("frmtmb")
orig <- get("check_nl_identified", ns)
wrap <- function(spec, frame, prior) {
  lps <- frame[["linpreds"]] %||% list()
  bodies <- Filter(function(lp) !is.null(lp[["nl_body"]]), lps)
  if (length(bodies)) {
    same_pair <- FALSE
    for (lp in bodies) {
      pars <- intersect(lp[["nl_pars"]], all.vars(lp[["nl_body"]]))
      der <- lapply(pars, function(q) tryCatch(stats::D(lp[["nl_body"]], q),
                                               error = function(e) NULL))
      der <- der[!vapply(der, is.null, NA)]
      if (length(der) >= 2L) {
        for (i in seq_along(der)) for (j in seq_along(der)) {
          if (i < j && identical(der[[i]], der[[j]])) same_pair <- TRUE
        }
      }
    }
    body_txt <- paste(vapply(bodies, function(lp) deparse1(lp[["nl_body"]]),
                             ""), collapse = " ; ")
    res <- tryCatch({orig(spec, frame, prior); "pass"},
                    error = function(e) paste("REFUSED:",
                                              substr(conditionMessage(e), 1,
                                                     160)))
    cat(sprintf("%s\t%s\tsame_D=%s\tprior=%s\t%s\n", basename(f), body_txt,
                same_pair, !is.null(prior), res), file = LOG, append = TRUE)
    if (startsWith(res, "REFUSED")) orig(spec, frame, prior)
  }
  invisible(NULL)
}
environment(wrap) <- ns
assignInNamespace("check_nl_identified", wrap, ns = "frmtmb")
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
      "\n", sep = "")
}
