# Reviewer of lane optima: one test file, one R process, on the lane's
# core (wt-optima-lib, then rellib-r6 for the extensions). Besides the
# RESULT line it logs every nlminb run that nlminb_best_par() changed
# (rejected_last), and for each top-level frm() that contained one, the
# reported objective against the objective re-evaluated at the fit's
# own par, so "every change is toward the reported objective" is a
# measurement over the suites.
#   Rscript dev/optima-rev-run-tests.R <package> <test file> <event log>
local({
LIB <- "C:/Users/adf44/source/r/wt-optima-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
suppressMessages(library(testthat))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
elog <- a[3]
suppressMessages(library(p, character.only = TRUE))
cat("lib: ", p, " ", as.character(utils::packageVersion(p)), " ",
    dirname(find.package(p)), "; frmtmb ",
    as.character(utils::packageVersion("frmtmb")), " ",
    dirname(find.package("frmtmb")), "\n", sep = "")
ns <- asNamespace("frmtmb")
rv <- new.env()
rv$fit <- 0L
rv$depth <- 0L
rv$runs <- 0L
rv$events <- 0L
rv$ev_fit <- 0L
rv$file <- basename(f)
rv$log <- elog
assign(".optrev", rv, envir = globalenv())
suppressMessages(trace(
  "nlminb_best_par", where = ns, print = FALSE,
  exit = quote({
    e <- get(".optrev", envir = globalenv())
    e$runs <- e$runs + 1L
    out <- returnValue()
    if (isTRUE(out$rejected_last)) {
      e$events <- e$events + 1L
      e$ev_fit <- e$ev_fit + 1L
      lst <- fnw$last()
      cat("EVENT\t", e$file, "\t", e$fit, "\t", e$depth, "\t",
          format(lst$val, digits = 17), "\t",
          format(out$objective, digits = 17), "\t",
          format(max(abs(as.numeric(lst$par) - as.numeric(out$par))),
                 digits = 4), "\t", length(out$par), "\t",
          res$convergence, "\t", gsub("[\t\n]", " ", res$message), "\n",
          sep = "", file = e$log, append = TRUE)
    }
  })))
suppressMessages(trace(
  "frm", where = ns, print = FALSE,
  tracer = quote({
    e <- get(".optrev", envir = globalenv())
    if (e$depth == 0L) {
      e$fit <- e$fit + 1L
      e$ev_fit <- 0L
    }
    e$depth <- e$depth + 1L
  }),
  exit = quote({
    e <- get(".optrev", envir = globalenv())
    e$depth <- e$depth - 1L
    if (e$depth == 0L && e$ev_fit > 0L) {
      ft <- returnValue()
      chk <- tryCatch({
        o <- ft$opt
        v <- as.numeric(ft$obj$fn(o$par))
        ll <- tryCatch(as.numeric(stats::logLik(ft)), error = function(x) NA)
        paste(format(o$objective, digits = 17), format(v, digits = 17),
              format(ll, digits = 17), o$convergence,
              isTRUE(o$rejected_last), sep = "\t")
      }, error = function(x) paste("ERR", conditionMessage(x)))
      cat("FIT\t", e$file, "\t", e$fit, "\t", e$ev_fit, "\t", chk, "\n",
          sep = "", file = e$log, append = TRUE)
    }
  })))
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
cat("TOTAL\t", rv$file, "\tfits ", rv$fit, "\tnlminb_runs ", rv$runs,
    "\tevents ", rv$events, "\n", sep = "", file = elog, append = TRUE)
if (!is.null(r)) {
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      " warn=", sum(r$warning), "\n", sep = "")
}
})
