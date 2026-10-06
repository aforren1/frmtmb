# Reviewer: one test file, one R process, against the trial merge
# (nanse-rev-lib, then rellib-r6 for the extensions). Records every
# frm_warning() with the index of the top-level frm() call it came from,
# so a fit that warns twice about its standard errors shows up.
#   Rscript dev/nanse-rev-run-tests.R <package> <test file> <warn log>
LIB <- "C:/Users/adf44/source/r/nanse-rev-lib"
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
wlog <- a[3]
suppressMessages(library(p, character.only = TRUE))
cat("lib: ", find.package("frmtmb"), " | ", find.package(p), "\n", sep = "")
ns <- asNamespace("frmtmb")
.rev <- new.env()
.rev$fit <- 0L
.rev$depth <- 0L
.rev$log <- wlog
assign(".rev", .rev, envir = globalenv())
suppressMessages(trace(
  "frm", where = ns, print = FALSE,
  tracer = quote({
    e <- get(".rev", envir = globalenv())
    if (e$depth == 0L) e$fit <- e$fit + 1L
    e$depth <- e$depth + 1L
  }),
  exit = quote({
    e <- get(".rev", envir = globalenv())
    e$depth <- e$depth - 1L
  })))
suppressMessages(trace(
  "frm_warning", where = ns, print = FALSE,
  tracer = quote({
    e <- get(".rev", envir = globalenv())
    m <- tryCatch(paste0(..., collapse = ""), error = function(err) "?")
    m <- gsub("[\r\n\t]+", " ", m)
    cat(e$fit, "\t", e$depth, "\t", substr(m, 1, 300), "\n", sep = "",
        file = e$log, append = TRUE)
  })))
r <- tryCatch(
  as.data.frame(test_file(f, package = p, env = testthat::test_env(p),
                          reporter = "silent")),
  error = function(e) {
    cat("RESULT ", basename(f), " LOADERROR ", conditionMessage(e), "\n",
        sep = "")
    NULL
  })
if (!is.null(r)) {
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      " warn=", sum(r$warning), "\n", sep = "")
  for (i in seq_len(nrow(r))) {
    if (r$failed[i] > 0 || isTRUE(r$error[i]) || r$warning[i] > 0) {
      cat(sprintf("  TEST %s | fail=%d err=%s warn=%d\n", r$test[i],
                  r$failed[i], r$error[i], r$warning[i]))
    }
  }
}
