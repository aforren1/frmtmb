# Which conditions predict() raises on lane drawdir's mix40 cumulative-
# mixture fits on a build (the 0.84.0 key-table compare moved 8 of 40 from
# a draw verdict to "frmtmb_warning"), and the fit's own conditions.
#   Rscript dev/rel084-kt/mix40-warn.R <lib> <seed ...>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
setwd("C:/Users/adf44/source/r/frmtmb-wt-rel084")
suppressPackageStartupMessages(library(frmtmb))
source("dev/drawdir-cases.R")
cat("lib:", as.character(packageVersion("frmtmb")), dirname(find.package("frmtmb")), "\n")
grab <- function(expr) {
  w <- list()
  withCallingHandlers(tryCatch(expr, error = function(e) w[[length(w) + 1L]] <<- e),
    warning = function(c) { w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning") },
    message = function(c) { w[[length(w) + 1L]] <<- c; invokeRestart("muffleMessage") })
  for (c in w) cat("    ", paste(head(class(c), 2), collapse = "/"), ":",
                   substr(gsub("[[:space:]]+", " ", conditionMessage(c)), 1, 300), "\n")
}
for (s in a[-1]) {
  id <- paste0("mix40:mixture_cumulative:interior:", s)
  cat("==", id, "\n  fit:\n")
  r <- NULL
  grab(r <- dc_fit(id))
  fit <- r$fit
  cat("  logLik", format(as.numeric(logLik(fit)), digits = 12), "code", fit$opt$convergence, "\n  predict:\n")
  set.seed(5)
  grab(predict(fit, ndraws = 200, summary = FALSE))
}
