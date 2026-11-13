# Common setup of lane nlstop's scripts (round of 2026-11-11).
# NLSTOP_ARM = before reads rellib-r20 (frmtmb 0.81.0 at e8abed0a);
# after (the default) puts the lane's library first.
arm <- Sys.getenv("NLSTOP_ARM", "after")
libs <- c(if (arm == "after") "C:/Users/adf44/source/r/wt-nlstop-lib",
          "C:/Users/adf44/source/r/rellib-r20",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(libs)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-nlstop/dev/stan-cache")
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "arm:", arm, "\n")
# run expr, keep its value, warnings (class and text) and messages
tt <- function(expr) {
  w <- character(0); wc <- character(0); m <- character(0); mc <- character(0)
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(c) {
      w <<- c(w, conditionMessage(c)); wc <<- c(wc, class(c)[1])
      invokeRestart("muffleWarning")
    },
    message = function(c) {
      m <<- c(m, conditionMessage(c)); mc <<- c(mc, class(c)[1])
      invokeRestart("muffleMessage")
    }),
    error = function(e) structure(list(msg = conditionMessage(e)),
                                  class = "tterr"))
  if (is.null(v)) v <- list()
  attr(v, "w") <- w; attr(v, "wc") <- wc
  attr(v, "m") <- m; attr(v, "mc") <- mc
  v
}
iserr <- function(x) inherits(x, "tterr")
llv <- function(r) {
  if (iserr(r) || is.null(r$opt)) return(NA_real_)
  v <- suppressWarnings(as.numeric(logLik(r)))
  if (is.finite(v)) v else NA_real_
}
llp <- function(r) {
  if (iserr(r) || is.null(r$opt)) return(NA_real_)
  v <- suppressWarnings(as.numeric(logLik_penalized(r)))
  if (is.finite(v)) v else NA_real_
}
wsum <- function(r) {
  wc <- attr(r, "wc")
  if (length(wc)) paste(unique(wc), collapse = "+") else "-"
}
