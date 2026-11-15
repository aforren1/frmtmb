# Common setup of lane nlstop2's scripts (round of 2026-11-13).
# NLSTOP2_ARM = before reads rellib-r21 (frmtmb 0.82.0 at ec8178d3);
# after (the default) puts the lane's library first.
arm <- Sys.getenv("NLSTOP2_ARM", "after")
libs <- c(if (arm == "after") "C:/Users/adf44/source/r/rellib-r22",
          if (arm != "after") "C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(libs)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/stan-cache")
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
msum <- function(r) {
  mc <- attr(r, "mc")
  if (length(mc)) paste(unique(mc), collapse = "+") else "-"
}
# the fresh generator of dev/nlstop-fresh.R (F, pos, B, gen, form)
fresh_defs <- function(env = parent.frame()) {
  src <- parse("dev/nlstop-fresh.R")
  keep <- vapply(src, function(e) {
    is.call(e) && identical(e[[1L]], as.name("<-")) &&
      as.character(e[[2L]])[1] %in% c("F", "pos", "B", "gen", "form")
  }, NA)
  for (e in src[keep]) eval(e, env)
  invisible(NULL)
}
