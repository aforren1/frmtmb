# Review of lane nlstart: library setup shared by dev/review-nlstart-*.R.
# REV_ARM = after (lane lib first), before (rellib-r18), or a path to a
# mutant library put first.
arm <- Sys.getenv("REV_ARM", "after")
libs <- c(switch(arm,
                 after = "C:/Users/adf44/source/r/wt-nlstart-lib",
                 before = character(0),
                 arm),
          "C:/Users/adf44/source/r/rellib-r18",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(libs)
Sys.setenv(FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-nlstart/dev/stan-cache")
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
