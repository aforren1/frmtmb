# Review of lane nlstop2: common setup. REV_ARM = before reads rellib-r21
# (0.82.0); after (default) puts the lane's library first.
arm <- Sys.getenv("REV_ARM", "after")
.libPaths(c(if (arm == "after") "C:/Users/adf44/source/r/wt-nlstop2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "arm:", arm, "\n")
ns <- asNamespace("frmtmb")
# value, warning classes and messages' classes of expr
tt <- function(expr) {
  w <- character(0); m <- character(0); wt <- character(0)
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(c) {
      w <<- c(w, class(c)[1]); wt <<- c(wt, conditionMessage(c))
      invokeRestart("muffleWarning")
    },
    message = function(c) {
      m <<- c(m, class(c)[1]); invokeRestart("muffleMessage")
    }),
    error = function(e) structure(list(msg = conditionMessage(e)),
                                  class = "tterr"))
  list(v = v, w = w, m = m, wt = wt)
}
show <- function(tag, r) {
  if (inherits(r$v, "tterr")) {
    cat("ROW", tag, "ERROR", substr(r$v$msg, 1, 220), "\n")
  } else {
    cat("ROW", tag, "ok", "warn", if (length(r$w)) paste(unique(r$w),
        collapse = "+") else "-", "msg", if (length(r$m))
        paste(unique(r$m), collapse = "+") else "-", "\n")
    if (length(r$wt)) cat("   W:", substr(r$wt[1], 1, 200), "\n")
  }
}
llv <- function(f) {
  if (inherits(f, "tterr")) return(NA_real_)
  suppressWarnings(as.numeric(logLik(f)))
}
