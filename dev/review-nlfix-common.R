# Common setup of the nlfix review scripts. FRM_ARM=before reads rellib-r16.
arm <- Sys.getenv("FRM_ARM", "after")
libs <- c(if (arm == "after") "C:/Users/adf44/source/r/wt-nlfix-lib2",
          if (arm == "mutant") Sys.getenv("REV_LIB", "C:/Users/adf44/source/r/rev-nlfix-lib"),
          "C:/Users/adf44/source/r/rellib-r16",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(libs)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-nlfix/dev/stan-cache")
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), as.character(packageVersion("frmtmb")), "\n")
tt <- function(expr) {
  w <- character(0)
  m <- character(0)
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") },
    message = function(c) { m <<- c(m, conditionMessage(c)); invokeRestart("muffleMessage") }),
    error = function(e) structure(list(msg = conditionMessage(e)), class = "tterr"))
  if (is.null(v)) v <- list()
  attr(v, "w") <- w
  attr(v, "m") <- m
  v
}
iserr <- function(x) inherits(x, "tterr")
show <- function(label, x) {
  if (iserr(x)) cat(label, ": ERROR ", substr(x$msg, 1, 300), "\n", sep = "")
  else {
    cat(label, ": ok", sep = "")
    if (is.numeric(x) && length(x) <= 8) cat(" ", format(as.numeric(x), digits = 10))
    cat("\n")
  }
  if (length(attr(x, "w"))) cat("  warn:", substr(attr(x, "w"), 1, 300), sep = "\n   ")
  if (length(attr(x, "m"))) cat("  msg:", substr(attr(x, "m"), 1, 300), sep = "\n   ")
  if (length(attr(x, "w")) || length(attr(x, "m"))) cat("\n")
  invisible(x)
}
