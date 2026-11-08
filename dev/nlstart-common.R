# Common setup of lane nlstart's scripts. NLSTART_ARM=before reads
# rellib-r18 (frmtmb 0.79.1 at 3cd2a0c8); "after" puts the lane's library
# first.
arm <- Sys.getenv("NLSTART_ARM", "after")
libs <- c(if (arm == "after") "C:/Users/adf44/source/r/wt-nlstart-lib",
          "C:/Users/adf44/source/r/rellib-r18",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(libs)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-nlstart/dev/stan-cache")
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "arm:", arm, "\n")
tt <- function(expr) {
  w <- character(0)
  wc <- character(0)
  m <- character(0)
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(c) {
      w <<- c(w, conditionMessage(c))
      wc <<- c(wc, class(c)[1])
      invokeRestart("muffleWarning")
    },
    message = function(c) {
      m <<- c(m, conditionMessage(c))
      invokeRestart("muffleMessage")
    }),
    error = function(e) structure(list(msg = conditionMessage(e)),
                                  class = "tterr"))
  if (is.null(v)) v <- list()
  attr(v, "w") <- w
  attr(v, "wc") <- wc
  attr(v, "m") <- m
  v
}
iserr <- function(x) inherits(x, "tterr")
ll <- function(r) {
  if (iserr(r)) return(paste("ERR", substr(r$msg, 1, 200)))
  paste0(format(as.numeric(logLik(r)), digits = 10),
         " conv=", r$opt$convergence,
         " beta=", paste(format(r$estimates$beta, digits = 5), collapse = ","),
         " w=", length(attr(r, "w")), " [", paste(attr(r, "wc"),
                                                  collapse = ","), "]",
         " m=", length(attr(r, "m")))
}
