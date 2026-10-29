## Review of lane rescorpp: library setup. First argument: lane | base | rv.
## lane: the lane library over rellib-r15; base: rellib-r15 alone;
## rv: the reviewer's library (lane frmtmb.sample only) over rellib-r15.
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "lane"
base <- c("C:/Users/adf44/source/r/rellib-r15",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-rescorpp-lib", base),
  rv = c("C:/Users/adf44/source/r/wt-rescorpp-rv-lib", base),
  base))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-rescorpp"
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(wt, "dev/stan-cache"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("arm:", arm, " lib:", dirname(find.package("frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "| sample",
    dirname(find.package("frmtmb.sample")),
    format(packageVersion("frmtmb.sample")), "\n")
options(width = 130)
rv_out <- file.path(wt, "dev/review-rescorpp-out")
dir.create(rv_out, showWarnings = FALSE)
show <- function(tag, expr) {
  cat("\n==", tag, "==\n")
  w <- character(0)
  r <- withCallingHandlers(tryCatch(expr, error = function(e) e),
                           warning = function(x) {
                             w <<- c(w, conditionMessage(x))
                             invokeRestart("muffleWarning")
                           },
                           message = function(m) invokeRestart("muffleMessage"))
  if (inherits(r, "error")) cat("ERROR:", conditionMessage(r), "\n")
  else print(r)
  if (length(w)) cat("warnings:", length(w), "\n")
  for (m in w) cat("WARNING:", substr(m, 1, 300), "\n")
  invisible(r)
}
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
