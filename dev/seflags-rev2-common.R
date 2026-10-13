# Reviewer of lane seflags (re-check): sourced by dev/seflags-rev2-*.R.
# First argument: lane | base | op2 (seflags + optrap current tree). Sets .libPaths(), loads frmtmb, defines cap().
rev_arm <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(rev_arm)) rev_arm <- "lane"
.libPaths(c(switch(rev_arm,
  lane = c("C:/Users/adf44/source/r/wt-seflags-lib",
           "C:/Users/adf44/source/r/rellib-r7"),
  base = "C:/Users/adf44/source/r/rellib-r7",
  op2 = c("C:/Users/adf44/source/r/seflags-rev-lib/op2",
          "C:/Users/adf44/source/r/rellib-r7")),
  "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("arm", rev_arm, "frmtmb from", find.package("frmtmb"),
    "| BLAS", basename(dirname(dirname(R.home("bin")))), R.home(), "\n")
cap <- function(expr) {
  w <- character()
  m <- character()
  v <- withCallingHandlers(tryCatch(expr, error = function(e) {
    structure(conditionMessage(e), class = "rev_err")
  }), warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(x) {
    m <<- c(m, conditionMessage(x))
    invokeRestart("muffleMessage")
  })
  list(v = v, w = w, m = m)
}
lost_of <- function(fit) {
  l <- suppressWarnings(ns$sdr_of(fit)$se_lost)
  if (!length(l)) return("-")
  paste(paste0(names(l), ":", l), collapse = " ")
}
se_outer <- function(fit) {
  v <- suppressWarnings(vcov(fit, full = TRUE))
  stats::setNames(sqrt(diag(as.matrix(v))), ns$outer_par_names(fit))
}
