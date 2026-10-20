# Release 0.72.0 consolidation: copy of lane adhess's battery script, with the "lane" arm on wt-release-lib (the merged build) and "base" on rellib-r9 (0.71.0); paths point at the release tree.
# Lane adhess, copied from lane illse2's copy (from dev/illse-sf-rev3-common.R): copy of dev/rel070-seflags-rev3-common.R by dev/illse-sf-make.R (lane arm: the lane library over
# rellib-r9; base arm: rellib-r9).
# Reviewer of lane seflags (final check): sourced by dev/seflags-rev3-*.R.
# First argument: lane | base | op3 (seflags + optrap tree of 2026-10-12 23:49) | m9 (the lane's own merge build). Sets .libPaths(), loads frmtmb, defines cap().
rev_arm <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(rev_arm)) rev_arm <- "lane"
.libPaths(c(switch(rev_arm,
  lane = c("C:/Users/adf44/source/r/wt-release-lib",
           "C:/Users/adf44/source/r/rellib-r9"),
  base = "C:/Users/adf44/source/r/rellib-r9",
  op3 = c("C:/Users/adf44/source/r/seflags-rev-lib/op3",
          "C:/Users/adf44/source/r/rellib-r7"),
  m9 = c("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-merge-out/lib",
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
