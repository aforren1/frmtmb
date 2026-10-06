# Lane fixes, punch round: what a user sees on a smooth fit whose
# fixed-effect standard errors are not finite (seed 3, s(x1) + s(x2)):
# the warnings of summary(), vcov() and diagnose().
#   Rscript dev/fixes-sx-conv5.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
set.seed(3)
d <- suppressMessages(mgcv::gamSim(eg = 6, n = 200, scale = 2,
                                   verbose = FALSE))
grab <- function(expr) {
  w <- character()
  withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  if (length(w)) substr(w, 1, 160) else "(none)"
}
for (k in 1:3) {
  fit <- suppressWarnings(frm(bf(y ~ s(x1) + s(x2)), data = d))
  cat("call", k, "\n")
  if (k == 1) cat("  fixef  :", grab(fixef(fit)), sep = "\n    ")
  if (k == 2) cat("  summary:", grab(summary(fit)), sep = "\n    ")
  if (k == 3) cat("  vcov   :", grab(vcov(fit)), sep = "\n    ")
}
cat("sdreport pdHess:", sdr_of(fit)$pdHess, "\n")
print(round(fit$estimates$theta, 3))
