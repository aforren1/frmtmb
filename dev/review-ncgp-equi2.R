# Item 2: the NaN standard errors that appear at y * 1000 and not at y, and the warnings.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
arm <- commandArgs(TRUE)[1]
rev_libs(arm)
one <- function(y, x, ysc, xsc, lab) {
  d <- data.frame(y = y * ysc, x = x * xsc)
  w <- character()
  f <- withCallingHandlers(frm(y ~ gp(x), d), warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  w2 <- character()
  ft <- withCallingHandlers(fitted(f), warning = function(c) { w2 <<- c(w2, conditionMessage(c)); invokeRestart("muffleWarning") })
  V <- as.matrix(vcov(f, full = TRUE))
  ses <- sqrt(diag(V))
  cat(sprintf("REQ %s ysc %g xsc %g conv %d msg %s ll+nlogc %.5f | se(full) %s | fitted NaN SE %d of %d | hess cond? | warnings fit %d fitted %d\n", lab, ysc, xsc, f$opt$convergence,
              f$opt$message, logLik(f) + length(y) * log(ysc), paste(names(ses), signif(ses, 4), collapse = " "), sum(is.nan(ft[, 2])), nrow(ft), length(w), length(w2)))
  for (m in unique(c(w, w2))) cat("REQW", substr(m, 1, 230), "\n")
}
for (cs in list(list(s = 1e-4, om = 1, seed = 14, n = 40), list(s = 1e-3, om = 2, seed = 13, n = 80))) {
  set.seed(cs$seed)
  x <- sort(runif(cs$n, 0, 8))
  y <- sin(cs$om * x) + cs$s * rnorm(cs$n)
  for (sc in list(c(1, 1), c(1000, 1), c(1e-3, 1), c(1e6, 1), c(1e-6, 1), c(1, 1000), c(1e6, 1e-6))) one(y, x, sc[1], sc[2], sprintf("s=%g", cs$s))
}
cat("DONE\n")
