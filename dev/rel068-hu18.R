# The one new SE warning of the ordmix false-alarm table on the merged
# tree: hurdle_huz seed 18 (dev/rel068-falsealarm.R's model and data).
# Base (rellib-r5) and release (rellib-r6): the fit, sdreport()'s SE of
# every outer parameter, and the release's verdict.
#   Rscript dev/rel068-hu18.R r5|r6
arm <- commandArgs(TRUE)[1]
lib <- paste0("C:/Users/adf44/source/r/rellib-", arm)
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-release/dev/rel068-falsealarm.R")
g0 <- grep("^cut4 <- ", src); g1 <- grep("^spec <- switch", src) - 1L
eval(parse(text = src[g0:g1]))
d <- gen(18)
w <- character()
fit <- withCallingHandlers(
  frm(bf(yh ~ x, hu1 ~ z, hu2 ~ z),
      family = mixture(hurdle_cumulative(), hurdle_cumulative()), data = d),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
cat(arm, "logLik", format(as.numeric(logLik(fit)), digits = 12), "code", fit$opt$convergence, "\n")
for (x in w) cat("  W:", substr(x, 1, 300), "\n")
sdr <- RTMB::sdreport(fit$obj, fit$opt$par)
se <- sqrt(diag(sdr$cov.fixed))
cat("sdreport SE by parameter:\n"); print(signif(setNames(se, make.unique(names(fit$opt$par))), 4))
H <- optimHess(fit$opt$par, fit$obj$fn, fit$obj$gr)
D <- 1 / sqrt(abs(diag(H)))
ev <- eigen(H * outer(D, D), symmetric = TRUE)$values
cat("unit-diagonal Hessian eigenvalues, smallest 3:", signif(tail(ev, 3), 4), " largest", signif(ev[1], 4), "\n")
fx <- suppressWarnings(fixef(fit))
print(signif(fx[grepl("Intercept", rownames(fx)), 1:2], 4))
