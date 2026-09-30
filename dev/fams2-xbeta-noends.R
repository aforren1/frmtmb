# xbeta on a response with no exact 0 or 1: where does kappa go, and
# does the fit say anything? Seed 31.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib", "C:/Users/adf44/source/r/rellib-r3", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(31)
d <- data.frame(x = rnorm(400))
mu <- plogis(0.2 + 0.5 * d$x)
d$y <- rbeta(400, mu * 5, (1 - mu) * 5)
w <- character()
f <- withCallingHandlers(frm(bf(y ~ x), family = xbeta(), data = d),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
cat("warnings:", length(w), "\n"); if (length(w)) cat(w, sep = "\n")
print(confint(f))
fb <- frm(bf(y ~ x), family = Beta(), data = d)
cat("logLik xbeta", logLik(f), " Beta", logLik(fb), "\n")
