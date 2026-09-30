# Reviewer: hurdle_cumulative with disc ~ 1 + z (an intercept the
# likelihood cannot separate from the threshold scale). Seed 4400 data
# of dev/fams2-rev-dens.R section D.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib", "C:/Users/adf44/source/r/rellib-r3", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(4400)
n <- 1500
dD <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(30, n / 30))
disc <- exp(0.4 * dD$z)
u <- rlogis(n) / disc + 0.8 * dD$x
yc <- 1L + (u > -1) + (u > 0.3) + (u > 1.5)
dD$y <- ifelse(runif(n) < plogis(-0.6 + 0.5 * dD$x), 0L, yc)
w <- character()
f <- withCallingHandlers(frm(bf(y ~ x, disc ~ 1 + z), data = dD, family = hurdle_cumulative()),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
cat("warnings:", length(w), "\n"); if (length(w)) cat(w, sep = "\n")
print(fixef(f))
f0 <- frm(bf(y ~ x, disc ~ 0 + z), data = dD, family = hurdle_cumulative())
cat("logLik disc ~ 1 + z", as.numeric(logLik(f)), " disc ~ 0 + z", as.numeric(logLik(f0)), "\n")
