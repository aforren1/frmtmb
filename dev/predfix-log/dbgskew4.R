source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
mk <- function(seed, s, n = 250) { set.seed(seed); xs <- -abs(rnorm(n)) * 3
  yy <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3); data.frame(y = yy, xs = xs * s, xr = xs) }
fo <- bf(y ~ xr, sigma ~ 1, alpha ~ 0 + xs)
for (seed in c(4, 10)) {
  d <- mk(seed, 1e-6)
  for (a in list(FALSE, TRUE)) {
    w <- character(0)
    f <- withCallingHandlers(frm(fo, family = skew_normal(), data = d, control = frmtmb_control(autoscale = a)),
      warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
    cat("seed", seed, "autoscale", a, "logLik", format(as.numeric(logLik(f)), digits = 10), "code", f$opt$convergence, "warnings:", length(w), "\n")
    for (x in w) cat("   ", substr(x, 1, 120), "\n")
  }
}
