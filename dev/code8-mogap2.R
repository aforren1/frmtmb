# The same gap fit stopped short by an iteration cap: does get_vcov()
# still take its "lost another direction" branch?
source("dev/code8-common.R")
set.seed(23)
n <- 300
d <- data.frame(x = sample(c(1, 2, 5), n, TRUE), z = rnorm(n),
                g = factor(sample(1:10, n, TRUE)))
d$y <- 1 + (0.4 + rnorm(10, 0, 0.2)[d$g]) * d$x + 0.5 * d$z + rnorm(n)
for (it in c(3, 5, 10, 20)) {
  r <- rfit(y ~ mo(x) + z + (mo(x) | g), data = d,
            control = frmtmb_control(optCtrl = list(iter.max = it,
                                                    eval.max = 1000)))
  cat("== iter.max", it, "code", r$f$opt$convergence, "\n")
  for (w in r$w) cat("  WARN:", substr(w, 1, 120), "\n")
  V <- withCallingHandlers(marginaleffects::get_vcov(r$f),
    warning = function(w) { cat("  GETVCOV WARN:", substr(conditionMessage(w),
                                                           1, 100), "\n")
      invokeRestart("muffleWarning") })
  cat("  zeta1_cum1 row NaN:", all(is.nan(V["zeta1_cum1", ])), "\n")
}
