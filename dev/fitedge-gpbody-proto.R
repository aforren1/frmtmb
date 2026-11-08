# Item 1 prototype: does a theta start alone carry a body's fit to the linear
# spelling's maximum? Oracle start = the linear fit's theta. Args: arm.
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
mk <- function(seed, n, curve) {
  set.seed(seed)
  x <- sort(runif(n, 0, 10)); z <- rnorm(n)
  f <- switch(curve, sin = sin(x), wig = 2 * sin(3 * x) + x,
              exp = exp(0.3 * x), step = sign(sin(1.5 * x)))
  data.frame(x, z, y = f + 0.3 * z + 0.1 * rnorm(n))
}
for (cs in list(list(1811, 60, "wig"), list(58001, 60, "wig"))) {
  d <- mk(cs[[1]], cs[[2]], cs[[3]])
  lin <- fe_fit(y ~ z + gp(x), data = d)
  nlf0 <- bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE)
  nl <- fe_fit(nlf0, data = d, start = list(theta = lin$fit$estimates$theta))
  # the working response a theta-only start would read: y itself
  dd <- data.frame(x = d$x, y = d$y)
  w <- fe_fit(y ~ gp(x), data = dd)
  nl2 <- fe_fit(nlf0, data = d, start = list(theta = w$fit$estimates$theta))
  cat(sprintf("P2 %s seed %d lin %.4f | nl oracle theta %.4f conv %d | nl y-only theta %.4f conv %d (theta %s)\n",
              arm, cs[[1]], logLik(lin$fit), logLik(nl$fit), nl$fit$opt$convergence,
              logLik(nl2$fit), nl2$fit$opt$convergence,
              paste(signif(w$fit$estimates$theta, 3), collapse = " ")))
}
cat("DONE\n")
