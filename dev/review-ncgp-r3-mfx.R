# Short check: predread's "0.46" (dev/review-predread-r1-check3.R, seed 7403): y ~ z + gp(x), predictions() against frm_linpred(se.fit),
# on this tree, on 0.78.0, and with gp(x, k = 15). Args: lane|base
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs(arm)
suppressMessages(library(marginaleffects))
set.seed(7403)
m <- 80
g <- data.frame(x = sort(runif(m, 0, 6)), z = rnorm(80))
g$y <- sin(g$x) + 0.3 * g$z + rnorm(m, 0, 0.25)
nd <- data.frame(x = c(g$x[5], 3.33, 7, 9), z = 0)
for (sp in list(bf(y ~ z + gp(x)), bf(y ~ z + gp(x, k = 15, c = 5/4)))) {
  f <- suppressWarnings(suppressMessages(frm(sp, family = gaussian(), data = g)))
  l <- frm_linpred(f, newdata = nd, se.fit = TRUE)
  p <- suppressWarnings(suppressMessages(predictions(f, newdata = nd)))
  cat(sprintf("RMX %s %s predictions se %s | frm_linpred se %s | conv %d\n", arm, deparse(sp$formula), paste(sprintf("%.3f", p$std.error), collapse = " "), paste(sprintf("%.3f", l$se.fit), collapse = " "), f$opt$convergence))
}
cat("DONE\n")
