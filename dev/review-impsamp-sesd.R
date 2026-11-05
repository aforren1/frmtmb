# nb_mushape seed 2 read a theta_2 SE 0.90 of the exact one at 500 draws. Noise of
# the Monte Carlo Hessian, or a bias a scalar curvature rate cannot remove with two
# blocks? Over draw seeds and draw counts, on the same data as the multid design.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
set.seed(2)
g <- rep(1:20, each = 6)
x <- rnorm(120)
u <- rnorm(20, 0, 0.5)
v <- rnorm(20, 0, 0.4)
y <- rnbinom(120, size = exp(0.7 + v[g]), mu = exp(1 + 0.4 * x + u[g]))
d <- data.frame(y, x, g = factor(g))
f <- bf(y ~ x + (1 | g), shape ~ 1 + (1 | g))
se <- function(N, s) {
  fi <- suppressWarnings(frm(f, family = negbinomial(), data = d, importance = N,
                             control = frmtmb_control(importance_seed = s)))
  c(sqrt(diag(frmtmb:::sdr_of(fi)$cov.fixed))[4:5], fi$opt$par[4:5], fi$importance$ess_min)
}
for (N in c(500L, 4000L)) {
  m <- t(vapply(1:6, function(s) se(N, s), numeric(5)))
  colnames(m) <- c("se_th1", "se_th2", "th1", "th2", "ess")
  cat("N", N, "\n")
  print(round(m, 4))
}
