# Noise-free y on the closed form at relative 1e-12, seed 20261102, against frm()'s
# stop, and frm() at growing noise down to where it stops reaching it.
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261102)
n <- 60L
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
D2 <- outer(x, x, "-")^2
marg <- function(p, y) {
  K <- exp(2 * p[2]) * exp(-D2 / (2 * exp(2 * p[3])))
  diag(K) <- diag(K) * (1 + 1e-12) + exp(2 * p[4])
  R <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(R)) return(1e12)
  r <- backsolve(R, y - p[1], transpose = TRUE)
  0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * n * log(2 * pi)
}
for (s in c(1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 1e-7, 1e-8, 0)) {
  y <- sin(2 * x) + s * e0
  d <- data.frame(x = x, y = y)
  f <- suppressWarnings(frm(bf(y ~ gp(x)), data = d, family = gaussian()))
  best <- NULL
  for (st in list(c(0, 1, 0.3, log(max(s, 1e-9))), c(0, 2, 0.5, log(max(s, 1e-9))), c(0, 3, 0.7, log(max(s, 1e-7))))) {
    o <- nlminb(st, marg, y = y, control = list(iter.max = 1000, eval.max = 2000, rel.tol = 1e-12))
    if (is.null(best) || o$objective < best$objective) best <- o
  }
  cat(sprintf("NCGPNF2 s=%.0e frm ll %.3f conv %d sigma %.3g sdgp %.3g grad %.2g | cf ll %.3f sigma %.3g sdgp %.3g | gap %.3g\n",
              s, -f$opt$objective, f$opt$convergence, sigma(f)[1], exp(f$estimates$theta[1]), max(abs(f$obj$gr(f$opt$par))),
              -best$objective, exp(best$par[4]), exp(best$par[2]), -f$opt$objective + best$objective))
}
