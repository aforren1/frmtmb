# Item 3: the condition bound gp_marginal_outer_gr() computes, on the design
# of dev/fitedge-hess-mpfr.R (seed 1, n 80, 2 sin(3x) + x, noise 1e-3), at
# y scale 1 and 1e-3, and on dev/fitedge-hess-time.R's 500-position design.
# Args: arm.
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
ns <- asNamespace("frmtmb")
kap <- function(f) {
  at <- f$opt$par
  pn <- names(at)
  s <- ns$gp_marginal_setup(f$frame, max_rows = Inf)
  sig2 <- exp(2 * at[[which(pn == "betad")]])
  noise <- (sig2 / s$w + s$se2) / s$nj
  vf <- diag(as.matrix(ns$gp_marginal_cov(s, unname(at[pn == "theta"]))))
  (sum(vf) + max(noise)) / min(noise)
}
set.seed(1)
x <- round(stats::runif(80L, 0, 10), 3)
y <- 2 * sin(3 * x) + x + 1e-3 * stats::rnorm(80L)
for (ys in c(1, 1e-3)) {
  f <- suppressWarnings(frm(y ~ gp(x), data = data.frame(x = x, y = y * ys)))
  cat(sprintf("KB mpfr design, y scale %g: bound %.3e, route %s\n", ys, kap(f),
              is.function(ns$gp_marginal_outer_gr(f, f$opt$par))))
}
set.seed(20261103)
pos <- sort(runif(500, 0, 10))
d <- data.frame(x = pos, y = sin(pos) + 0.3 * rnorm(500))
f <- suppressWarnings(frm(y ~ gp(x), data = d))
cat(sprintf("KB 500 positions, noise 0.3: bound %.3e\n", kap(f)))
cat("DONE\n")
