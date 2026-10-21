# Reviewer: lgamma_shift_diff() against a 256-bit Rmpfr reference over
# the shifts a dirichlet-multinomial count can take (s to 1e6, beyond
# the lane's 2000), next to the plain lgamma difference brms writes.
# Then the whole density: frmtmb against Rmpfr at large T, where the
# shift form's loop of log1p(s / (a + j)) carries counts of 1e5.
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R")
suppressMessages(library(Rmpfr))
sd_ <- frmtmb:::lgamma_shift_diff
ref <- function(a, s) {
  A <- mpfr(a, 256)
  as.numeric(lgamma(A + s) - lgamma(A))
}
g <- expand.grid(a = 10^c(-10, -6, -3, -1, 0, 0.5, 1, 2, 3, 6, 9, 12, 15),
                 s = c(1, 3, 17, 99, 2000, 1e4, 1e5, 1e6))
g$ref <- mapply(ref, g$a, g$s)
g$shift <- mapply(sd_, g$a, g$s)
g$plain <- lgamma(g$a + g$s) - lgamma(g$a)
g$rel_shift <- abs(g$shift / g$ref - 1)
g$rel_plain <- abs(g$plain / g$ref - 1)
cat("max relative error by a:\n")
print(aggregate(cbind(rel_shift, rel_plain) ~ a, data = g, FUN = max),
      digits = 3)
cat("max relative error by s:\n")
print(aggregate(cbind(rel_shift, rel_plain) ~ s, data = g, FUN = max),
      digits = 3)
cat(sprintf("overall: shift %.3g at (a %g, s %g); plain %.3g at (a %g, s %g)\n",
            max(g$rel_shift), g$a[which.max(g$rel_shift)],
            g$s[which.max(g$rel_shift)], max(g$rel_plain),
            g$a[which.max(g$rel_plain)], g$s[which.max(g$rel_plain)]))
cat("points where shift is worse than plain by > 4x and > 1e-14:",
    sum(g$rel_shift > 4 * g$rel_plain & g$rel_shift > 1e-14), "of",
    nrow(g), "\n")
print(g[g$rel_shift > 4 * g$rel_plain & g$rel_shift > 1e-14,
        c("a", "s", "rel_shift", "rel_plain")], digits = 3)

# the density at microbiome-like totals: K = 6, T up to 2e5
dm <- dirichlet_multinomial(K = 6)
dmref <- function(y, alpha) {
  A <- mpfr(alpha, 256)
  ph <- sum(A)
  tot <- sum(y)
  as.numeric(lgamma(ph) - lgamma(ph + tot) + lgamma(mpfr(tot + 1, 256)) -
               sum(lgamma(mpfr(y + 1, 256))) +
               sum(lgamma(A + y) - lgamma(A)))
}
set.seed(11)
worst <- c(frm = 0, brms = 0)
cases <- 0
for (tot in c(1e3, 3e4, 2e5)) {
  for (phi in c(0.05, 2, 60, 1e4, 1e8)) {
    for (r in 1:3) {
      eta <- rnorm(5, 0, 2)
      p <- exp(c(0, eta)) / sum(exp(c(0, eta)))
      y <- as.vector(rmultinom(1, tot, p))
      ym <- matrix(y, 1)
      dp <- c(stats::setNames(as.list(eta), dm$primary_dpars),
              list(phi = phi))
      v <- as.numeric(frmtmb:::row_lpdf(dm, ym, ym, dp, list(trials = tot),
                                         NULL))
      a <- phi * p
      b <- lgamma(phi) + lgamma(tot + 1) - lgamma(tot + phi) +
        sum(lgamma(y + a) - lgamma(a) - lgamma(y + 1))
      e <- dmref(y, a)
      worst <- pmax(worst, c(abs(v - e), abs(b - e)) /
                      (abs(e) * .Machine$double.eps))
      cases <- cases + 1
    }
  }
}
cat("density at T in {1e3, 3e4, 2e5}, phi 0.05..1e8,", cases,
    "cases: worst error in ulps of |log P|: frmtmb", signif(worst[1], 4),
    " brms formula", signif(worst[2], 4), "\n")
