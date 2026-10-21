# Reviewer, fix round 1: windows (lo, hi] with lo in an upper tail for
# the families whose lF comes from a SUM of masses (poisson for
# 1 <= q <= 30, beta_binomial, compois, and the zi / hurdle forms on
# them): log_window_two() takes the F form there, F(hi) - F(lo), which
# needs lF to relative accuracy near 0. References by mass sums.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
F <- function(f) ns$as_frmtmb_family(f)
rl <- function(fam, y, dp, av) ns$row_lpdf(F(fam), y, y, dp, av, NULL)
lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
out <- function(lab, got, ref) {
  e <- abs(got - ref) / pmax(1, abs(ref))
  cat(sprintf("%-56s got %-14.9g ref %-14.9g rel %.1e%s\n", lab, got, ref, e,
              if (!is.finite(e) || e > 1e-8) "  <--" else ""))
}
for (w in list(c(10, 12), c(20, 25), c(25, 30), c(30, 31), c(31, 33))) {
  k <- (w[1] + 1):w[2]
  out(sprintf("poisson mu 3 interval (%g, %g]", w[1], w[2]),
      rl(poisson(), w[1], list(mu = 3), list(cens = 2, cens_y2 = w[2])), lse(dpois(k, 3, log = TRUE)))
  out(sprintf("poisson mu 3 trunc(lb = %g, ub = %g), row at lb", w[1] + 1, w[2]),
      rl(poisson(), w[1] + 1, list(mu = 3), list(trunc_lb = w[1] + 1, trunc_ub = w[2])),
      dpois(w[1] + 1, 3, log = TRUE) - lse(dpois(k, 3, log = TRUE)))
  out(sprintf("zero_inflated_poisson mu 3 interval (%g, %g]", w[1], w[2]),
      rl(zero_inflated_poisson(), w[1], list(mu = 3, zi = 0.3), list(cens = 2, cens_y2 = w[2])),
      log(0.7) + lse(dpois(k, 3, log = TRUE)))
  out(sprintf("hurdle_poisson mu 3 interval (%g, %g]", w[1], w[2]),
      rl(hurdle_poisson(), w[1], list(mu = 3, hu = 0.3), list(cens = 2, cens_y2 = w[2])),
      log(0.7) + lse(dpois(k, 3, log = TRUE)) - log1p(-exp(-3)))
  lc <- RTMB::dcompois2(0:200, 3, 1.5, log = TRUE)
  out(sprintf("compois mu 3 nu 1.5 interval (%g, %g]", w[1], w[2]),
      rl(brmsfamily("com_poisson"), w[1], list(mu = 3, nu = 1.5), list(cens = 2, cens_y2 = w[2])),
      lse(lc[k + 1]))
}
for (w in list(c(30, 35), c(40, 45), c(45, 50))) {
  k <- (w[1] + 1):w[2]
  lb <- extraDistr::dbbinom(k, 50, 0.2 * 5, 0.8 * 5, log = TRUE)
  out(sprintf("beta_binomial(50) mu .2 phi 5 interval (%g, %g]", w[1], w[2]),
      rl(beta_binomial(), w[1], list(mu = 0.2, phi = 5), list(cens = 2, cens_y2 = w[2], trials = 50)),
      lse(lb))
  out(sprintf("binomial(50) p .2 interval (%g, %g]", w[1], w[2]),
      rl(binomial(), w[1], list(mu = 0.2), list(cens = 2, cens_y2 = w[2], trials = 50)),
      lse(dbinom(k, 50, 0.2, log = TRUE)))
}
cat("== the pair's lF near 0 for 1 <= q <= 30: relative accuracy of log F\n")
for (q in c(10, 20, 30, 31)) {
  p <- ns$poisson_cdf_log(q, 3)
  cat(sprintf("q %2g: lF %.6g  ref %.6g  (lS %.6g)\n", q, p$lF, ppois(q, 3, log.p = TRUE), p$lS))
}
