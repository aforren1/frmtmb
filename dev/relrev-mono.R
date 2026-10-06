# Reviewer: is each ordinal CDF monotone in floating point near its
# branch points? A non-monotone step would give a cumulative() row
# without cs() a negative category probability, which ord_sim() now
# draws as NA. 2e5 adjacent doubles around each point.
chk <- function(f, a, nm) {
  x <- a + (seq_len(200000) - 100000) * abs(a) * .Machine$double.eps
  if (a == 0) x <- (seq_len(200000) - 100000) * 1e-300
  dd <- diff(f(x))
  cat(sprintf("%-8s at %-12g negative steps %d, min step %g\n", nm, a, sum(dd < 0), min(dd)))
}
for (a in c(-38, -37.5193, -5.656854, -0.67448975, 0, 0.67448975, 5.656854, 8.2924)) chk(pnorm, a, "pnorm")
for (a in c(-36, -18, 0, 18, 33)) chk(plogis, a, "plogis")
for (a in c(-1e8, -1, 0, 1, 1e8)) chk(pcauchy, a, "pcauchy")
cll <- function(x) 1 - exp(-exp(x))
for (a in c(-36, -1, 0, 1, 3.6)) chk(cll, a, "cloglog")
