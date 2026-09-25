# The same-time tolerance of frm_ode_records(), measured. A repeat
# instant is computed as t0 + k * ii, the form both frm_ode_records()
# and frm_ode()'s own addl expansion use, while a user writes the same
# instant as a decimal literal. The error between the two, in units of
# .Machine$double.eps * max(|t0|, k * |ii|), is what the tolerance has
# to cover.
#
# Exact values come from integer arithmetic: t0 = a / 10^d and
# ii = b / 10^d with integers a, b, so the instant is (a + k b) / 10^d,
# formed with one rounding. 200000 draws, seed 20260924: d in 0..6,
# |t0| up to 2e9 (epoch seconds) and ii up to 1e5, k in 0..1000.
set.seed(20260924)
n <- 200000
d <- sample(0:6, n, TRUE)
sc <- 10^d
a <- round(stats::runif(n, -1, 1) * 10^sample(0:9, n, TRUE) * sc)
a <- pmax(pmin(a, 2e9 * sc), -2e9 * sc)
b <- round(stats::runif(n, 0, 1) * 10^sample(-2:5, n, TRUE) * sc)
k <- sample(0:1000, n, TRUE)
ok <- abs(a + k * b) < 2^53 & abs(a) < 2^53 & b > 0
a <- a[ok]; b <- b[ok]; k <- k[ok]; sc <- sc[ok]
t0 <- a / sc
ii <- b / sc
computed <- t0 + k * ii
literal <- (a + k * b) / sc
scale <- pmax(abs(t0), k * ii)
# t0 = 0 and k = 0 is the instant 0 on both sides, exactly
stopifnot(all(computed[scale == 0] == literal[scale == 0]))
keep <- scale > 0
units <- abs(computed - literal)[keep] / (.Machine$double.eps * scale[keep])
cat(sprintf("draws %d\n", length(units)))
cat(sprintf("error in eps * max(|t0|, k ii): max %.3f, 99.99%% %.3f, median %.3f\n",
            max(units), stats::quantile(units, 0.9999), stats::median(units)))
cat(sprintf("share above 1: %.4f; above 2: %.4f\n", mean(units > 1),
            mean(units > 2)))
# the chosen bound, 64 eps * max(|t0|, k |ii|), against the worst seen
cat(sprintf("margin of the 64-eps bound over the worst draw: %.1fx\n",
            64 / max(units)))
# the spacing it still separates at the scales the review named
for (s in c(1, 1e3, 1e7, 1.7e9)) {
  cat(sprintf("|t| = %-8g tolerance %.3g\n", s,
              64 * .Machine$double.eps * s))
}
