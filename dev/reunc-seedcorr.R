# Lane wt-reunc: predict() takes the group-effect draw from one seed
# per replicate and the simulation from another. If two streams seeded
# by two independent integers were correlated, the predictive variance
# would come out below Var(b) + sigma^2, which is what
# dev/reunc-log/plugvar.txt shows at about 1 percent. This measures the
# correlation directly, in base R, with no frmtmb on the path.
N <- 200000L
set.seed(99)
a <- sample.int(.Machine$integer.max, N)
b <- sample.int(.Machine$integer.max, N)
za <- numeric(N); zb <- numeric(N)
for (i in seq_len(N)) {
  set.seed(a[i]); za[i] <- stats::rnorm(1)
  set.seed(b[i]); zb[i] <- stats::rnorm(1)
}
r <- stats::cor(za, zb)
cat(sprintf("N %d  cor %.5f  se %.5f  z %.2f\n", N, r, 1 / sqrt(N),
            r * sqrt(N)))
cat(sprintf("var(za) %.5f var(zb) %.5f var(za + 0.4 * zb) %.5f (want %.5f)\n",
            var(za), var(zb), var(za + 0.4 * zb), 1 + 0.16))
