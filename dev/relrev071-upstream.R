# Release reviewer of 0.71.0: rerun the repros of dev/upstream-bugs.md
# brms-24, brms-26 and brms-27 on brms 2.23.0.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
cat("brms", as.character(packageVersion("brms")), "\n")
b <- asNamespace("brms")
# brms-24: R-side pmf mean at mu 3, shape 2, against the Stan-side pmf
x <- 0:200
pr <- b$dcom_poisson(x, mu = 3, shape = 2)
cat("brms-24 R side: mass", sum(pr), "mean", sum(x * pr) / sum(pr), "\n")
ls <- x * log(3) - 2 * lgamma(x + 1)
ps <- exp(ls - max(ls)); ps <- ps / sum(ps)
cat("brms-24 Stan side (mu^y / (y!)^nu): mean", sum(x * ps), "\n")
# brms-26
set.seed(1)
d <- data.frame(f = factor(sample(c("a", "b", "c"), 60, TRUE)),
                x = rnorm(60), z = rnorm(60))
r <- tryCatch({ brms::stancode(brms::bf(f ~ x, mub ~ z), data = d,
                               family = brms::categorical()); "ok" },
              error = function(e) conditionMessage(e))
cat("brms-26:", substr(r, 1, 160), "\n")
# brms-27: posterior_epred_wiener formula against z = bias * bs
f <- function(mu, bs, bias, ndt) {
  ndt - bias / mu + bs / mu * (exp(-2 * mu * bias) - 1) /
    (exp(-2 * mu * bs) - 1)
}
g <- function(mu, bs, bias, ndt) {
  z <- bias * bs
  ndt - z / mu + bs / mu * (exp(-2 * mu * z) - 1) / (exp(-2 * mu * bs) - 1)
}
cat("brms-27 at bs 2.5, mu 1.5, bias 0.7, ndt 0: brms formula", f(1.5, 2.5,
    0.7, 0), "with z = bias * bs", g(1.5, 2.5, 0.7, 0), "\n")
src <- deparse(b$posterior_epred_wiener)
cat("brms-27 source:", paste(trimws(src), collapse = " "), "\n")
