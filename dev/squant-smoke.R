# Lane squant: a first look at quantile_elf() on the lane build.
# Run: Rscript dev/squant-smoke.R > dev/squant-out/smoke.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb")

# d/p/q/r on doubles
lam <- 0.25; p <- 0.9; mu <- 1; s <- 0.7
y <- seq(-5, 5, by = 0.5)
lp <- ns$elf_lpdf(y, mu, s, p, 1 - p, lam)
cp <- ns$elf_cdf_log(y, mu, s, p, 1 - p, lam)
num <- vapply(y, function(q) integrate(function(z) {
  exp(ns$elf_lpdf(z, mu, s, p, 1 - p, lam))
}, -Inf, q, rel.tol = 1e-12)$value, 0)
print(cbind(y, lp, F = exp(cp$lF), Fnum = num, S = exp(cp$lS)))
qq <- ns$elf_quantile(c(1e-8, 0.1, 0.5, 0.9, 1 - 1e-8), mu, s, p, lam)
print(qq)
print(exp(ns$elf_cdf_log(qq, mu, s, p, 1 - p, lam)$lF))
set.seed(1)
dr <- ns$elf_draw(1e5, mu, s, p, 1 - p, lam)
cat("mean draw", mean(dr), "mean fn",
    mu + lam * s * (digamma(lam * (1 - p)) - digamma(lam * p)), "\n")
cat("P(draw < mu)", mean(dr < mu), "\n")

# a fit against rq
set.seed(2)
n <- 400
dd <- data.frame(x = runif(n))
dd$y <- 1 + 2 * dd$x + (0.5 + dd$x) * rnorm(n)
for (h in c(1, 0.25, 0.05, 0.01)) {
  fit <- frm(bf(y ~ x, quantile = 0.9), family = quantile_elf(smoothing = h),
             data = dd)
  cat("h", h, "conv", fit$opt$convergence, "\n")
  print(fixef(fit))
}
rq <- quantreg::rq(y ~ x, tau = 0.9, data = dd)
print(coef(rq))
fit <- frm(bf(y ~ x, quantile = 0.9), family = quantile_elf(), data = dd)
print(summary(fit))
dd$row <- seq_len(n)
print(sqrt(diag(vcov(fit))))
print(sqrt(diag(vcov_cluster(fit, ~ row, "CR0"))))
print(summary(rq, se = "nid"))
