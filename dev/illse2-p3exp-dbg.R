# Lane illse2, item 3: where the exp() fit's covariance error comes from
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
c0 <- 20
set.seed(1)
d <- data.frame(x = c0 + stats::runif(200))
d$y <- 2 * exp(0.8 * (d$x - c0)) + stats::rnorm(200, 0, 0.3)
f <- bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE)
fit <- frm(f, data = d, start = list(beta = c(2 * exp(-0.8 * c0), 0.8)))
p <- fit$opt$par
print(p)
xm <- m(d$x); ym <- m(d$y)
nll <- function(th) {
  r <- ym - th[1] * exp(th[2] * xm)
  s2 <- exp(2 * th[3])
  (length(d$y) * log(s2) + sum(r^2) / s2) / 2
}
He <- mp_hessian(nll, p)
Ve <- mp_inv(He)
dg <- sqrt(num(diag(Ve)))
err <- function(V) max(abs(num(mM(V) - Ve)) / outer(dg, dg))
Had <- fit$obj$he(p)
cat("AD Hessian rel err:", max(abs(Had / num(He) - 1)), "\n")
cat("solve(AD) cor err:", err(solve(Had)), "\n")
cat("vcov(full) cor err:", err(vcov(fit, full = TRUE)), "\n")
sdr <- ns$sdr_of(fit)
cat("sdreport cov.fixed cor err:", err(sdr$cov_fixed_raw %||% sdr$cov.fixed), "\n")
oh <- stats::optimHess(p, fit$obj$fn, fit$obj$gr)
cat("optimHess cor err:", err(solve(oh)), "\n")
cat("par_units:", fit$par_units, "\n")
