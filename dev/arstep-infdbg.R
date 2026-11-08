# Lane arstep, fix round: the one Inf of dev/arstep-armarg.R (step 160 2):
# which random point, and which piece fails. Same points (seed 9100).
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
n <- 160L; seed <- 2L
set.seed(74000 + seed * 31 + n)
x <- sort(round(runif(n, 0, 10), 4))
d <- data.frame(x = x, t = seq_len(n), y = 2 * (x > 5) - 1 + as.numeric(arima.sim(list(ar = 0.4), n, sd = 0.03)))
fit <- frm(y ~ gp(x) + ar(time = t), data = d, family = gaussian())
s <- ns$gp_marginal_setup(fit$frame)
set.seed(9100)
p0 <- fit$opt$par
for (k in 1:20) {
  p <- p0 + rnorm(length(p0), sd = 0.3)
  pl <- fit$obj$env$parList(p)
  v <- ns$gp_ar_nll(s, list(beta = pl$beta, betad = pl$betad, theta = pl$theta, thetaac = pl$thetaac))
  if (!is.finite(v)) {
    V <- ns$gp_marginal_cov(s, pl$theta)
    pr <- ns$gp_ar_proj(s, pl$thetaac)
    a <- s$ar
    M <- t(ns$gp_ar_left(t(ns$gp_ar_left(V, a, pr$phi)), a, pr$phi))
    M <- 0.5 * (M + t(M)); sig2 <- exp(2 * pl$betad)
    ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
    evV <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    cat(sprintf("INF k %d sdgp %.4g lscale %.4g sigma %.4g ar %.4g | eig M min %.3e max %.3e | eig V min %.3e max %.3e | obj %.4f\n",
                k, exp(pl$theta[1]), exp(pl$theta[2]), sqrt(sig2), pl$thetaac, min(ev), max(ev), min(evV), max(evV), fit$obj$fn(p)))
  }
}
cat("DONE\n")
