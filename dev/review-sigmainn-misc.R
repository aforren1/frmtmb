# Reviewer check, frmtmb alone (lane build or rellib-r14):
#  mv     mvbf(y1 ~ x + ar(cov = TRUE), y2 ~ x, rescor = FALSE) against
#         the two univariate fits: logLik additive, sigma and pearson
#         residuals of y1 equal the univariate ones
#  resx   RE + ma(cov = TRUE) + sigma ~ x2: pearson residual of each row
#         is (y - fitted) / (sigma_i * c), c from the dense ACF
#  ar2sx  ar(p = 2, cov = TRUE) + sigma ~ x2: logLik against a dense MVN
#         with brms's construction diag(sigma_i) Gamma diag(sigma_i),
#         Gamma the unit-innovation autocovariance (ARMAacf times the
#         psi variance), ragged groups
#  dev    deviance residuals under cov = TRUE
# Usage: Rscript dev/review-sigmainn-misc.R after|before
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-sigmainn-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(20261031)
mk <- function(ng, nt, spec) {
  d <- do.call(rbind, lapply(seq_len(ng), function(i) {
    k <- sample(4:nt, 1)
    data.frame(g = i, t = seq_len(k),
               e = as.numeric(arima.sim(spec, k, sd = 1)))
  }))
  d$g <- factor(d$g)
  d$x <- rnorm(nrow(d)); d$x2 <- rnorm(nrow(d))
  d$h <- factor(sample(1:5, nrow(d), TRUE))
  d$y <- 1 + 0.5 * d$x + exp(-0.5 + 0.4 * d$x2) * d$e + c(-.5, -.2, 0, .3, .4)[d$h]
  d$y1 <- d$y
  d$y2 <- 2 - 0.3 * d$x + rnorm(nrow(d))
  d
}
psi_var <- function(ar, ma) sum(c(1, ARMAtoMA(ar = ar, ma = ma, lag.max = 5000))^2)
dense_ll <- function(d, mu, sg, ar, ma) {
  tot <- 0
  for (gi in split(seq_len(nrow(d)), d$g)) {
    k <- length(gi)
    G <- psi_var(ar, ma) * toeplitz(as.numeric(ARMAacf(ar = ar, ma = ma, lag.max = k - 1)))
    S <- diag(sg[gi], k) %*% G %*% diag(sg[gi], k)
    L <- chol(S)
    z <- backsolve(L, d$y[gi] - mu[gi], transpose = TRUE)
    tot <- tot - 0.5 * sum(z^2) - sum(log(diag(L))) - k / 2 * log(2 * pi)
  }
  tot
}
# mv
d <- mk(25, 9, list(ar = 0.7))
fmv <- frm(mvbf(bf(y1 ~ x + ar(time = t, gr = g, cov = TRUE)), bf(y2 ~ x),
                rescor = FALSE), data = d)
f1 <- frm(bf(y1 ~ x + ar(time = t, gr = g, cov = TRUE)), data = d)
f2 <- frm(bf(y2 ~ x), data = d)
cat(sprintf("mv: logLik %.10f sum of univariate %.10f rel %.2e\n",
            as.numeric(logLik(fmv)), as.numeric(logLik(f1) + logLik(f2)),
            abs(as.numeric(logLik(fmv)) / as.numeric(logLik(f1) + logLik(f2)) - 1)))
cat("mv: sigma", sprintf("%.8f", sigma(fmv)), "| univariate", sprintf("%.8f", c(sigma(f1), sigma(f2))), "\n")
r_mv <- as.numeric(residuals(fmv, type = "pearson", resp = "y1"))
r_1 <- residuals(f1, type = "pearson")
cat(sprintf("mv: pearson y1 max rel diff to univariate %.2e (sd %.5f)\n",
            max(abs(r_mv - as.numeric(r_1))) / max(abs(r_1)), sd(r_1)))
# resx
d2 <- mk(25, 9, list(ma = 0.6))
fr <- frm(bf(y ~ x + (1 | h) + ma(time = t, gr = g, cov = TRUE), sigma ~ x2), data = d2)
ac <- fr$frame$autocor[[1]]
ma <- frmtmb:::autocor_natural(fr$estimates$thetaac[ac$theta_idx], ac)[[1]]
sgi <- as.numeric(frm_linpred(fr, dpar = "sigma", type = "response"))
rp <- residuals(fr, type = "pearson")
rr <- residuals(fr, type = "response")
cat(sprintf("resx: ma %.5f; max |pearson - response/(sigma_i sqrt(1+ma^2))| rel %.2e; vs response/sigma_i rel %.2e\n",
            ma, max(abs(rp - rr / (sgi * sqrt(1 + ma^2)))) / max(abs(rp)),
            max(abs(rp - rr / sgi)) / max(abs(rp))))
rd <- residuals(fr, type = "deviance")
cat(sprintf("dev: deviance residuals: %d NA of %d; max |dev - pearson| %.2e\n",
            sum(is.na(rd)), length(rd), max(abs(rd - rp), na.rm = TRUE)))
# ar2sx at the optimum and at a moved point
d3 <- mk(30, 9, list(ar = c(0.5, 0.25)))
fa <- frm(bf(y ~ x + ar(time = t, gr = g, p = 2, cov = TRUE), sigma ~ x2), data = d3)
for (shift in c(0, 0.07)) {
  par <- fa$obj$par + shift * seq_along(fa$obj$par) / length(fa$obj$par)
  nll <- fa$obj$fn(par)
  pl <- fa$obj$env$parList(par)
  ac <- fa$frame$autocor[[1]]
  arv <- frmtmb:::autocor_natural(pl$thetaac[ac$theta_idx], ac)
  X <- cbind(1, d3$x)
  mu <- as.numeric(X %*% pl$beta)
  sg <- exp(as.numeric(cbind(1, d3$x2) %*% pl$betad))
  ref <- dense_ll(d3, mu, sg, arv, numeric(0))
  cat(sprintf("ar2sx shift %.2f: -fn %.10f dense(brms construction) %.10f rel %.2e\n",
              shift, -nll, ref, abs(-nll / ref - 1)))
}
