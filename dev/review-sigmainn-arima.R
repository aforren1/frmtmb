# Reviewer check: higher-order cov = TRUE terms against stats::arima()
# exact ML (its sigma2 is the INNOVATION variance), one long series, and
# against a dense MVN built from ARMAacf() times the psi-weight process
# variance (sum of psi^2 to lag 5000) at a fixed point.
# Usage: Rscript dev/review-sigmainn-arima.R
.libPaths(c("C:/Users/adf44/source/r/wt-sigmainn-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
set.seed(20261028)
n <- 300
cases <- list(
  ar2 = list(spec = list(ar = c(0.5, 0.3)), order = c(2, 0, 0),
             term = "ar(time = t, gr = g, p = 2, cov = TRUE)"),
  ar2neg = list(spec = list(ar = c(1.2, -0.6)), order = c(2, 0, 0),
                term = "ar(time = t, gr = g, p = 2, cov = TRUE)"),
  ma2 = list(spec = list(ma = c(0.5, -0.3)), order = c(0, 0, 2),
             term = "ma(time = t, gr = g, q = 2, cov = TRUE)"),
  arma21 = list(spec = list(ar = c(0.6, 0.2), ma = 0.4), order = c(2, 0, 1),
                term = "arma(time = t, gr = g, p = 2, q = 1, cov = TRUE)"),
  arma12 = list(spec = list(ar = 0.7, ma = c(-0.3, 0.4)), order = c(1, 0, 2),
                term = "arma(time = t, gr = g, p = 1, q = 2, cov = TRUE)"),
  ar1 = list(spec = list(ar = 0.8), order = c(1, 0, 0),
             term = "ar(time = t, gr = g, cov = TRUE)"),
  arma11 = list(spec = list(ar = 0.8, ma = -0.5), order = c(1, 0, 1),
                term = "arma(time = t, gr = g, cov = TRUE)"))
psi_var <- function(ar, ma) {
  ps <- c(1, ARMAtoMA(ar = ar, ma = ma, lag.max = 5000))
  sum(ps^2)
}
for (nm in names(cases)) {
  cs <- cases[[nm]]
  x <- rnorm(n)
  e <- as.numeric(arima.sim(cs$spec, n = n, sd = 0.7))
  d <- data.frame(t = seq_len(n), g = 1L, x = x, y = 1 + 0.5 * x + e)
  ref <- arima(d$y, order = cs$order, xreg = d$x, method = "ML",
               optim.control = list(reltol = 1e-14, maxit = 5000))
  fit <- frm(as.formula(paste("y ~ x +", cs$term)), data = d)
  ac <- fit$frame$autocor[[1]]
  th <- fit$estimates$thetaac[ac$theta_idx]
  nat <- frmtmb:::autocor_natural(th, ac)
  # internal consistency: closed form (order 1) or acv against psi weights
  phi <- nat[grepl("^ar", names(nat))]
  mth <- nat[grepl("^ma", names(nat))]
  lm_lane <- frmtmb:::autocor_log_marg(th, ac)
  lm_psi <- 0.5 * log(psi_var(phi, mth))
  cat(sprintf(paste0("%-7s sigma frm %.8f arima %.8f rel %.2e | logLik frm",
                     " %.8f arima %.8f diff %.2e | log c lane %.12f psi %.12f",
                     " diff %.1e\n"),
              nm, sigma(fit), sqrt(ref$sigma2),
              abs(sigma(fit) / sqrt(ref$sigma2) - 1),
              as.numeric(logLik(fit)), ref$loglik,
              as.numeric(logLik(fit)) - ref$loglik, lm_lane, lm_psi,
              lm_lane - lm_psi))
  cat("        coefs frm:", sprintf("%.5f", nat), " arima:",
      sprintf("%.5f", coef(ref)[seq_along(nat)]), "\n")
  # dense MVN at a fixed point away from the optimum
  par <- fit$obj$par
  par <- par + 0.05 * seq_along(par) / length(par)
  fit$obj$fn(par)
  pl <- fit$obj$env$parList(par)
  th2 <- pl$thetaac[ac$theta_idx]
  nat2 <- frmtmb:::autocor_natural(th2, ac)
  phi2 <- nat2[grepl("^ar", names(nat2))]
  mth2 <- nat2[grepl("^ma", names(nat2))]
  sg2 <- exp(pl$betad[1])
  mu2 <- as.numeric(cbind(1, d$x) %*% pl$beta)
  acf <- ARMAacf(ar = phi2, ma = mth2, lag.max = n - 1)
  S <- sg2^2 * psi_var(phi2, mth2) * toeplitz(as.numeric(acf))
  L <- chol(S)
  z <- backsolve(L, d$y - mu2, transpose = TRUE)
  ld <- -0.5 * sum(z^2) - sum(log(diag(L))) - n / 2 * log(2 * pi)
  cat(sprintf("        fixed point: -fn %.10f dense %.10f rel %.2e\n",
              -fit$obj$fn(par), ld, abs(-fit$obj$fn(par) / ld - 1)))
}
