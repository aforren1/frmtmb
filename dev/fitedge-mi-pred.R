# Fix round, m3: a read mi() predictor beside a latent residual,
# `ln ~ mi(xm) + ar(t, g, cov = TRUE)` with `xm | mi() ~ 1`, against an exact
# reference with no frmtmb code. On lognormal, log ln and xm are jointly
# gaussian, so the missing xm integrate out in closed form: per group,
# log ln ~ N(b0 + b1 m, sderr^2 C + sigma^2 I + b1^2 sx^2 D), with m the
# observed xm or mx where it is missing, D the indicator of a missing row,
# and each observed xm ~ N(mx, sx). The exact maximum is found by BFGS then
# Nelder-Mead from frm()'s point. Seeds 73001 to 73005; 10 groups of 6, 6
# missing xm. Args: arm.
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
ar1C <- function(k, ar) ar^abs(outer(seq_len(k), seq_len(k), "-")) / (1 - ar^2)
for (seed in 73001:73005) {
  set.seed(seed)
  G <- 10L; K <- 6L
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$xm <- rnorm(n, 0.5, 1)
  e <- as.vector(replicate(G, arima.sim(list(ar = 0.6), K, sd = 0.5)))
  d$ln <- exp(0.6 + 0.4 * d$xm + e + rnorm(n, 0, 0.3))
  miss <- sort(sample(n, 6))
  d$xm[miss] <- NA
  r <- fe_fit(bf(ln ~ mi(xm) + ar(t, g, cov = TRUE), family = lognormal()) +
                bf(xm | mi() ~ 1, family = gaussian()) + set_rescor(FALSE),
              data = d)
  if (inherits(r$fit, "error")) {
    cat(sprintf("MP seed %d ERROR %s\n", seed, conditionMessage(r$fit)))
    next
  }
  f <- r$fit
  ex <- function(p) {
    b0 <- p[1]; b1 <- p[2]; sg <- exp(p[3]); s2 <- exp(2 * p[4])
    ar <- tanh(p[5]); mx <- p[6]; sx <- exp(p[7])
    o <- !is.na(d$xm)
    ll <- sum(dnorm(d$xm[o], mx, sx, log = TRUE))
    m <- ifelse(o, d$xm, mx)
    for (g in levels(d$g)) {
      i <- which(d$g == g)
      C <- s2 * ar1C(K, ar) + diag(sg^2 + b1^2 * sx^2 * (!o[i]), K)
      R <- chol(C)
      z <- backsolve(R, log(d$ln[i]) - b0 - b1 * m[i], transpose = TRUE)
      ll <- ll - sum(log(diag(R))) - 0.5 * K * log(2 * pi) - 0.5 * sum(z^2) -
        sum(log(d$ln[i]))
    }
    -ll
  }
  est <- f$estimates
  # beta holds ln's intercept, then mi(xm)'s slope, then xm's intercept (by
  # name below); the reference reads them by name too
  b <- est$beta
  nb <- names(b) %||% rep("", length(b))
  p0 <- c(0.6, 0.4, log(0.3), log(0.5), atanh(0.5), 0.5, 0)
  o1 <- optim(p0, ex, method = "BFGS", control = list(reltol = 1e-14, maxit = 5000))
  o1 <- optim(o1$par, ex, method = "Nelder-Mead",
              control = list(reltol = 1e-15, maxit = 20000))
  o1 <- optim(o1$par, ex, method = "BFGS", control = list(reltol = 1e-15, maxit = 5000))
  cat(sprintf("MP seed %d frm ll %.8f conv %d nwarn %d | exact max ll %.8f | diff %.3e\n",
              seed, as.numeric(logLik(f)), f$opt$convergence, length(r$warn),
              -o1$value, as.numeric(logLik(f)) + o1$value))
  cat(sprintf("MP   exact b0 %.6f b1 %.6f sigma %.6f sderr %.6f ar %.6f mx %.6f sx %.6f\n",
              o1$par[1], o1$par[2], exp(o1$par[3]), exp(o1$par[4]),
              tanh(o1$par[5]), o1$par[6], exp(o1$par[7])))
  cat("MP   frm beta", paste(signif(b, 7), collapse = " "), "| betad",
      paste(signif(exp(est$betad), 7), collapse = " "), "\n")
  for (w in r$warn) cat("   W:", substr(w, 1, 150), "\n")
}
cat("DONE\n")
