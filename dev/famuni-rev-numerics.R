# Reviewer probe: discrete Weibull mean cost and accuracy for heavy
# tails; where a left-censored frechet / GEV row's log F underflows; the
# gradient of an outside GEV row near the exp(700) cap.
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
dwm <- frmtmb:::dw_moments

cat("\n=== A. discrete Weibull mean: direct sum against dw_moments()\n")
g <- expand.grid(mu = c(0.3, 0.9, 0.99, 0.999), shape = c(0.05, 0.1, 0.2, 0.5))
for (i in seq_len(nrow(g))) {
  mu <- g$mu[i]; k <- g$shape[i]
  t0 <- proc.time()[[3]]
  m <- dwm(mu, k)$mean
  el <- proc.time()[[3]] - t0
  # the direct sum to 1e7 terms, and how much of S is left there
  y <- seq_len(1e7)
  ds <- sum(exp(log(mu) * y^k))
  left <- exp(log(mu) * 1e7^k)
  cat(sprintf("mu %.3f shape %.2f: dw_moments %.10g (%.3f s) | 1e7-term sum %.10g, S(1e7) = %.3g | rel %.2e\n",
              mu, k, m, el, ds, left, abs(m - ds) / abs(m)))
}
cat("\ncost per row of fitted(): 1e4 rows at mu 0.9, shape 0.5\n")
t0 <- proc.time()[[3]]
invisible(dwm(rep(0.9, 1e4), rep(0.5, 1e4)))
cat("  elapsed", proc.time()[[3]] - t0, "s\n")
t0 <- proc.time()[[3]]
invisible(dwm(runif(1e4, 0.5, 0.999), runif(1e4, 0.05, 2)))
cat("  mixed parameters, elapsed", proc.time()[[3]] - t0, "s\n")
cat("extremes: shape 1e-3, mu 0.5:", dwm(0.5, 1e-3)$mean,
    "; shape 50:", dwm(0.5, 50)$mean, "; mu 1 - 1e-12, shape 0.5:",
    dwm(1 - 1e-12, 0.5)$mean, "\n")

cat("\n=== B. left censoring far in the left tail: log F\n")
fam <- frechet()
for (nu in c(2, 4, 8)) {
  s <- 1
  mu <- s * gamma(1 - 1 / nu)
  # log F = -(q/s)^-nu; the exact value is known in closed form
  for (lf in c(-100, -700, -708, -720, -745, -800)) {
    q <- s * (-lf)^(-1 / nu)
    got <- log(fam$lcdf(q, list(mu = mu, nu = nu, .eta_mu = log(mu)), list()))
    cat(sprintf("frechet nu %g: q/s = %.4g exact log F %g, family %.10g\n",
                nu, q / s, lf, got))
  }
}
famg <- gen_extreme_value()
for (xi in c(0.3, 0, -0.3)) {
  for (lf in c(-100, -700, -720, -745, -800)) {
    # log F = -t^(-1/xi), t = (-lf)^(-xi); x = (t - 1) / xi
    x <- if (xi == 0) -log(-lf) else ((-lf)^(-xi) - 1) / xi
    got <- log(famg$lcdf(x, list(mu = 0, sigma = 1, xi = xi, .eta_sigma = 0),
                         list()))
    cat(sprintf("GEV xi %g: x = %.4g exact log F %g, family %.10g\n",
                xi, x, lf, got))
  }
}
cat("\nA realistic fit: frechet nu = 4, 10 percent left censored at a limit\n")
set.seed(1101)
n <- 300
y <- 2 / gamma(1 - 1 / 4) * rexp(n)^(-1 / 4)
lim <- quantile(y, 0.1)
d <- data.frame(y = pmax(y, lim), cc = ifelse(y < lim, -1, 0))
f <- frm(y | cens(cc) ~ 1, family = frechet(), data = d)
cat("conv", f$opt$convergence, "logLik", as.numeric(logLik(f)),
    " limit / scale at the fit:",
    lim / (exp(fixef(f)[1, 1]) / gamma(1 - 1 / summary(f)$spec_pars["nu", 1])),
    "\n")
# the start, where the family's mu = 1.2 median, nu = 3
s0 <- 1.2 * median(d$y) / gamma(1 - 1 / 3)
cat("log F of the limit at the start:", -(lim / s0)^(-3), "\n")

cat("\n=== C. gradient of an outside GEV row near the exp(700) cap (xi 0.06, lower end at x = -16.7)\n")
lp_row <- function(p, y) {
  famg$lpdf(y, list(mu = p[1], sigma = exp(p[2]), xi = expm1(p[3]),
                    .eta_sigma = p[2], .eta_xi = p[3]), list())
}
for (sg in c(1, 0.01, 1e-3)) {
  for (xx in c(-19, -20.2, -20.3, -21)) {
    p <- c(0, log(sg), log1p(0.06))
    y <- xx * sg
    tp <- RTMB::MakeTape(function(p) lp_row(p, y), p)
    v <- tp(p)
    gr <- tp$jacobian(p)
    cat(sprintf("sigma %g x %g: lpdf %.4g, gradient %s\n", sg, xx, v,
                paste(format(gr, digits = 4), collapse = " ")))
  }
}
