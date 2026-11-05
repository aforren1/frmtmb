# Which kind misses the exact marginal on a gaussian inner problem, and why.
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(167)
G <- 6L; K <- 4L
d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
e <- as.vector(replicate(G, arima.sim(list(ar = 0.6), K, sd = 1)))
eta <- 0.6 + 0.4 * d$x + e
d$pos <- exp(eta + rnorm(n, 0, 0.4))
target <- function(kind, k, ar = NA, ma = NA, rho = NA) {
  lag <- abs(outer(seq_len(k), seq_len(k), "-"))
  switch(kind,
    ar1 = ar^lag / (1 - ar^2),
    ma1 = (1 + ma^2) * (lag == 0) + ma * (lag == 1),
    arma1 = {
      g1 <- (1 + ar * ma) * (ar + ma)
      m <- g1 * ar^pmax(lag - 1, 0)
      diag(m) <- 1 + ma^2 + 2 * ar * ma
      m / (1 - ar^2)
    },
    cosy = rho + (1 - rho) * (lag == 0))
}
marg <- function(par, fit) {
  # at any outer point: sigma, beta, theta from the opt vector
  nm <- names(par)
  beta <- par[nm == "beta"]
  sig <- exp(par[nm == "betad"])
  th <- par[nm == "theta"]
  bk <- ns$frame_latres(fit$frame)[[1L]]
  cf <- ns$latres_coefs(th[-1L], bk$aux)
  kind <- bk$aux$kind
  C <- exp(2 * th[1L]) * target(kind, K, ar = cf$ar %||% NA, ma = cf$ma %||% NA,
                                rho = cf$rho %||% NA) + diag(sig^2, K)
  R <- chol(C)
  mu <- as.numeric(cbind(1, d$x) %*% beta)
  ll <- 0
  for (g in levels(d$g)) {
    r <- which(d$g == g)
    res <- backsolve(R, log(d$pos[r]) - mu[r], transpose = TRUE)
    ll <- ll - sum(log(diag(R))) - K / 2 * log(2 * pi) - 0.5 * sum(res^2) -
      sum(log(d$pos[r]))
  }
  ll
}
for (tm in c("ar(t, g, cov = TRUE)", "ma(t, g, cov = TRUE)",
             "arma(t, g, cov = TRUE)", "cosy(t, g)")) {
  fit <- suppressWarnings(frm(bf(as.formula(paste("pos ~ x +", tm))),
                              family = lognormal(), data = d))
  par <- fit$opt$par
  cat(tm, "conv", fit$opt$convergence, "logLik", format(as.numeric(logLik(fit)), digits = 12),
      "at opt", format(-fit$obj$fn(par), digits = 12),
      "marg", format(marg(par, fit), digits = 12), "\n")
  print(signif(par, 6))
  set.seed(2)
  p2 <- par + rnorm(length(par), 0, 0.2)
  cat("  away: tape", format(-fit$obj$fn(p2), digits = 12),
      "marg", format(marg(p2, fit), digits = 12), "\n")
}
