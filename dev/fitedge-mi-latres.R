# Item 4: mi() on the response beside a latent residual, where every missing
# row is unread, and where a row is read. Args: arm (lane | base). Seeds 12301
# (the test's Gamma data) and 40401 to 40405 (lognormal against the exact
# marginal of log y on the observed rows).
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
ar1C <- function(k, ar) ar^abs(outer(seq_len(k), seq_len(k), "-")) / (1 - ar^2)
say <- function(lab, r) {
  if (inherits(r$fit, "error")) {
    cat(lab, "ERROR", substr(conditionMessage(r$fit), 1, 160), "\n")
    return(invisible(NULL))
  }
  cat(sprintf("%s ll %.8f conv %d nwarn %d\n", lab, as.numeric(logLik(r$fit)),
              r$fit$opt$convergence, length(r$warn)))
  for (w in r$warn) cat("   W:", substr(w, 1, 160), "\n")
}
# 1. lognormal: the exact marginal of log y over the observed rows of each
# group is N(mu, sigma^2 I + sderr^2 C[obs, obs])
for (seed in 40401:40405) {
  set.seed(seed)
  G <- 10L; K <- 6L
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  d <- d[order(d$g, d$t), ]
  d$x <- rnorm(nrow(d))
  e <- as.vector(replicate(G, arima.sim(list(ar = 0.6), K, sd = 0.5)))
  d$y <- exp(0.6 + 0.4 * d$x + e + rnorm(nrow(d), 0, 0.3))
  miss <- sort(sample(nrow(d), 8))
  d$y[miss] <- NA
  r <- fe_fit(bf(y | mi() ~ x + ar(t, g, cov = TRUE)), family = lognormal(),
              data = d)
  say(sprintf("LN seed %d mi()", seed), r)
  if (inherits(r$fit, "error")) next
  ex <- function(p) {
    b <- p[1:2]; sg <- exp(p[3]); s2 <- exp(2 * p[4]); ar <- tanh(p[5])
    mu <- b[1] + b[2] * d$x
    ll <- 0
    for (g in levels(d$g)) {
      i <- which(d$g == g)
      o <- !is.na(d$y[i])
      C <- s2 * ar1C(K, ar)[o, o, drop = FALSE] + diag(sg^2, sum(o))
      R <- chol(C)
      z <- backsolve(R, log(d$y[i][o]) - mu[i][o], transpose = TRUE)
      ll <- ll - sum(log(diag(R))) - 0.5 * sum(o) * log(2 * pi) -
        0.5 * sum(z^2) - sum(log(d$y[i][o]))
    }
    -ll
  }
  est <- r$fit$estimates
  p0 <- c(est$beta, est$betad, est$theta[1], atanh(0.5))
  o <- optim(p0, ex, method = "BFGS", control = list(reltol = 1e-14, maxit = 2000))
  o <- optim(o$par, ex, method = "Nelder-Mead", control = list(reltol = 1e-15, maxit = 5000))
  cat(sprintf("LN seed %d exact max ll %.8f  frm %.8f  diff %.3e\n", seed,
              -o$value, as.numeric(logLik(r$fit)),
              as.numeric(logLik(r$fit)) + o$value))
  # the same model as a weight of 0 on the missing rows
  dw <- d; dw$w <- ifelse(is.na(d$y), 0, 1); dw$y[is.na(dw$y)] <- 1
  rw <- fe_fit(bf(y | weights(w) ~ x + ar(t, g, cov = TRUE)),
               family = lognormal(), data = dw)
  say(sprintf("LN seed %d weights 0", seed), rw)
  if (!inherits(rw$fit, "error")) {
    cat(sprintf("LN seed %d beta mi %s | w0 %s\n", seed,
                paste(signif(r$fit$estimates$beta, 7), collapse = " "),
                paste(signif(rw$fit$estimates$beta, 7), collapse = " ")))
  }
}
# 2. the test's Gamma data: every missing row unread
set.seed(12301)
dg <- expand.grid(t = 1:8, g = factor(1:10))
dg <- dg[order(dg$g, dg$t), ]
dg$x <- rnorm(nrow(dg))
dg$y <- rgamma(nrow(dg), 5, 5 / exp(0.3 * dg$x))
dg$y[c(3L, 20L, 41L)] <- NA
rg <- fe_fit(bf(y | mi() ~ x + ar(t, g, cov = TRUE)),
             family = Gamma(link = "log"), data = dg)
say("GA mi()", rg)
dw <- dg; dw$w <- ifelse(is.na(dg$y), 0, 1); dw$y[is.na(dw$y)] <- 1
rgw <- fe_fit(bf(y | weights(w) ~ x + ar(t, g, cov = TRUE)),
              family = Gamma(link = "log"), data = dw)
say("GA weights 0", rgw)
if (!inherits(rg$fit, "error") && !inherits(rgw$fit, "error")) {
  cat(sprintf("GA estimates mi %s\nGA estimates w0 %s\n",
              paste(signif(unlist(rg$fit$estimates[c("beta", "betad", "theta")]), 8), collapse = " "),
              paste(signif(unlist(rgw$fit$estimates[c("beta", "betad", "theta")]), 8), collapse = " ")))
}
# 3. a read row: another response reads y's missing values
dg$w2 <- rnorm(nrow(dg), 0.5 * log(ifelse(is.na(dg$y), 1, dg$y)), 0.3)
rr <- fe_fit(bf(y | mi() ~ x + ar(t, g, cov = TRUE), family = Gamma("log")) +
               bf(w2 ~ mi(y), family = gaussian()) + set_rescor(FALSE),
             data = dg)
say("GA read by w2", rr)
cat("DONE\n")
