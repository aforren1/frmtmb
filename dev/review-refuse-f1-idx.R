# Short check, item 1: a reader through idx whose rows miss every missing value:
# no latent is left, and the reader still names mi(y). Also a read row whose
# response has a mix, with fitted() and predict() on the reader. Seed 11600.
source("dev/review-refuse-lib.R")
set.seed(11600)
n <- 120; x <- rnorm(n)
y0 <- rgamma(n, 2, 2 / exp(0.2 + 0.4 * x))
d <- data.frame(x = x, y = y0, id = seq_len(n)); d$z <- 0.5 + 0.6 * y0 + rnorm(n, 0, 0.7)
M <- seq(4, n, by = 12); d$y[M] <- NA
Gam <- stats::Gamma(link = "log")
for (k in list(c(all = 10L), c(some = 6L))) {
  dd <- d; dd$ok <- TRUE; dd$ok[M[seq_len(k)]] <- FALSE; dd$z2 <- ifelse(dd$ok, dd$z, 0)
  a <- tr(frm(bf(y | mi() + index(id) ~ x) + Gam + bf(z2 | subset(ok) ~ mi(y, idx = id) + x) + gaussian(), data = dd))
  pr("reader skips %d of 10: %s; latent %s", k, msg(a, 120), if (is_err(a)) "-" else length(a$estimates$miss))
  if (is_err(a)) next
  d3 <- dd[dd$ok, ]
  r <- tr(frm(bf(y | mi() ~ x) + Gam + bf(z2 ~ mi(y) + x) + gaussian(), data = d3))
  if (!is_err(r)) pr("  minus hand deletion %.2e", as.numeric(logLik(a)) - as.numeric(logLik(r)))
  for (me in c("fitted", "predict", "residuals")) {
    v <- tr(do.call(me, list(a, resp = "z2"))); pr("  %s(resp = z2): %s", me, msg(v, 130))
    v <- tr(do.call(me, list(a, resp = "y"))); pr("  %s(resp = y): %s", me, msg(v, 130))
  }
  v <- tr(simulate(a, nsim = 1, seed = 1)); pr("  simulate: %s", msg(v, 130))
  v <- tr(fitted(a, newdata = d3[1:4, ], resp = "z2")); pr("  fitted(newdata, z2): %s", msg(v, 130))
}
