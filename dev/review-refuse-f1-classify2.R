# Short check, item 1: autocor reads, the plain column on both builds, the chain
# on a family where the Laplace error is large. Seed 11001.
source("dev/review-refuse-lib.R")
pr("build %s", dirname(find.package("frmtmb")))
set.seed(11001)
n <- 120; x <- rnorm(n)
d <- data.frame(x = x, g = rep(1:12, each = 10), t0 = rep(1:10, 12), yn = arima.sim(list(ar = 0.5), n) + 0.5 * x, w = rnorm(n), y = rgamma(n, 2, 1))
M <- seq(4, n, by = 12); d$yn[M] <- NA
d$z <- 0.5 + 0.6 * ifelse(is.na(d$yn), 0, d$yn) + rnorm(n)
a <- tr(frm(bf(yn | mi() ~ x + ar(time = t0, gr = g, p = 1)), data = d))
pr("ar(cov = FALSE) with mi: %s; latent %s of 10", msg(a, 130), if (is_err(a)) "-" else length(a$estimates$miss))
b <- tr(brms::make_stancode(brms::bf(yn | mi() ~ x + ar(time = t0, gr = g, p = 1)), data = d)); pr("  brms: %s", msg(b, 100))
d2 <- d; d2$yy <- d$yn
for (fm in list(bf(z ~ yy + x), bf(z | se(abs(yy) + 0.1) ~ x), bf(z | weights(abs(yy) + 0.1) ~ x))) {
  a <- tr(frm(bf(yy | mi() ~ x) + gaussian() + fm + gaussian(), data = d2))
  pr("%-34s: %s", deparse1(fm$formula), msg(a, 160))
}
a <- tr(frm(bf(z ~ yy + x), data = d2)); pr("alone, z ~ yy + x: %s; nobs %s", msg(a, 80), if (is_err(a)) "-" else nobs(a))
