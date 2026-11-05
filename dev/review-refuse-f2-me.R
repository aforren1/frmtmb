# Item 6: unread mi() rows beside me(): the me() values follow the mi() ones in
# miss, so their indices shift. Hand deletion and fitted() on this build. Seed 11900.
source("dev/review-refuse-lib.R")
set.seed(11900); n <- 120; x <- rnorm(n); w <- rnorm(n)
d <- data.frame(x = x, y = rgamma(n, 2, 2 / exp(0.2 + 0.4 * x)), wo = w + rnorm(n, 0, 0.3), s = 0.3, g = factor(rep(1:12, 10)))
d$z <- 0.5 + 0.7 * w + rnorm(n, 0, 0.5); M <- seq(4, n, by = 12); d$y[M] <- NA
a <- frm(bf(y | mi() ~ x) + stats::Gamma(link = "log") + bf(z ~ me(wo, s) + x) + gaussian(), data = d)
b1 <- frm(bf(y ~ x), family = stats::Gamma(link = "log"), data = d[-M, ])
b2 <- frm(bf(z ~ me(wo, s) + x), data = d)
pr("unread beside me(): logLik %.6f; sum of the two separate fits %.6f; diff %.2e", as.numeric(logLik(a)), as.numeric(logLik(b1)) + as.numeric(logLik(b2)), as.numeric(logLik(a)) - as.numeric(logLik(b1)) - as.numeric(logLik(b2)))
fa <- fitted(a, resp = "z"); fb <- fitted(b2)
pr("fitted(z) against the separate fit, max %.2e; attr nobs %s", max(abs(fa[, 1] - fb[, 1])), attr(logLik(a), "nobs"))
lab <- frmtmb:::brms_miss_labels(a); pr("labels: %s", paste(head(lab, 3), collapse = " "))
