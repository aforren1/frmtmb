# me() against the closed-form marginal likelihood (lane me, 2026-09-25).
#
# With a gaussian response and linear me() terms, (y, x_obs, ...) is
# multivariate normal once the latent values are integrated out, so the
# Laplace approximation is exact and logLik() must equal the
# mvtnorm density at the fitted parameters. The same constructions as
# tests/testthat/test-me.R, printed with the numbers. For the one-term
# model an independent optimizer (optim BFGS on the closed form) is run
# from a perturbed start to compare the maxima.
#
#   Rscript dev/me-closed-form.R > dev/me-closed-form.txt
.libPaths(c("/opt/rlib/lane-me", "/opt/rlib/base", "/opt/rlib/deps",
            "/opt/r/lib/R/library"))
suppressMessages(library(frmtmb))

me_data <- function(seed = 31, n = 90) {
  set.seed(seed)
  tx <- rnorm(n, 1, 0.8)
  tz <- 0.5 * tx + rnorm(n, 0, 0.7)
  sx <- runif(n, 0.2, 0.5)
  d <- data.frame(x = tx + rnorm(n, 0, sx), sx = sx,
                  z = tz + rnorm(n, 0, 0.3), sz = 0.3, w = rnorm(n))
  d$y <- 2 + 0.7 * tx - 0.4 * tz + 0.2 * d$w + rnorm(n, 0, 0.5)
  d
}
me_closed_ll <- function(d, b, bw, mu, sdv, R, sigma, me_vars) {
  M <- length(me_vars)
  S <- diag(sdv, M) %*% R %*% diag(sdv, M)
  out <- 0
  for (i in seq_len(nrow(d))) {
    noise <- vapply(me_vars, function(v) d[[paste0("s", v)]][i], 0)
    m <- c(b[1] + sum(b[-1] * mu) + bw * d$w[i], mu)
    Sb <- S %*% b[-1]
    V <- rbind(c(drop(t(b[-1]) %*% Sb) + sigma^2, Sb),
               cbind(Sb, S + diag(noise^2, M)))
    obs <- c(d$y[i], vapply(me_vars, function(v) d[[v]][i], 0))
    out <- out + mvtnorm::dmvnorm(obs, m, V, log = TRUE)
  }
  out
}
sig <- function(fit) exp(fit$estimates[["betad"]][["sigma_(Intercept)"]])
report <- function(lab, fit, ll) {
  a <- as.numeric(logLik(fit))
  cat(sprintf("%-22s logLik %.12f closed form %.12f diff %.3e rel %.3e\n",
              lab, a, ll, a - ll, abs(a - ll) / abs(ll)))
}

d <- me_data()
fit <- frm(bf(y ~ me(x, sx) + w) + gaussian(), data = d)
fe <- fixef(fit)[, "Estimate"]
report("one term", fit,
       me_closed_ll(d, fe[c("Intercept", "mexsx")], fe[["w"]],
                    fit$estimates$meanme[[1]],
                    exp(fit$estimates$logsdme[[1]]), diag(1), sig(fit),
                    "x"))
nll <- function(p) {
  -me_closed_ll(d, p[1:2], p[3], p[4], exp(p[5]), diag(1), exp(p[6]), "x")
}
p0 <- c(fe[["Intercept"]], fe[["mexsx"]], fe[["w"]],
        fit$estimates$meanme[[1]], fit$estimates$logsdme[[1]],
        log(sig(fit)))
op <- stats::optim(p0 + 0.05, nll, method = "BFGS",
                   control = list(reltol = 1e-14, maxit = 2000))
se <- fixef(fit)[c("Intercept", "mexsx", "w"), "Est.Error"]
cat(sprintf(paste0("one term, optim max %.12f frmtmb max %.12f diff %.3e; ",
                   "max |coef diff| / SE %.3e\n"),
            -op$value, as.numeric(logLik(fit)),
            -op$value - as.numeric(logLik(fit)),
            max(abs(op$par[1:3] - fe[c("Intercept", "mexsx", "w")]) / se)))

d <- me_data(seed = 32)
fit <- frm(bf(y ~ me(x, sx) + me(z, sz) + w) + gaussian(), data = d)
fe <- fixef(fit)[, "Estimate"]
hy <- summary(fit)$me
r <- hy["corme__mex__mez", "Estimate"]
report("two terms, corme", fit,
       me_closed_ll(d, fe[c("Intercept", "mexsx", "mezsz")], fe[["w"]],
                    hy[c("meanme_mex", "meanme_mez"), 1],
                    hy[c("sdme_mex", "sdme_mez"), 1],
                    matrix(c(1, r, r, 1), 2), sig(fit), c("x", "z")))
f0 <- frm(bf(y ~ me(x, sx) + me(z, sz) + w) + set_mecor(FALSE) +
            gaussian(), data = d)
fe <- fixef(f0)[, "Estimate"]
hy <- summary(f0)$me
report("two terms, mecor FALSE", f0,
       me_closed_ll(d, fe[c("Intercept", "mexsx", "mezsz")], fe[["w"]],
                    hy[c("meanme_mex", "meanme_mez"), 1],
                    hy[c("sdme_mex", "sdme_mez"), 1], diag(2), sig(f0),
                    c("x", "z")))

set.seed(33)
ng <- 25
g <- factor(sprintf("s%02d", rep(seq_len(ng), each = 4)))
tx <- rnorm(ng, 0.5, 1)
sxg <- runif(ng, 0.2, 0.6)
xg <- tx + rnorm(ng, 0, sxg)
d <- data.frame(g = g, xg = xg[g], sxg = sxg[g])
d$y <- 1 + 0.8 * tx[g] + rnorm(nrow(d), 0, 0.4)
fit <- frm(bf(y ~ me(xg, sxg, gr = g)) + gaussian(), data = d)
fe <- fixef(fit)[, "Estimate"]
b0 <- fe[["Intercept"]]
b1 <- fe[["mexgsxggrEQg"]]
mu <- fit$estimates$meanme[[1]]
s <- exp(fit$estimates$logsdme[[1]])
sg <- sig(fit)
ll <- 0
for (j in seq_len(ng)) {
  rows <- which(as.integer(g) == j)
  k <- length(rows)
  V <- rbind(cbind(b1^2 * s^2 + diag(sg^2, k), b1 * s^2),
             c(rep(b1 * s^2, k), s^2 + sxg[j]^2))
  ll <- ll + mvtnorm::dmvnorm(c(d$y[rows], xg[j]),
                              c(rep(b0 + b1 * mu, k), mu), V, log = TRUE)
}
report("gr = g, 25 levels", fit, ll)
