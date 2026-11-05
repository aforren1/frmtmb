# Item 2: every method on an mi() fit, per family. A placeholder must never
# read as data. Seed 9500. Hand references: fitted against the family's mean,
# residuals NA at the missing rows, simulate's row, predict at newdata.
source("dev/review-refuse-lib.R")
set.seed(9500)
n <- 150
x <- rnorm(n)
eta <- 0.3 + 0.5 * x
sets <- list(
  gaussian = list(f = gaussian(), y = rnorm(n, eta, 1), inv = identity),
  student = list(f = student(), y = eta + rt(n, 5), inv = identity),
  skew_normal = list(f = skew_normal(), y = eta + abs(rnorm(n)) + rnorm(n, 0, 0.5), inv = identity),
  beta = list(f = Beta(), y = rbeta(n, plogis(eta) * 5, (1 - plogis(eta)) * 5), inv = plogis),
  Gamma = list(f = stats::Gamma(link = "log"), y = rgamma(n, 2, 2 / exp(eta)), inv = exp),
  weibull = list(f = weibull(), y = rweibull(n, 2, exp(eta)), inv = exp),
  exponential = list(f = exponential(), y = rexp(n, 1 / exp(eta)), inv = exp),
  frechet = list(f = frechet(), y = 1 / rweibull(n, 3, exp(-eta)), inv = exp),
  inverse.gaussian = list(f = stats::inverse.gaussian(link = "log"), y = rgamma(n, 4, 4 / exp(eta)), inv = exp),
  lognormal = list(f = lognormal(), y = rlnorm(n, eta, 0.6), inv = identity),
  exgaussian = list(f = exgaussian(), y = eta + rnorm(n, 0, 0.5) + rexp(n, 1), inv = identity),
  gen_extreme_value = list(f = gen_extreme_value(), y = eta + 0.5 * (rexp(n)^(-0.1) - 1) / 0.1, inv = identity),
  von_mises = list(f = von_mises(), y = ((2 * atan(eta) + rnorm(n, 0, 0.5) + pi) %% (2 * pi)) - pi, inv = identity))
mrows <- seq(4, 148, by = 12)
ok_na <- function(v) all(is.na(v[mrows])) && all(!is.na(v[-mrows]))
for (nm in names(sets)) {
  s <- sets[[nm]]
  d <- data.frame(x = x, y = s$y); d$y[mrows] <- NA
  fit <- tr(frm(bf(y | mi() ~ x), family = s$f, data = d))
  if (is_err(fit)) { pr("%-18s FIT ERROR %s", nm, msg(fit)); next }
  cc <- frm(bf(y ~ x), family = s$f, data = d[-mrows, ])
  res <- list()
  for (ty in c("response", "pearson", "deviance")) {
    r <- tr(residuals(fit, type = ty))
    res[[ty]] <- if (is_err(r)) msg(r, 50) else {
      v <- if (is.matrix(r)) r[, 1] else r
      paste0(ok_na(v), "/", length(v))
    }
  }
  ft <- tr(fitted(fit))
  fv <- if (!is_err(ft)) ft[, 1]
  # the fitted mean at a missing row reads the same parameters as any row
  pnew <- tr(predict(fit, newdata = d[mrows, c("x", "y"), drop = FALSE], re_formula = NA))
  pnew2 <- tr(fitted(fit, newdata = d[mrows, "x", drop = FALSE]))
  dpn <- if (!is_err(pnew2) && !is.null(fv)) max(abs(pnew2[, 1] - fv[mrows])) else NA
  sims <- tr(simulate(fit, nsim = 2, seed = 1))
  simd <- if (is_err(sims)) msg(sims, 60) else {
    v <- as.matrix(sims)[, 1]
    sprintf("NAatmiss=%s", sum(is.na(v[mrows])))
  }
  rg <- tr(c(nobs(fit), nobs(cc)))
  ll <- tr(logLik(fit))
  sm <- tr(summary(fit))
  vc <- tr(vcov(fit))
  ci <- tr(confint(fit))
  an <- tr(AIC(fit))
  up <- tr(update(fit, . ~ . - x))
  ll_up <- if (is_err(up)) msg(up, 60) else sprintf("%.2f", as.numeric(logLik(up)))
  mf <- tr(model.frame(fit))
  cook <- tr(cooks.distance(fit))
  ckd <- if (is_err(cook)) msg(cook, 50) else paste0(sum(is.na(cook[mrows])), "NA")
  hv <- tr(hatvalues(fit))
  hvd <- if (is_err(hv)) msg(hv, 50) else "ok"
  pr("%-18s resid %s fitted %s newdata-vs-fit %.1e sim %s nobs %s logLik %s sum %s vcov %s ci %s aic %s update %s cook %s hat %s",
     nm, paste(unlist(res), collapse = ","), if (is_err(ft)) msg(ft, 50) else "ok",
     dpn, simd,
     if (is_err(rg)) msg(rg, 40) else paste(rg, collapse = "/"),
     msg(ll, 40), msg(sm, 40), msg(vc, 40), msg(ci, 40), msg(an, 40), ll_up,
     ckd, hvd)
}

# the core log_lik(), posterior_predict(), and frm_curve on one per chart
for (nm in c("Gamma", "beta", "von_mises")) {
  s <- sets[[nm]]
  d <- data.frame(x = x, y = s$y); d$y[mrows] <- NA
  fit <- frm(bf(y | mi() ~ x), family = s$f, data = d)
  cc <- frm(bf(y ~ x), family = s$f, data = d[-mrows, ])
  a <- tr(frm_curve(fit, "x", at = data.frame(x = c(-1, 0, 1))))
  pr("%-10s frm_curve: %s", nm, msg(a, 90))
  lp <- tr(frm_linpred(fit))
  pr("%-10s frm_linpred: %s", nm, msg(lp, 90))
  ce <- tr(conditional_effects(fit))
  pr("%-10s conditional_effects: %s", nm, msg(ce, 90))
  # an imputation is what the reader takes: mode on the chart vs density mean,
  # median and mode at the missing rows
  v <- frmtmb:::mi_values(fit, "y")[mrows]
  fv <- fitted(fit)[mrows, 1]
  pr("%-10s imputed value / fitted mean at a missing row: %s", nm,
     paste(sprintf("%.4f", (v / fv)[1:4]), collapse = " "))
}

# mi() read by another response, new data: the reader needs a value
set.seed(9501)
d <- data.frame(x = x, w = rgamma(n, 2, 2 / exp(eta)))
d$z <- 1 + 0.5 * d$w + rnorm(n, 0, 0.5)
d$w[mrows] <- NA
fit <- frm(bf(z ~ mi(w) + x) + gaussian() +
             bf(w | mi() ~ x) + stats::Gamma(link = "log"), data = d)
bz <- fixef(fit)
print(round(bz, 4))
nd <- data.frame(x = c(0, 1), w = c(1, 2))
a <- tr(predict(fit, newdata = nd, resp = "z"))
pr("predict(newdata) with mi(w) observed: %s", if (is_err(a)) msg(a) else paste(round(a[, 1], 4), collapse = " "))
a <- tr(predict(fit, newdata = data.frame(x = c(0, 1), w = c(1, NA)), resp = "z"))
pr("predict(newdata) with mi(w) NA: %s", if (is_err(a)) msg(a) else paste(round(a[, 1], 4), collapse = " "))
