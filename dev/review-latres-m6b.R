# Short check, item 4: response_mean_gh(). Exact mean on families with their own
# mean, nodes at a large variance, its SEs through fit_fd_se(), and marginaleffects
# against fitted() and brms.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
suppressPackageStartupMessages(library(rstan))
hel <- new.env(parent = asNamespace("frmtmb"))
hel$`%||%` <- `%||%`
sys.source("dev/equiv-lib.R", envir = hel)
gh <- fi("response_mean_gh")
rm_ <- fi("response_mean")
dpn <- fi("dpars_natural")
set.seed(167)
d <- expand.grid(t = 1:6, g = factor(1:25))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
e <- as.vector(replicate(25, arima.sim(list(ar = .5), 6, sd = .5)))
mu <- exp(.8 + .3 * d$x + e)
d$cnt <- rpois(n, mu)
d$zi <- d$cnt * rbinom(n, 1, .7)
d$hp <- ifelse(rbinom(n, 1, .3) == 1, 0, pmax(rpois(n, mu), 1))
d$pos <- rgamma(n, 4, 4 / mu)
d$hg <- d$pos * rbinom(n, 1, .7)
d$hl <- d$pos * rbinom(n, 1, .7)
d$un <- pmin(pmax(rbeta(n, 8 * plogis(e), 8 * (1 - plogis(e))), 1e-3), 1 - 1e-3)
d$zb <- d$un * rbinom(n, 1, .8)
d$bi <- rbinom(n, 8, plogis(e))
d$tr <- 8L
d$ord <- factor(findInterval(.3 * d$x + e + rlogis(n), c(-1, 0, 1)) + 1L, ordered = TRUE)
d$ex <- rep(c(1, 2), length.out = n)
nd <- data.frame(t = 1:3, g = factor(rep("new", 3)), x = c(-1, 0, 1), tr = 8L, ex = 2)
# closed forms by integrate(): E[mean | eta + Z], per family, at sd
intg <- function(f, eta, s) sapply(eta, function(e0) integrate(function(u) f(e0 + s * u) * dnorm(u), -12, 12, rel.tol = 1e-10, subdivisions = 2000L)$value)
cases <- list(
  poisson = list(y ~ x + ar(t, g, cov = TRUE), poisson(), function(e, p) exp(e)),
  lognormal = list(pos ~ x + ar(t, g, cov = TRUE), lognormal(), function(e, p) exp(e + p$sigma^2 / 2)),
  zip = list(zi ~ x + ar(t, g, cov = TRUE), zero_inflated_poisson(), function(e, p) (1 - p$zi) * exp(e)),
  hurdle_poisson = list(hp ~ x + ar(t, g, cov = TRUE), hurdle_poisson(), function(e, p) (1 - p$hu) * exp(e) / (-expm1(-exp(e)))),
  hurdle_gamma = list(hg ~ x + ar(t, g, cov = TRUE), hurdle_gamma(), function(e, p) (1 - p$hu) * exp(e)),
  hurdle_lognormal = list(hl ~ x + ar(t, g, cov = TRUE), hurdle_lognormal(), function(e, p) (1 - p$hu) * exp(e + p$sigma^2 / 2)),
  zero_inflated_beta = list(zb ~ x + ar(t, g, cov = TRUE), zero_inflated_beta(), function(e, p) (1 - p$zi) * plogis(e)),
  binomial = list(bi | trials(tr) ~ x + ar(t, g, cov = TRUE), binomial(), function(e, p) 8 * plogis(e)),
  rate = list(cnt | rate(ex) ~ x + ar(t, g, cov = TRUE), poisson(), function(e, p) 2 * exp(e)),
  shifted_lognormal = list(pos ~ x + ar(t, g, cov = TRUE), shifted_lognormal(), function(e, p) exp(e + p$sigma^2 / 2) + p$ndt))
names(cases)[1] <- "poisson"
cases$poisson[[1]] <- cnt ~ x + ar(t, g, cov = TRUE)
for (nm in names(cases)) {
  cc <- cases[[nm]]
  m <- msg_of(frm(cc[[1]], family = cc[[2]], data = d))
  if (!is.null(m$err)) { cat(nm, "ERR", substr(m$err, 1, 150), "\n"); next }
  f <- m$val
  bk <- fi("frame_latres")(f$frame)[[1]]
  th <- f$estimates$theta
  ms <- fi("latres_marginal_sd")(th[bk$theta_idx], bk$aux)
  rspec <- f$spec$responses[[1]]
  dp <- dpn(f, rspec, nd, NULL, TRUE)
  eta <- as.numeric(frm_linpred(f, newdata = nd, type = "link", dpar = "mu", allow_new_levels = TRUE))
  fv <- fitted(f, newdata = nd, allow_new_levels = TRUE)
  p <- lapply(dp, function(v) v[1])
  ex <- intg(function(e) cc[[3]](e, p), eta, ms)
  cat(sprintf("%-20s marg sd %.3f  fitted %s  exact %s  max rel %.2e  SE %s\n", nm, ms,
              paste(sprintf("%.5f", fv[, 1]), collapse = " "), paste(sprintf("%.5f", ex), collapse = " "), max(abs(fv[, 1] / ex - 1)),
              paste(sprintf("%.3f", fv[, 2]), collapse = " ")))
}
cat("\n== nodes at large variance, against integrate and closed forms\n")
for (s in c(1, 2.3, 4, 6)) {
  fam <- poisson()
  ff <- frm(cnt ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d)$spec$responses[[1]]$family
  v <- gh(0, s^2, ff, list(mu = 0), mean_fn = function(fam, dpars, aterms) dpars$mu)
  cat(sprintf("log link, sd %.1f: gh %.6g closed %.6g rel %.2e\n", s, v, exp(s^2 / 2), v / exp(s^2 / 2) - 1))
  fb <- frm(bi | trials(tr) ~ x + ar(t, g, cov = TRUE), family = binomial(), data = d)$spec$responses[[1]]$family
  v2 <- gh(1, s^2, fb, list(mu = 0), mean_fn = function(fam, dpars, aterms) dpars$mu)
  x2 <- intg(plogis, 1, s)
  cat(sprintf("logit, sd %.1f: gh %.8f integrate %.8f rel %.2e\n", s, v2, x2, v2 / x2 - 1))
  # a hurdle poisson at low mean, sharp near 0
  v3 <- gh(-2, s^2, ff, list(mu = 0), mean_fn = function(fam, dpars, aterms) dpars$mu / (1 - exp(-dpars$mu)))
  x3 <- intg(function(e) exp(e) / (-expm1(-exp(e))), -2, s)
  cat(sprintf("truncated poisson mean, eta -2, sd %.1f: gh %.8g integrate %.8g rel %.2e\n", s, v3, x3, v3 / x3 - 1))
}
cat("\n== cumulative category means, brms's per-draw, truncation\n")
fo <- frm(ord ~ x + ar(t, g, cov = TRUE), family = cumulative(), data = d)
fvo <- msg_of(fitted(fo, newdata = nd, allow_new_levels = TRUE))
cat("cumulative fitted(new) dims:", dim(fvo$val), " err:", fvo$err, "\n"); print(round(fvo$val[, , 1, drop = FALSE], 4))
bk <- fi("frame_latres")(fo$frame)[[1]]
ms <- fi("latres_marginal_sd")(fo$estimates$theta[bk$theta_idx], bk$aux)
ftr <- msg_of(frm(cnt | trunc(ub = 100) ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d))
if (is.null(ftr$err)) {
  fv <- fitted(ftr$val, newdata = transform(nd), allow_new_levels = TRUE)
  cat("trunc fitted new:", round(fv[, 1], 4), "\n")
}

cat("\n== SE through fit_fd_se() against the delta method on the closed form (poisson)\n")
f <- frm(cnt ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d)
nd2 <- data.frame(t = 1:6, g = factor(rep("new", 6)), x = 0)
fv <- fitted(f, newdata = nd2, allow_new_levels = TRUE)
V <- as.matrix(vcov(f, full = TRUE))
th <- f$estimates$theta; b <- f$estimates$beta
# E = exp(b0 + b1 x + s^2 (1 + t^2) / 2)
msd2 <- exp(2 * th[1]) * (1 + th[2]^2)
gr <- c(1, 0, msd2, exp(2 * th[1]) * th[2])
E <- exp(b[1] + msd2 / 2)
cat("value", fv[1, 1], E, " SE fd", fv[1, 2], " SE closed", E * sqrt(drop(t(gr) %*% V %*% gr)), " interval", fv[1, 3], fv[1, 4], "\n")
cat("at a boundary fit, SE:\n")
set.seed(1)
db <- expand.grid(t = 1:6, g = factor(1:10))
db$x <- rnorm(nrow(db))
db$y <- rpois(nrow(db), exp(.5 + .3 * db$x + rep(rnorm(10, 0, .5), each = 6)))
fb <- msg_of(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = db))$val
r <- msg_of(fitted(fb, newdata = db[1:3, ]))
print(r$val); cat(substr(r$warn, 1, 200), "\n")

cat("\n== marginaleffects against fitted() and brms\n")
nd3 <- data.frame(t = 1:6, g = factor(1), x = 0)
pr <- marginaleffects::predictions(f, newdata = nd3, vcov = FALSE)
cat("predictions", round(pr$estimate[1:3], 5), " fitted", round(fitted(f, newdata = nd3)[1:3, 1], 5), " linpred at 0", round(exp(frm_linpred(f, newdata = nd3, type = "link"))[1:3], 5), "\n")
pall <- marginaleffects::avg_predictions(f, vcov = FALSE)
cat("avg_predictions", pall$estimate, " mean fitted(newdata = d)", mean(fitted(f, newdata = d)[, 1]), " mean fitted()", mean(fitted(f)[, 1]), " mean y", mean(d$cnt), "\n")
sl <- marginaleffects::avg_slopes(f, variables = "x", vcov = FALSE)
cat("avg_slopes", sl$estimate, " closed form b1 * mean(E)", b[2] * mean(fitted(f, newdata = d)[, 1]), "\n")
slse <- msg_of(marginaleffects::avg_slopes(f, variables = "x"))
cat("avg_slopes se", slse$val$std.error, " err", slse$err, " warn", substr(slse$warn, 1, 120), "\n")
pse <- msg_of(marginaleffects::predictions(f, newdata = nd3))
cat("predictions se", round(pse$val$std.error[1:3], 4), " fitted's SE", round(fitted(f, newdata = nd3)[1:3, 2], 4), "\n")
hyp <- msg_of(marginaleffects::predictions(f, newdata = nd3, type = "link", vcov = FALSE))
cat("type link", round(hyp$val$estimate[1:3], 4), "\n")
fb2 <- brms::bf(cnt ~ x + ar(t, g, cov = TRUE))
code <- brms::stancode(fb2, d, family = brms::brmsfamily("poisson"))
sdat <- brms::standata(fb2, d, family = brms::brmsfamily("poisson"))
mod <- hel$brms_stan_model(code)
sf <- suppressWarnings(rstan::sampling(mod, data = sdat, chains = 4, iter = 1500, warmup = 750, seed = 11, refresh = 0))
scaf <- brms::brm(fb2, d, family = brms::brmsfamily("poisson"), empty = TRUE)
scaf$fit <- sf
bfit <- brms::rename_pars(scaf)
ebn <- brms::posterior_epred(bfit, newdata = nd3)
cat("brms epred new mean", round(mean(ebn), 4), " frm fitted", round(mean(fitted(f, newdata = nd3)[, 1]), 4), " marginaleffects", round(mean(pr$estimate), 4), " b0 etc: brms sderr", round(median(rstan::extract(sf, "sderr")[[1]]), 3), " frm", round(exp(th[1]), 3), "\n")
mb <- msg_of(marginaleffects::predictions(bfit, newdata = nd3))
cat("marginaleffects on brms:", if (is.null(mb$err)) round(mb$val$estimate[1:3], 4) else mb$err, "\n")
mbs <- msg_of(marginaleffects::avg_slopes(bfit, variables = "x"))
cat("avg_slopes on brms:", if (is.null(mbs$err)) mbs$val$estimate else mbs$err, "  frm", sl$estimate, "\n")
ce <- conditional_effects(f, "x")$x
cat("conditional_effects at err = 0, frm:", round(ce$estimate__[50], 4), " exp(eta):", round(exp(b[1] + b[2] * ce$x[50]), 4), "\n")
em <- msg_of(summary(emmeans::emmeans(f, ~ x, at = list(x = 0), type = "response")))
cat("emmeans at 0:", if (is.null(em$err)) em$val[1, 2] else em$err, " exp(b0)", exp(b[1]), "\n")
