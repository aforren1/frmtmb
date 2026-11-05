# Item 2: the chart at its edges. Each case states a design that is hard for a
# mode: wrap on von_mises, a U-shaped beta, a GEV support bound, a tiny
# exgaussian sigma, a kink, RE beside the latent, weights, the argument's
# validation, a family outside the 14. Seeds in each block.
source("dev/review-refuse-lib.R")
cmp <- function(lab, fj, dat, ...) {
  a <- tr(frm(fj$form, family = fj$fam, data = dat, ...))
  b <- tr(frm(fj$cc, family = fj$fam, data = dat[!is.na(dat$y), ], ...))
  if (is_err(a) || is_err(b)) {
    pr("%-34s joint: %s | cc: %s", lab, msg(a, 110), msg(b, 60)); return(invisible(list(a, b)))
  }
  v <- tr(frmtmb:::mi_values(a, "y"))
  pr("%-34s conv %d/%d gap/miss %+.4f joint [%s] cc [%s] maxgrad %.1e warn %s", lab,
     a$opt$convergence, b$opt$convergence,
     (as.numeric(logLik(a)) - as.numeric(logLik(b))) / sum(is.na(dat$y)),
     paste(sprintf("%.3f", c(a$estimates$beta, a$estimates$betad)), collapse = ","),
     paste(sprintf("%.3f", c(b$estimates$beta, b$estimates$betad)), collapse = ","),
     max(abs(a$obj$gr(a$opt$par))),
     substr(paste(attr(a, "w"), collapse = "; "), 1, 90))
  invisible(list(a, b))
}
F0 <- function(fam) list(form = bf(y | mi() ~ x), cc = bf(y ~ x), fam = fam)

# 1. von_mises with the mean near the window's end: the density is high at
# both ends of (-pi, pi), and a chart with a logit has two modes
set.seed(9601); n <- 200; x <- rnorm(n)
for (m0 in c(0, 2, 2.8, 3.05, 3.14)) for (k in c(2, 20)) {
  th <- m0 + 0.1 * x
  y <- atan2(sin(th + rnorm(n, 0, 1 / sqrt(k))), cos(th + rnorm(n, 0, 1 / sqrt(k))))
  d <- data.frame(x = x, y = y); d$y[seq(5, 195, 10)] <- NA
  r <- cmp(sprintf("von_mises mu=%.2f kappa~%g", m0, k), F0(von_mises()), d)
  if (!is.null(r) && !is_err(r[[1]])) {
    fit <- r[[1]]
    v <- frmtmb:::mi_values(fit, "y")[seq(5, 195, 10)]
    # the exact marginal of a missing value is 1, whatever the mode found
    pv <- fit$obj$env$last.par.best; idx <- which(names(pv) == "miss")
    b <- fit$obj$env$f(pv)
    pr("   imputed range [%.3f, %.3f]; distance to the window's end min %.4f; max |u| %.1f",
       min(v), max(v), pi - max(abs(v)), max(abs(pv[idx])))
  }
}

# 2. beta: a U shape and means near the bounds
set.seed(9602)
for (cfg in list(c(0, 0.3), c(0, 0.05), c(4, 2), c(4, 0.3), c(-5, 5))) {
  et <- cfg[1] + 0.5 * x; phi <- cfg[2]; mu <- plogis(et)
  y <- rbeta(n, mu * phi, (1 - mu) * phi)
  y <- pmin(pmax(y, 1e-10), 1 - 1e-10)
  d <- data.frame(x = x, y = y); d$y[seq(5, 195, 10)] <- NA
  r <- cmp(sprintf("beta eta0=%g phi=%g", cfg[1], cfg[2]), F0(Beta()), d)
  if (!is.null(r) && !is_err(r[[1]])) {
    u <- r[[1]]$estimates$miss
    pr("   max |u| %.2f", max(abs(u)))
  }
}

# 3. gen_extreme_value on the identity chart: a support that moves with the
# parameters. Negative xi: an upper bound. Positive xi: a lower bound
set.seed(9603)
for (xi in c(-0.8, -0.4, -0.1, 0.4, 0.9)) {
  e0 <- 0.3 + 0.5 * x
  y <- e0 + ((rexp(n))^(-xi) - 1) / xi
  d <- data.frame(x = x, y = y); d$y[seq(5, 195, 10)] <- NA
  cmp(sprintf("gev xi=%g", xi), F0(gen_extreme_value()), d)
}

# 4. exgaussian near its kink, skew_normal with a large alpha, student nu low
set.seed(9604)
for (sg in c(0.5, 0.05, 0.005)) {
  y <- 0.3 + 0.5 * x + rnorm(n, 0, sg) + rexp(n, 1)
  d <- data.frame(x = x, y = y); d$y[seq(5, 195, 10)] <- NA
  cmp(sprintf("exgaussian sigma=%g", sg), F0(exgaussian()), d)
}
for (al in c(20, 100)) {
  del <- al / sqrt(1 + al^2)
  y <- 0.3 + 0.5 * x + del * abs(rnorm(n)) + sqrt(1 - del^2) * rnorm(n)
  d <- data.frame(x = x, y = y); d$y[seq(5, 195, 10)] <- NA
  cmp(sprintf("skew_normal alpha=%g", al), F0(skew_normal()), d)
}

# 5. Gamma: random effects beside the latent; weights; REML; data with a
# missing row at an extreme x; nearly all missing
set.seed(9605)
g <- factor(rep(1:20, each = 10))
e0 <- 0.3 + 0.5 * x + rnorm(20, 0, 0.5)[g]
y <- rgamma(n, 2, 2 / exp(e0))
d <- data.frame(x = x, y = y, g = g, w = sample(1:3, n, TRUE)); d$y[seq(5, 195, 10)] <- NA
r <- tr(frm(bf(y | mi() ~ x + (1 | g)), family = stats::Gamma(link = "log"), data = d))
r2 <- tr(frm(bf(y ~ x + (1 | g)), family = stats::Gamma(link = "log"), data = d[!is.na(d$y), ]))
pr("Gamma RE: %s; logLik gap/miss %.4f (cc %.3f); sd %s vs %s", msg(r),
   (as.numeric(logLik(r)) - as.numeric(logLik(r2))) / 20, as.numeric(logLik(r2)),
   format(VarCorr(r)$g$sd[1], digits = 4), format(VarCorr(r2)$g$sd[1], digits = 4))
r <- tr(frm(bf(y | mi() + weights(w) ~ x), family = stats::Gamma(link = "log"), data = d))
pr("Gamma mi + weights: %s", msg(r))
if (!is_err(r)) {
  # with weight w the integrand in u is f^w |dy/du|, which does not integrate
  # to 1: a missing row's weight has a meaning in brms too (its Ymi has f^w)
  u <- r$estimates$miss; pr("   u range %.3f %.3f", min(u), max(u))
}
r <- tr(frm(bf(y | mi() ~ x), family = stats::Gamma(link = "log"), data = d, REML = TRUE))
pr("Gamma mi REML: %s", msg(r))
d2 <- d; d2$y[1:150] <- NA
r <- tr(frm(bf(y | mi() ~ x), family = stats::Gamma(link = "log"), data = d2))
pr("Gamma 75%% missing: %s", msg(r))
d3 <- d[1:3, ]; d3$y <- c(NA, NA, NA)
r <- tr(frm(bf(y | mi() ~ x), family = stats::Gamma(link = "log"), data = rbind(d3, d[10:11, ])))
pr("Gamma, 3 of 5 missing: %s", msg(r))
d4 <- d; d4$y <- NA_real_
r <- tr(frm(bf(y | mi() ~ x), family = stats::Gamma(link = "log"), data = d4))
pr("Gamma, all missing: %s", msg(r))

# 6. mi(sd) with a chart: a measured value outside the support, and with se
set.seed(9606)
d <- data.frame(x = x, y = rgamma(n, 3, 3 / exp(0.5 + 0.2 * x)), s = 0.3)
d2 <- d; d2$y[3] <- -0.1
r <- tr(frm(bf(y | mi(s) ~ x), family = stats::Gamma(link = "log"), data = d2))
pr("Gamma mi(s), a measured y below 0: %s", msg(r, 150))
rb <- tr(brms::make_standata(brms::bf(y | mi(s) ~ x), family = stats::Gamma(link = "log"), data = d2))
pr("  brms on the same: %s", msg(rb, 150))
d3 <- d; d3$y[3] <- 0
r <- tr(frm(bf(y | mi(s) ~ x), family = stats::Gamma(link = "log"), data = d3))
pr("Gamma mi(s), a measured y = 0: %s", msg(r, 150))
rb <- tr(brms::make_standata(brms::bf(y | mi(s) ~ x), family = stats::Gamma(link = "log"), data = d3))
pr("  brms on the same: %s", msg(rb, 150))
d4 <- data.frame(x = x, y = runif(n, -3, 3), s = 0.3); d4$y[3] <- pi
r <- tr(frm(bf(y | mi(s) ~ x), family = von_mises(), data = d4))
pr("von_mises mi(s), a measured y = pi: %s", msg(r, 150))
rb <- tr(brms::make_standata(brms::bf(y | mi(s) ~ x), family = brms::von_mises(), data = d4))
pr("  brms on the same: %s", msg(rb, 150))
d5 <- d4; d5$y[3] <- NA
r <- tr(frm(bf(y | mi() ~ x), family = von_mises(), data = d5))
pr("von_mises mi(), observed y = 3.5 elsewhere? %s", msg(tr(frm(bf(y | mi() ~ x), family = von_mises(),
   data = transform(d5, y = ifelse(seq_len(n) == 9, 3.5, y)))), 150))

# 7. mi_chart's validation and defaults
fam_args <- list(family = "tfam", dpars = "mu", links = list(mu = "log"),
  lpdf = function(y, dpars, aterms) {
    RTMB::dgamma(y, 2, dpars[["mu"]] / 2, log = TRUE)
  })
for (v in list("foo", c("log", "logit"), NA_character_, 1, "lo", "c", "log", "no_mode", "identity", TRUE)) {
  r <- tr(do.call(frmtmb_family, c(fam_args, list(accepts_aterms = "mi", mi_chart = v))))
  pr("frmtmb_family(mi_chart = %s): %s%s", deparse(v), msg(r, 120),
     if (!is_err(r)) paste0(" -> ", r[["mi_chart"]]) else "")
}
set.seed(9607)
d <- data.frame(x = x, y = rgamma(n, 2, 2 / exp(0.3 + 0.5 * x))); d$y[seq(5, 195, 10)] <- NA
# a positive family of its own: it declares mi and no chart, the default
fp <- do.call(frmtmb_family, c(fam_args, list(accepts_aterms = "mi",
  valid_y = function(y) all(y > 0))))
r <- tr(frm(bf(y | mi() ~ x), family = fp, data = d))
pr("custom positive family, mi and no chart: %s; est %s", msg(r, 150),
   if (is_err(r)) "" else paste(round(r$estimates$beta, 3), collapse = ","))
fp <- do.call(frmtmb_family, c(fam_args, list(accepts_aterms = "mi")))
r <- tr(frm(bf(y | mi() ~ x), family = fp, data = d))
pr("custom positive family no valid_y, mi, no chart: %s; est %s", msg(r, 150),
   if (is_err(r)) "" else paste(round(r$estimates$beta, 3), collapse = ","))
fp <- do.call(frmtmb_family, c(fam_args, list(accepts_aterms = "weights", mi_chart = "log")))
r <- tr(frm(bf(y | mi() ~ x), family = fp, data = d))
pr("custom family with a chart and no mi: %s", msg(r, 150))
fp <- do.call(frmtmb_family, c(fam_args, list(mi_chart = "log")))
r <- tr(frm(bf(y | mi() ~ x), family = fp, data = d))
pr("custom family, no accepts_aterms: %s", msg(r, 150))
cf <- tr(custom_family("tcf", dpars = c("mu"), links = "log",
   lpdf = function(y, dpars, aterms) RTMB::dgamma(y, 2, dpars[["mu"]] / 2, log = TRUE)))
pr("custom_family: %s", msg(cf, 120))
if (!is_err(cf)) {
  r <- tr(frm(bf(y | mi() ~ x), family = cf, data = d))
  pr("custom_family with mi: %s", msg(r, 150))
}

# 8. a family outside the 14, against brms
set.seed(9608)
d <- data.frame(x = x, y = rpois(n, 2)); d$y[3] <- NA
for (ft in c("poisson()", "negbinomial()", "shifted_lognormal()", "hurdle_gamma()",
             "zero_inflated_beta()", "student()")) {
  dd <- d
  if (ft == "zero_inflated_beta()") dd$y <- ifelse(is.na(dd$y), NA, (dd$y + 0.5) / 12)
  if (ft == "student()") dd$y <- ifelse(is.na(dd$y), NA, dd$y + 0.1)
  if (ft %in% c("shifted_lognormal()", "hurdle_gamma()")) dd$y <- dd$y + 0.5
  a <- tr(frm(bf(y | mi() ~ x), family = eval(parse(text = ft), envir = asNamespace("frmtmb")), data = dd, dry_run = "frame"))
  b <- tr(brms::make_stancode(brms::bf(y | mi() ~ x), family = eval(parse(text = ft), envir = asNamespace("brms")), data = dd))
  pr("%-22s frm: %s | brms: %s", ft, msg(a, 70), msg(b, 70))
}

# 9. mi() x subset(), mi() x trials(): a charted row outside the subset
set.seed(9609)
d <- data.frame(x = x, y = rgamma(n, 2, 2 / exp(0.3 + 0.5 * x)), z = rnorm(n), s = TRUE)
d$y[seq(5, 195, 10)] <- NA; d$s[seq(5, 95, 10)] <- FALSE
a <- tr(frm(bf(y | mi() + subset(s) ~ x) + stats::Gamma(link = "log") + bf(z ~ x) + gaussian(), data = d))
b <- tr(brms::make_standata(brms::bf(y | mi() + subset(s) ~ x, family = stats::Gamma(link = "log")) + brms::bf(z ~ x), data = d))
pr("mi + subset: frm %s (nmiss %s) | brms %s (Nmi %s)", msg(a, 110),
   if (is_err(a)) "-" else length(a$estimates$miss), msg(b, 70), if (is_err(b)) "-" else b$Nmi_y)
