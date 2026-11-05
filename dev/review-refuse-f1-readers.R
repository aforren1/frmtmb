# Short check, item 1: the readers on a fit with dropped rows, and the row
# accounting against brms. Item 3's placeholder cases. Seed 11100.
source("dev/review-refuse-lib.R")
library(frmtmb.sample)
set.seed(11100)
n <- 120
x <- rnorm(n); g <- factor(rep(1:12, 10))
d <- data.frame(x = x, z = rnorm(n), g = g, y = rgamma(n, 2, 2 / exp(0.2 + 0.4 * x)), w = 1 + (x > 0))
M <- seq(4, n, by = 12); d$y[M] <- NA
Gam <- stats::Gamma(link = "log")
a <- frm(bf(y | mi() ~ x), family = Gam, data = d)
b <- frm(bf(y ~ x), family = Gam, data = d[-M, ])
pr("nobs(fit) %d; attr(logLik, nobs) %s; complete-case %d, %s; df %s, %s", nobs(a), attr(logLik(a), "nobs"), nobs(b), attr(logLik(b), "nobs"), attr(logLik(a), "df"), attr(logLik(b), "df"))
pr("logLik %.8f %.8f; AIC %.6f %.6f; BIC %.6f %.6f; estimates diff %.1e", logLik(a), logLik(b), AIC(a), AIC(b), BIC(a), BIC(b), max(abs(unlist(a$estimates[c("beta", "betad")]) - unlist(b$estimates[c("beta", "betad")]))))
fx <- tr(E$brms_fixed_fit(brms::bf(y | mi() ~ x), brms::brmsfamily("gamma", "log"), d, a, ndraws = 3))
pr("brms Fixed_param: %s", msg(fx, 200))
if (!is_err(fx)) {
  pr("brms nobs %d; log_lik dims %s", stats::nobs(fx), tryCatch(paste(dim(brms::log_lik(fx)), collapse = "x"), error = function(e) substr(conditionMessage(e), 1, 90)))
  pr("brms loo: %s", tryCatch({l <- suppressWarnings(brms::loo(fx)); paste(dim(l$pointwise), collapse = "x")}, error = function(e) substr(conditionMessage(e), 1, 90)))
  pr("brms posterior_predict dims %s; fitted dims %s", paste(dim(brms::posterior_predict(fx)), collapse = "x"), paste(dim(brms::posterior_epred(fx)), collapse = "x"))
}
for (me in c("residuals", "fitted", "predict")) {
  v <- tr(do.call(me, list(a)))
  m <- if (is_err(v)) msg(v, 80) else {
    u <- if (is.matrix(v)) v[, 1] else v
    sprintf("length %d, NA at missing %d, NA elsewhere %d", length(u), sum(is.na(u[M])), sum(is.na(u[-M])))
  }
  pr("%-10s %s", me, m)
}
for (ty in c("pearson", "deviance", "standardized", "quantile")) {
  v <- tr(residuals(a, type = ty)); pr("residuals(%s): %s", ty, if (is_err(v)) msg(v, 80) else sprintf("NA at missing %d", sum(is.na((if (is.matrix(v)) v[, 1] else v)[M]))))
}
sm <- tr(simulate(a, nsim = 2, seed = 1)); pr("simulate: %s, length %s, NA at missing %s", msg(sm, 80), if (is_err(sm)) "-" else nrow(as.matrix(sm)), if (is_err(sm)) "-" else sum(is.na(as.matrix(sm)[M, 1])))
ck <- cooks.distance(a); ck0 <- cooks.distance(b)
pr("cooks: at missing %s; observed rows against complete case, max diff %.1e", paste(signif(range(ck[M]), 3), collapse = " "), max(abs(ck[-M] - ck0)))
inf <- tr(influence(a)); pr("influence: %s", msg(inf, 100))
ww <- tr(weights(a)); pr("weights(): %s", if (is_err(ww)) msg(ww, 80) else sprintf("length %d, at missing %s", length(ww), paste(unique(ww[M]), collapse = ",")))
a2 <- frm(bf(y | mi() ~ 1), family = Gam, data = d)
an <- tr(anova(a2, a)); pr("anova: %s", msg(an, 150)); if (!is_err(an)) print(an)
an0 <- anova(frm(bf(y ~ 1), family = Gam, data = d[-M, ]), b); print(an0)
ax <- tr(anova(a, b)); pr("anova against the complete-case fit: %s", msg(ax, 150))
dp <- tr(drop1(a)); pr("drop1: %s", msg(dp, 150)); if (!is_err(dp)) print(dp)
stp <- tr(step(a, trace = 0)); pr("step: %s", msg(stp, 150))
up <- tr(update(a, . ~ . + z)); pr("update + z: %s; nobs %s/%s; latent %d", msg(up, 150), if (is_err(up)) "-" else nobs(up), if (is_err(up)) "-" else attr(logLik(up), "nobs"), if (is_err(up)) -1L else length(up$estimates$miss))
up2 <- tr(update(a, newdata = transform(d, y = ifelse(seq_len(n) %in% M[1:3], 1, y)))); pr("update newdata, 3 filled: %s; nobs attr %s", msg(up2, 150), if (is_err(up2)) "-" else attr(logLik(up2), "nobs"))
s1 <- tr(vcov(a, type = "sandwich")); pr("vcov sandwich: %s", msg(s1, 170))
s2 <- tr(vcov(a, cluster = ~ g)); pr("vcov cluster: %s", msg(s2, 170))
s0 <- tr(vcov(b, cluster = ~ g)); pr("  on the complete-case fit: %s", msg(s0, 100))
bo <- tr(frm_bootstrap(a, nsim = 6, seed = 1)); pr("frm_bootstrap: %s", msg(bo, 170))
imp <- tr(frm(bf(y | mi() ~ x + (1 | g)), family = Gam, data = d, importance = 50)); pr("importance, all unread: %s", msg(imp, 170))
qd <- tr(frm(bf(y | mi() ~ x + (1 | g)), family = Gam, data = d, quadrature = TRUE)); pr("quadrature, all unread: %s", msg(qd, 170))
qr0 <- tr(frm(bf(y ~ x + (1 | g)), family = Gam, data = d[-M, ], quadrature = TRUE))
if (!is_err(qd) && !is_err(qr0)) pr("  logLik diff to the complete case %.2e", as.numeric(logLik(qd)) - as.numeric(logLik(qr0)))
rm <- tr(frm(bf(y | mi() ~ x), family = Gam, data = d, REML = TRUE)); rm0 <- frm(bf(y ~ x), family = Gam, data = d[-M, ], REML = TRUE)
pr("REML: %s; logLik diff %.2e", msg(rm, 100), if (is_err(rm)) NA else as.numeric(logLik(rm)) - as.numeric(logLik(rm0)))
rw <- tr(frm(bf(y | mi() + weights(w) ~ x), family = Gam, data = d)); rw0 <- frm(bf(y | weights(w) ~ x), family = Gam, data = d[-M, ])
pr("weights: logLik diff %.2e; nobs attr %s vs %s", as.numeric(logLik(rw)) - as.numeric(logLik(rw0)), attr(logLik(rw), "nobs"), attr(logLik(rw0), "nobs"))
cf <- tr(frm_curve(a, "x")); pr("frm_curve: %s", msg(cf, 120))
ce <- tr(conditional_effects(a)); pr("conditional_effects: %s", msg(ce, 120))
rr <- tr(bayes_R2(a)); pr("bayes_R2 on a fit: %s", msg(rr, 120))
sa <- tr(frm_sample(a, chains = 1, iter = 200, seed = 1, refresh = 0)); pr("frm_sample(fit): %s", msg(sa, 150))
if (!is_err(sa)) {
  sv <- names(posterior::as_draws_df(sa)); pr("  draws: %s", paste(head(sv, 8), collapse = " "))
  lo <- tr(loo(sa)); pr("  loo(draws): %s", msg(lo, 200))
}
# the same on brms: loo with a missing row
# item 3: the placeholder's corners
pr("--- placeholder")
d1 <- data.frame(x = x, y = 1.5); d1$y[M] <- NA
a <- tr(frm(bf(y | mi() ~ x), family = Gam, data = d1)); pr("all-equal column, Gamma: %s", msg(a, 170))
b <- tr(frm(bf(y ~ x), family = Gam, data = d1[-M, ])); pr("  complete case: %s", msg(b, 170))
a <- tr(frm(bf(y | mi() ~ x), family = gaussian(), data = d1)); pr("all-equal column, gaussian: %s", msg(a, 170))
b <- tr(frm(bf(y ~ x), family = gaussian(), data = d1[-M, ])); pr("  complete case: %s", msg(b, 170))
a <- tr(frm(bf(y | mi() ~ x), family = Beta(), data = transform(d1, y = y / 3))); pr("all-equal column, beta: %s", msg(a, 170))
d2 <- data.frame(x = x, g = g, y = rgamma(n, 2, 2 / exp(0.2 + 0.4 * x)))
d2$y[d2$g == "3"] <- NA
a <- tr(frm(bf(y | mi() ~ x + (1 | g)), family = Gam, data = d2)); b <- tr(frm(bf(y ~ x + (1 | g)), family = Gam, data = d2[d2$g != "3", ]))
pr("a group with no observed value: %s; logLik diff %.2e; levels in ranef %s", msg(a, 150), if (is_err(a) || is_err(b)) NA else as.numeric(logLik(a)) - as.numeric(logLik(b)), if (is_err(a)) "-" else nrow(ranef(a)$g))
pr("  brms: Nmi %s, group levels %s", brms::make_standata(brms::bf(y | mi() ~ x + (1 | g)), family = brms::brmsfamily("gamma", "log"), data = d2)$Nmi, brms::make_standata(brms::bf(y | mi() ~ x + (1 | g)), family = brms::brmsfamily("gamma", "log"), data = d2)$N_1)
if (!is_err(a)) {
  ft <- fitted(a, newdata = d2[d2$g == "3", ][1:2, ]); pr("  fitted in that group: %s", paste(round(ft[, 1], 3), collapse = " "))
  r3 <- ranef(a)$g["3", 1]; pr("  its ranef %.4e; se? %s", r3, tryCatch(format(ranef(a, condVar = TRUE)$g["3", 2], digits = 3), error = function(e) "-"))
}
d3 <- transform(d2, f = factor(g == "3"))
a <- tr(frm(bf(y | mi() ~ x + f), family = Gam, data = d3)); pr("a factor level all missing, fixed: %s %s", msg(a, 170), if (is_err(a)) "" else paste(round(a$estimates$beta, 3), collapse = ","))
b <- tr(brms::make_standata(brms::bf(y | mi() ~ x + f), family = brms::brmsfamily("gamma", "log"), data = d3)); pr("  brms %s", msg(b, 90))
if (!is_err(a)) pr("  SE of that level: %s", paste(format(sqrt(diag(vcov(a))), digits = 3), collapse = " "))
# a fixed effect that only the missing rows inform: flat now, on the new route
d4 <- data.frame(x = x, u = as.numeric(seq_len(n) %in% M), y = rgamma(n, 2, 2)); d4$y[M] <- NA
a <- tr(frm(bf(y | mi() ~ x + u), family = Gam, data = d4)); pr("a column only the missing rows carry: %s %s; warn %s", msg(a, 170), if (is_err(a)) "" else paste(round(a$estimates$beta, 3), collapse = ","), paste(substr(attr(a, "w"), 1, 100), collapse = "; "))
b <- tr(frm(bf(y ~ x + u), family = Gam, data = d4[-M, ])); pr("  complete case says: %s", msg(b, 170))
# a mixture's chart from its components, at the sampling route
kk <- rbinom(n, 1, 0.4)
d5 <- data.frame(x = x, y = rlnorm(n, ifelse(kk == 1, 0.8, -0.2), 0.6)); d5$y[M] <- NA
a <- tr(frm(bf(y | mi() ~ x), family = mixture(lognormal(), lognormal()), data = d5)); b <- tr(frm(bf(y ~ x), family = mixture(lognormal(), lognormal()), data = d5[-M, ]))
pr("mixture(lognormal x 2), unread, frm: %s; logLik diff %.2e", msg(a, 150), if (is_err(a) || is_err(b)) NA else as.numeric(logLik(a)) - as.numeric(logLik(b)))
o <- tr(frm(bf(y | mi() ~ x), family = mixture(lognormal(), lognormal()), data = d5, dry_run = "objective"))
pr("  objective: %s; chart %s", msg(o, 100), if (is_err(o)) "-" else o$frame$mi_map$y$chart)
o2 <- tr(frm(bf(y | mi() ~ x), family = mixture(lognormal(), Gam), data = d5, dry_run = "objective"))
pr("  mixture(lognormal, Gamma) chart %s", if (is_err(o2)) msg(o2, 100) else o2$frame$mi_map$y$chart)
o3 <- tr(frm(bf(y | mi() ~ x), family = mixture(gaussian(), lognormal()), data = d5, dry_run = "objective"))
pr("  mixture(gaussian, lognormal) chart %s", if (is_err(o3)) msg(o3, 100) else o3$frame$mi_map$y$chart)
o4 <- tr(frm(bf(y | mi() ~ x), family = mixture(Beta(), lognormal()), data = transform(d5, y = plogis(log(y))), dry_run = "objective"))
pr("  mixture(beta, lognormal) chart %s", if (is_err(o4)) msg(o4, 100) else o4$frame$mi_map$y$chart)
# read mixture stays refused for its reason
d6 <- transform(d5, z = rnorm(n)); a <- tr(frm(bf(y | mi() ~ x) + mixture(lognormal(), lognormal()) + bf(z ~ mi(y)) + gaussian(), data = d6))
pr("mixture, read: %s", msg(a, 150))
# the bootstrap: does a refit keep the unread rows out?
set.seed(11101)
fb <- getAnywhere("frm_bootstrap")$objs[[1]]
pr("frm_bootstrap formals: %s", paste(names(formals(fb)), collapse = " "))
a <- frm(bf(y | mi() ~ x), family = Gam, data = d)
b <- frm(bf(y ~ x), family = Gam, data = d[-M, ])
ba <- tr(frm_bootstrap(a, nsim = 60, seed = 5)); bb <- tr(frm_bootstrap(b, nsim = 60, seed = 5))
pr("bootstrap: %s | %s", msg(ba, 100), msg(bb, 100))
if (!is_err(ba)) { print(class(ba)); print(head(ba$table %||% summary(ba), 3)); print(head(bb$table %||% summary(bb), 3)) }
