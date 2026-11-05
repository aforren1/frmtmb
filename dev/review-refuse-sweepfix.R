# Item 6 and 7: the sweep's fixes on designs of mine, and the refusals' stated
# reasons. Seed 10000.
source("dev/review-refuse-lib.R")
show <- function(lab, r) {
  if (is_err(r)) { pr("%-38s ERROR %s", lab, msg(r, 170)); return(invisible()) }
  pr("%-38s verdict=%s const=%s pts=[%s] grad=%s frm[%s] brms[%s] l1[%s]", lab, r$verdict,
     format(r$l1_const %||% NA, digits = 3),
     paste(format(unlist(r$l1_extra_const) %||% NA, digits = 3), collapse = ","),
     format(r$l1_grad %||% NA, digits = 3),
     substr(r$frm_error %||% "", 1, 110), substr(r$brms_error %||% "", 1, 110),
     substr(r$l1_error %||% "", 1, 110))
}
set.seed(10000)
n <- 180
x <- rnorm(n); z <- rnorm(n); s <- 0.1 + 0.4 * runif(n)
g <- factor(sample(1:10, n, TRUE))

# se() x skew_normal: varying alpha, a dpar formula, an RE, large skew
del <- 4 / sqrt(17)
d <- data.frame(x = x, z = z, s = s, g = g,
                y = 0.5 + 0.7 * x + 0.8 * (del * abs(rnorm(n)) + sqrt(1 - del^2) * rnorm(n)) + rnorm(n, 0, s))
show("se(s) skew_normal, RE", tr(meas(list(bf = "bf(y | se(s) ~ x + (1 | g))", family = "skew_normal()", data = d))))
show("se(s, sigma = TRUE) skew_normal", tr(meas(list(bf = "bf(y | se(s, sigma = TRUE) ~ x)", family = "skew_normal()", data = d))))
show("se(s) skew_normal, alpha ~ z", tr(meas(list(bf = "bf(y | se(s) ~ x, alpha ~ z)", family = "skew_normal()", data = d))))
show("se(s, sigma = TRUE), sigma ~ z", tr(meas(list(bf = "bf(y | se(s, sigma = TRUE) ~ x, sigma ~ z)", family = "skew_normal()", data = d))))
# simulate, predict: the known sd reaches both
ft <- frm(bf(y | se(s) ~ x), family = skew_normal(), data = d)
sm <- simulate(ft, nsim = 400, seed = 7)
resid_sd <- apply(as.matrix(sm), 1, sd)
# a skew normal has sd = omega sqrt(1 - 2 delta^2 / pi), by the fitted alpha
al <- ft$estimates$betad[[1]]; dl <- al / sqrt(1 + al^2)
want <- d$s * sqrt(1 - 2 * dl^2 / pi)
pr("skew_normal se(): simulate sd / family sd at se: median %.3f, range %.3f %.3f (alpha %.3f)", median(resid_sd / want), min(resid_sd / want), max(resid_sd / want), al)
pp <- tr(predict(ft)); pr("  predict: %s", msg(pp, 100))
ft2 <- fitted(ft); mh <- ft$estimates$beta[[1]] + ft$estimates$beta[[2]] * d$x + d$s * dl * sqrt(2 / pi)
pr("  fitted vs hand mean (mu + se delta sqrt(2/pi)): max diff %.3e; vs mu only %.3e", max(abs(ft2[, 1] - mh)), max(abs(ft2[, 1] - (ft$estimates$beta[[1]] + ft$estimates$beta[[2]] * d$x))))
# brms's fitted for the same parameter values on the same rows, through a
# Fixed_param brmsfit, for mu and epred
bfit <- tr(E$brms_fixed_fit(brms::bf(y | se(s) ~ x), brms::skew_normal(), d, ft, ndraws = 3))
if (!is_err(bfit)) {
  be <- brms::posterior_epred(bfit)[1, ]
  pr("  brms posterior_epred at these values vs frmtmb fitted: %.3e", max(abs(be - ft2[, 1])))
} else pr("  brms fixed fit: %s", msg(bfit, 150))
nd <- d[1:5, ]
a <- tr(predict(ft, newdata = nd)); pr("  predict(newdata): %s", msg(a, 90))
nd2 <- d[1:5, c("x", "g")]
a <- tr(predict(ft, newdata = nd2)); pr("  predict(newdata without s): %s", msg(a, 150))
a <- tr(frm(bf(y | se(s) ~ x), family = skew_normal(), data = transform(d, s = ifelse(seq_len(n) == 1, 0, s))))
pr("  se = 0 on a row: %s", msg(a, 150))
# compat row reads
a <- tr(residuals(ft, type = "pearson")); pr("  residuals(pearson): %s", msg(a, 100))
a <- tr(frmtmb::log_lik)
la <- tr(logLik(ft)); pr("  logLik %s", format(as.numeric(la), digits = 8))
hand <- sum(log(2 / d$s) + dnorm((d$y - ft$estimates$beta[[1]] - ft$estimates$beta[[2]] * d$x) / d$s, log = TRUE) +
              pnorm(al * (d$y - ft$estimates$beta[[1]] - ft$estimates$beta[[2]] * d$x) / d$s, log.p = TRUE))
pr("  logLik - hand %.3e", as.numeric(la) - hand)

# mi() x se(): a missing row, the known sd; read by a predictor; mi(sd) x se()
d2 <- d; d2$y[seq(4, n, 11)] <- NA
show("mi() + se(s), read by mi(y)", tr(meas(list(
  bf = "mvbf(bf(y | mi() + se(s) ~ x), bf(z ~ mi(y)), rescor = FALSE)", family = NA, data = d2))))
show("mi() + se(s, sigma = TRUE)", tr(meas(list(bf = "bf(y | mi() + se(s, sigma = TRUE) ~ x)", data = d2))))
show("mi() + se(s) student", tr(meas(list(bf = "bf(y | mi() + se(s) ~ x)", family = "student()", data = d2))))
show("mi() + se(s) skew_normal", tr(meas(list(bf = "bf(y | mi() + se(s) ~ x)", family = "skew_normal()", data = d2))))
d3 <- d2; d3$sd <- 0.2
show("mi(sd) + se(s) gaussian", tr(meas(list(bf = "bf(y | mi(sd) + se(s) ~ x)", data = d3))))
d4 <- d2
d4$s[seq(4, n, 11)] <- NA
a <- tr(frm(bf(y | mi() + se(s) ~ x), data = d4)); b <- tr(brms::make_standata(brms::bf(y | mi() + se(s) ~ x), data = d4))
pr("se NA at a missing y: frm %s | brms %s", msg(a, 150), msg(b, 150))
a <- tr(frm(bf(y | mi() + se(s) ~ x), data = d2))
if (!is_err(a)) {
  p <- tr(fitted(a)); r <- tr(residuals(a)); si <- tr(simulate(a, nsim = 1, seed = 1))
  pr("mi+se methods: fitted %s residuals %s simulate %s", msg(p, 60), msg(r, 60), msg(si, 60))
  # the exact: a gaussian mi()+se() is complete-case, Laplace exact
  cc <- frm(bf(y | se(s) ~ x), data = d2[!is.na(d2$y), ])
  pr("  mi+se vs complete case: logLik gap %.3e, beta gap %.3e", as.numeric(logLik(a)) - as.numeric(logLik(cc)), max(abs(a$estimates$beta - cc$estimates$beta)))
}

# gr(id =, cor = FALSE): by hand, a key with a different term that does
# correlate. brms: an id with cor = FALSE on one and cor = TRUE on the other
set.seed(10001)
d5 <- data.frame(x = x, w = z, g = g, y = rnorm(n) + 0.5 * x)
fr <- function(lab, form) {
  a <- tr(frm(eval(parse(text = form), envir = asNamespace("frmtmb")), data = d5))
  b <- tr(brms::make_stancode(eval(parse(text = form), envir = asNamespace("brms")), data = d5))
  blocks <- if (is_err(a)) "-" else paste(vapply(a$frame$re_blocks, function(k) paste0(k$covstruct, k$dim), ""), collapse = ",")
  ll <- if (is_err(a)) NA else as.numeric(logLik(a))
  zb <- if (is_err(b)) "-" else paste(regmatches(b, gregexpr("(z_[0-9]+;|cor_[0-9]+;|L_[0-9]+;)", b))[[1]], collapse = " ")
  pr("%-34s frm %s blocks %s ll %.4f | brms %s [%s]", lab, msg(a, 110), blocks, ll, msg(b, 110), substr(zb, 1, 90))
}
fr("id cor F, one term", "bf(y ~ x + (1 + w | gr(g, id = 'p', cor = FALSE)))")
fr("|p| cor F", "bf(y ~ x + (1 + w | p | gr(g, cor = FALSE)))")
fr("id cor F on both, 2 dpars", "bf(y ~ x + (1 + w | gr(g, id = 'p', cor = FALSE)), sigma ~ (1 | gr(g, id = 'p', cor = FALSE)))")
fr("mixed: cor F + cor T, one key", "bf(y ~ x + (1 | gr(g, id = 'p', cor = FALSE)), sigma ~ (1 | gr(g, id = 'p')))")
fr("mixed: cor T + cor F, one key", "bf(y ~ x + (1 | gr(g, id = 'p')), sigma ~ (1 | gr(g, id = 'p', cor = FALSE)))")
fr("mm key cor F", "bf(y ~ x + (1 + w | p | mm(g, g, cor = FALSE)))")
# what the mixed spelling means in the two: a correlation exists or not
mx <- tr(frm(bf(y ~ x + (1 | gr(g, id = "p", cor = FALSE)), sigma ~ (1 | gr(g, id = "p"))), data = d5))
if (!is_err(mx)) print(sapply(mx$frame$re_blocks, function(b) paste(b$covstruct, b$dim, b$term_label %||% "")))
bc <- tr(brms::make_stancode(brms::bf(y ~ x + (1 | gr(g, id = "p", cor = FALSE)), sigma ~ (1 | gr(g, id = "p"))), data = d5))
if (!is_err(bc)) writeLines(grep("cor_|z_1|L_1|sd_1", strsplit(bc, "\n")[[1]], value = TRUE)[1:6])
bc <- tr(brms::make_standata(brms::bf(y ~ x + (1 | gr(g, id = "p", cor = FALSE)), sigma ~ (1 | gr(g, id = "p"))), data = d5))

# cens() x mixture: the gradient, the tails by hand, mixed components, by the
# log_prob at several points; the lane's 2.6e-4 gradient against the control's
set.seed(10002)
k <- rbinom(n, 1, 0.35)
y <- ifelse(k == 1, 2.5 + 0.5 * x, -2 - 0.3 * x) + rnorm(n)
ce <- ifelse(y > 3, 1L, ifelse(y < -3, -1L, 0L))
yy <- pmax(pmin(y, 3), -3)
d6 <- data.frame(x = x, y = yy, ce = ce, y2 = yy + 1)
d6$y2 <- ifelse(ce == 2L, 1, 1)
show("control: no cens, mix none", tr(meas(list(bf = "bf(y ~ 1, mu1 ~ x, mu2 ~ x)", family = "mixture(gaussian(), gaussian(), order = 'none')", data = d6))))
show("cens mix right+left", tr(meas(list(bf = "bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x)", family = "mixture(gaussian(), gaussian(), order = 'none')", data = d6))))
d7 <- d6; d7$ce2 <- ce; d7$ce2[seq(3, n, 20)] <- 2L; d7$ub <- d7$y + 1
show("cens mix interval", tr(meas(list(bf = "bf(y | cens(ce2, ub) ~ 1, mu1 ~ x, mu2 ~ x)", family = "mixture(gaussian(), gaussian(), order = 'none')", data = d7))))
d8 <- d6; d8$y <- exp(d8$y / 2)
show("cens mix lognormal + Gamma", tr(meas(list(bf = "bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x)", family = "mixture(lognormal(), Gamma(link = 'log'), order = 'none')", data = d8))))
show("cens mix student + gaussian", tr(meas(list(bf = "bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x)", family = "mixture(student(), gaussian(), order = 'none')", data = d6))))
show("cens mix ordered (default)", tr(meas(list(bf = "bf(y | cens(ce) ~ x)", family = "mixture(gaussian(), gaussian())", data = d6))))
show("cens mix with theta ~ x", tr(meas(list(bf = "bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x, theta2 ~ x)", family = "mixture(gaussian(), gaussian(), order = 'none')", data = d6))))
show("cens mix weights", tr(meas(list(bf = "bf(y | cens(ce) + weights(w) ~ 1, mu1 ~ x, mu2 ~ x)", family = "mixture(gaussian(), gaussian(), order = 'none')", data = transform(d6, w = 1 + (x > 0))))))

# the gradient itself: nlminb's own stopping (control) against a tighter stop
cm <- frm(bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x), family = mixture(gaussian(), gaussian(), order = "none"), data = d6)
pr("cens mix: frm's own max|gr| %.2e at convergence %d, msg %s", max(abs(cm$obj$gr(cm$opt$par))), cm$opt$convergence, cm$opt$message)
# hand density of a censored mixture row at the fit
est <- cm$estimates
print(unlist(lapply(est[c("beta", "betad")], function(v) round(v, 4))))
# independent evaluation: logLik by R's pnorm
fe <- fixef(cm)
print(round(fe[, 1], 4))
mu1 <- fe["mu1_Intercept", 1] + fe["mu1_x", 1] * x
mu2 <- fe["mu2_Intercept", 1] + fe["mu2_x", 1] * x
s1 <- exp(est$betad[grep("sigma1", names(est$betad))]); s2 <- exp(est$betad[grep("sigma2", names(est$betad))])
th <- est$betad[grep("theta", names(est$betad))]
pr("theta chart names: %s", paste(names(est$betad), collapse = ","))
w2 <- plogis(th); w1 <- 1 - w2
ld <- ifelse(ce == 0, log(w1 * dnorm(yy, mu1, s1) + w2 * dnorm(yy, mu2, s2)),
      ifelse(ce == 1, log(w1 * pnorm(yy, mu1, s1, FALSE) + w2 * pnorm(yy, mu2, s2, FALSE)),
             log(w1 * pnorm(yy, mu1, s1) + w2 * pnorm(yy, mu2, s2))))
pr("cens mix: logLik - hand %.3e", as.numeric(logLik(cm)) - sum(ld))
for (me in c("predict", "fitted", "residuals", "simulate", "pp_mixture")) {
  v <- tr(do.call(me, list(cm))); pr("  %s: %s", me, msg(v, 100))
}
# group mixture refused; cens + trunc on a mixture
a <- tr(frm(bf(y | cens(ce) ~ x), family = mixture(gaussian(), gaussian(), groups = ~ gr), data = transform(d6, gr = g)))
pr("mixture(groups=) + cens: %s", msg(a, 100))
a <- tr(frm(bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x), family = mixture(poisson(), poisson(), order = "none"), data = transform(d6, y = rpois(n, 3))))
pr("cens on mixture(poisson, poisson): %s", msg(a, 120))
a <- tr(frm(bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x), family = mixture(gaussian(), Beta(), order = "none"), data = transform(d6, y = plogis(y))))
pr("cens on mixture(gaussian, beta): %s", msg(a, 120))
a <- tr(brms::make_stancode(brms::bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x), family = brms::mixture(gaussian(), brms::Beta(), order = "none"), data = transform(d6, y = plogis(y))))
pr("  brms: %s", msg(a, 120))

# thres() x mvbf
set.seed(10003)
o <- cut(0.7 * x + rnorm(n), 4, labels = FALSE)
d9 <- data.frame(x = x, z = z, o = o, w = rnorm(n, 1 + 0.3 * x), g = g, nt = 3L)
show("thres(3) + gaussian", tr(meas(list(bf = "mvbf(bf(o | thres(3) ~ x, family = cumulative()), bf(w ~ x), rescor = FALSE)", family = NA, data = d9))))
show("thres(gr=) mv", tr(meas(list(bf = "mvbf(bf(o | thres(3, gr = gg) ~ x, family = cumulative()), bf(w ~ x), rescor = FALSE)", family = NA, data = transform(d9, gg = factor(z > 0))))))
show("thres + (1|p|g) across responses", tr(meas(list(bf = "mvbf(bf(o | thres(3) ~ x + (1 | p | g), family = cumulative()), bf(w ~ x + (1 | p | g)), rescor = FALSE)", family = NA, data = d9))))
show("two ordinal", tr(meas(list(bf = "mvbf(bf(o | thres(3) ~ x, family = cumulative()), bf(o2 | thres(3) ~ x, family = sratio()), rescor = FALSE)", family = NA, data = transform(d9, o2 = cut(x + rnorm(n), 4, labels = FALSE))))))
fit <- tr(frm(mvbf(bf(o | thres(3) ~ x, family = cumulative()), bf(w ~ x), rescor = FALSE), data = d9))
pr("thres mv fit: %s", msg(fit, 150))
if (!is_err(fit)) {
  f1 <- frm(bf(o | thres(3) ~ x), family = cumulative(), data = d9); f2 <- frm(bf(w ~ x), data = d9)
  pr("  independent responses: logLik sum gap %.3e", as.numeric(logLik(fit)) - as.numeric(logLik(f1)) - as.numeric(logLik(f2)))
  for (me in c("predict", "fitted", "residuals")) {
    v <- tr(do.call(me, list(fit, resp = "o"))); pr("  %s(resp=o): %s", me, msg(v, 120))
  }
  v <- tr(simulate(fit, nsim = 1)); pr("  simulate: %s", msg(v, 120))
}

# refusals that remain: each reason, by running brms
# mi() x cens(): "brms reads Y, Inf at a missing row, log_prob -Inf at a finite point"
set.seed(10004)
d10 <- data.frame(x = x, y = rnorm(n, 0.5 * x), ce = 0L)
d10$y[seq(5, n, 9)] <- NA
bb <- brms::bf(y | mi() + cens(ce) ~ x)
code <- brms::make_stancode(bb, data = d10)
sd10 <- brms::make_standata(bb, data = d10)
writeLines(grep("Ymi|Yl|Y\\[|Y,|Y |Jmi|target", strsplit(code, "\n")[[1]], value = TRUE))
pr("brms's Y at a missing row: %s", paste(unique(sd10$Y[seq(5, n, 9)]), collapse = ","))
pri <- E$brms_flat_prior(bb, data = d10, family = stats::gaussian())
code <- brms::make_stancode(bb, data = d10, prior = pri)
sd10 <- E$brms_standata(bb, data = d10, prior = pri)
sf <- suppressMessages(rstan::sampling(E$brms_stan_model(code), data = sd10, chains = 0))
pts <- list(b = array(0.5, 1L), Intercept = 0, sigma = 1, Ymi = array(rep(0.2, sd10$Nmi)))
up <- rstan::unconstrain_pars(sf, pts)
pr("brms mi()+cens() log_prob at a finite point: %s", format(rstan::log_prob(sf, up, adjust_transform = FALSE)))
g <- tr(rstan::grad_log_prob(sf, up, adjust_transform = FALSE))
pr("  grad finite: %s", if (is_err(g)) msg(g, 100) else any(is.finite(g)))
s1 <- tr(suppressWarnings(rstan::sampling(E$brms_stan_model(code), data = sd10, chains = 1, iter = 200, refresh = 0, seed = 1)))
pr("  can brms sample it? %s", msg(s1, 200))
b1 <- tr(suppressWarnings(brms::brm(bb, data = d10, chains = 1, iter = 100, refresh = 0, seed = 1, silent = 2, backend = "rstan")))
pr("  brm(): %s", msg(b1, 200))
b2 <- tr(brms::make_stancode(brms::bf(y | mi() + cens(ce) ~ x), data = d10[!is.na(d10$y), ]))
d11 <- d10; d11$y[is.na(d11$y)] <- 0
b3 <- tr(brms::make_standata(brms::bf(y | mi() + cens(ce) ~ x), data = d11))
pr("  no NA rows: stancode %s, Nmi %s", msg(b2, 60), if (is_err(b3)) msg(b3) else b3$Nmi %||% "none")
