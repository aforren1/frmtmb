# Items 2, 5: other links under a chart; prior = "flat" and class sds on the
# sampling route; a smooth in an extension's reader; the frm route with a bound
# on sds; sds with a bound box. Seed 10700.
source("dev/review-refuse-lib.R")
library(frmtmb.sample)
set.seed(10700)
n <- 150; x <- rnorm(n)
show <- function(lab, r) {
  if (is_err(r)) { pr("%-34s ERROR %s", lab, msg(r, 170)); return(invisible()) }
  pr("%-34s verdict=%s const=%s pts=[%s] grad=%s frm[%s] l1[%s]", lab, r$verdict, format(r$l1_const %||% NA, digits = 3),
     paste(format(unlist(r$l1_extra_const) %||% NA, digits = 3), collapse = ","), format(r$l1_grad %||% NA, digits = 3),
     substr(r$frm_error %||% "", 1, 120), substr(r$l1_error %||% "", 1, 110))
}
dg <- data.frame(x = x, y = rgamma(n, 2, 2 / exp(0.3 + 0.2 * x))); dg$y[seq(4, n, 12)] <- NA
show("Gamma(identity)", tr(meas(list(bf = "bf(y | mi() ~ x)", family = "Gamma(link = 'identity')", data = dg))))
show("Gamma(inverse)", tr(meas(list(bf = "bf(y | mi() ~ x)", family = "Gamma(link = 'inverse')", data = dg))))
show("inverse.gaussian default link", tr(meas(list(bf = "bf(y | mi() ~ x)", family = "inverse.gaussian()", data = dg))))
show("exponential(identity)", tr(meas(list(bf = "bf(y | mi() ~ x)", family = "exponential(link = 'identity')", data = dg))))
show("weibull(identity)", tr(meas(list(bf = "bf(y | mi() ~ x)", family = "weibull(link = 'identity')", data = dg))))
db <- data.frame(x = x, y = rbeta(n, plogis(0.2 + 0.3 * x) * 5, (1 - plogis(0.2 + 0.3 * x)) * 5)); db$y[seq(4, n, 12)] <- NA
show("Beta(probit)", tr(meas(list(bf = "bf(y | mi() ~ x)", family = "Beta(link = 'probit')", data = db))))
show("Gamma, shape ~ x", tr(meas(list(bf = "bf(y | mi() ~ x, shape ~ x)", family = "Gamma(link = 'log')", data = dg))))
dg$w <- 1 + rbinom(n, 2, 0.5)
show("Gamma, weights", tr(meas(list(bf = "bf(y | mi() + weights(w) ~ x)", family = "Gamma(link = 'log')", data = dg))))
dg$g <- factor(rep(1:15, 10))
show("Gamma, (1 | g)", tr(meas(list(bf = "bf(y | mi() ~ x + (1 | g))", family = "Gamma(link = 'log')", data = dg))))
dg$z <- 1 + 0.4 * log(ifelse(is.na(dg$y), 1, dg$y)) + rnorm(n)
show("two-way: Gamma read, z mi too", tr(meas(list(bf = "mvbf(bf(y | mi() ~ x, family = Gamma(link = 'log')), bf(z | mi() ~ mi(y)), rescor = FALSE)", family = NA,
  data = transform(dg, z = ifelse(seq_len(n) %% 17 == 0, NA, z))))))
dg2 <- dg; dg2$s <- 0.2
dg2$y[seq(4, n, 12)] <- 1
show("Gamma mi(s)", tr(meas(list(bf = "bf(y | mi(s) ~ x)", family = "Gamma(link = 'log')", data = dg2))))
fbr <- brms::Beta()
a <- tr(frm(bf(y | mi() ~ x), family = brms::brmsfamily("beta"), data = db)); pr("brmsfamily('beta') object: %s", msg(a, 120))
a <- tr(frm(bf(y | mi() ~ x), family = "beta", data = db)); pr("family = 'beta' string: %s", msg(a, 120))

# class sds and prior = "flat" on the sampling route
set.seed(10701)
n <- 150
d <- data.frame(x = runif(n), z = runif(n)); d$y <- sin(6 * d$x) + rnorm(n, 0, 0.4)
fb <- brms::bf(y ~ s(x, k = 5))
pri <- E$brms_flat_prior(fb, data = d, family = stats::gaussian())
code <- brms::make_stancode(fb, data = d, prior = pri)
sdat <- E$brms_standata(fb, data = d, prior = pri)
sb <- suppressWarnings(rstan::sampling(E$brms_stan_model(code), data = sdat, chains = 4, iter = 2000, refresh = 0, seed = 11, cores = 1))
sf <- tr(frm_sample(bf(y ~ s(x, k = 5)), data = d, prior = "flat", chains = 4, iter = 2000, seed = 11, refresh = 0, cores = 1))
if (is_err(sf)) pr("frm_sample flat: %s", msg(sf, 200)) else {
  eb <- rstan::extract(sb); dm <- posterior::as_draws_df(sf)
  q <- c(0.1, 0.5, 0.9)
  pr("flat sds quantiles: brms %s; frm_sample %s", paste(round(quantile(eb$sds_1, q), 3), collapse = " "), paste(round(quantile(dm[["sds_sx_1"]], q), 3), collapse = " "))
  pr("flat sigma mean: brms %.4f; frm_sample %.4f", mean(eb$sigma), mean(dm[["sigma"]]))
}
# set_prior with a bound on class sds, on each route
a <- tr(frm(bf(y ~ s(x, k = 5)), data = d, prior = set_prior("normal(0, 1)", class = "sds", lb = 0.1, ub = 1)))
pr("frm sds with lb, ub: %s; sd %s", msg(a, 150), if (is_err(a)) "" else format(VarCorr(a)[[1]]$sd[1], digits = 4))
a <- tr(frm_sample(bf(y ~ s(x, k = 5)), data = d, prior = set_prior("normal(0, 1)", class = "sds", lb = 0.1, ub = 1), chains = 1, iter = 300, seed = 2, refresh = 0))
if (is_err(a)) pr("frm_sample sds bounded: %s", msg(a, 200)) else {
  v <- posterior::as_draws_df(a)[["sds_sx_1"]]
  pr("frm_sample sds with lb, ub: range %.3f %.3f", min(v), max(v))
}
a <- tr(frm_sample(bf(y ~ s(x, k = 5)), data = d, prior = set_prior("normal(0, 1)", class = "sds", coef = "s(x, k = 5)", lb = 0.1), chains = 1, iter = 100, seed = 2, refresh = 0))
pr("frm_sample sds coef with lb (brms refuses): %s", msg(a, 170))
# mixed draws of a smooth prior_draws, sample_prior only
a <- tr(frm_sample(bf(y ~ s(x, k = 5)), data = d, prior = set_prior("exponential(2)", class = "sds") + set_prior("normal(0, 2)", class = "b") + set_prior("exponential(1)", class = "sigma") + set_prior("normal(0, 2)", class = "Intercept"), chains = 1, iter = 500, seed = 2, refresh = 0, sample_prior = "only"))
pr("sample_prior = only: %s; sds mean %s", msg(a, 200), if (is_err(a)) "" else format(mean(posterior::as_draws_df(a)[["sds_sx_1"]]), digits = 3))
# pairs with a gp and a car: class sd still reaches a car, as NEWS leaves it
# frm_simulate: its real signature
print(names(formals(frm_simulate)))
a <- tr(frm_simulate(bf(y ~ s(x, k = 5)), data = d, prior = set_prior("exponential(2)", class = "sds") + set_prior("normal(0, 1)", class = "b") + set_prior("normal(0, 1)", class = "Intercept") + set_prior("exponential(1)", class = "sigma"), nsim = 2, seed = 3))
pr("frm_simulate, sds prior: %s", msg(a, 240))
a0 <- tr(frm_simulate(bf(y ~ s(x, k = 5)), data = d, prior = set_prior("exponential(2)", class = "sd") + set_prior("normal(0, 1)", class = "b") + set_prior("normal(0, 1)", class = "Intercept") + set_prior("exponential(1)", class = "sigma"), nsim = 2, seed = 3))
pr("frm_simulate, sd prior on a smooth: %s", msg(a0, 240))
