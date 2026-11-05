# Item 2 and 3: the sampling route. brms and frm_sample() on one mi() model per
# chart, with the complete-case posterior as the reference: an unread response
# has the complete-case posterior for every other parameter. The draws a user
# reads: miss_<k> against the family's own predictive at the fitted parameters.
# Usage: Rscript dev/review-refuse-sample.R <case>. Seed 9700.
case <- commandArgs(trailingOnly = TRUE)[1L]
source("dev/review-refuse-lib.R")
library(frmtmb.sample)
set.seed(9700)
n <- 120; x <- rnorm(n); eta <- 0.3 + 0.5 * x
cfg <- switch(case,
  Gamma = list(y = rgamma(n, 1.5, 1.5 / exp(eta)), f = stats::Gamma(link = "log"),
               fb = brms::brmsfamily("gamma", "log"), dp = "shape", inv = exp),
  beta = list(y = rbeta(n, plogis(eta) * 3, (1 - plogis(eta)) * 3), f = Beta(),
              fb = brms::Beta(), dp = "phi", inv = plogis),
  von_mises = list(y = ((2 * atan(eta) + rnorm(n, 0, 0.6) + pi) %% (2 * pi)) - pi,
                   f = von_mises(), fb = brms::von_mises(), dp = "kappa",
                   inv = function(u) 2 * pi * plogis(u) - pi),
  student = list(y = eta + rt(n, 3), f = student(), fb = brms::student(), dp = "nu",
                 inv = identity))
d <- data.frame(x = x, y = cfg$y)
mrows <- seq(5, n, by = 8)
d$y[mrows] <- NA
dc <- d[-mrows, ]
one_brms <- function(form, dat) {
  pri <- brms::default_prior(form, data = dat, family = cfg$fb)
  code <- brms::make_stancode(form, data = dat, family = cfg$fb, prior = pri)
  sdat <- E$brms_standata(form, data = dat, family = cfg$fb, prior = pri)
  rstan::sampling(E$brms_stan_model(code), data = sdat, chains = 4, iter = 2000,
                  seed = 11, refresh = 0, cores = 1)
}
fj <- suppressWarnings(one_brms(brms::bf(y | mi() ~ x), d))
fc <- suppressWarnings(one_brms(brms::bf(y ~ x), dc))
ej <- rstan::extract(fj); ec <- rstan::extract(fc)
sum2 <- function(v) sprintf("%.4f (%.4f)", mean(v), sd(v))
pr("%s brms joint: b %s dp %s Int %s", case, sum2(ej$b), sum2(ej[[cfg$dp]]), sum2(ej$Intercept))
pr("%s brms cc   : b %s dp %s Int %s", case, sum2(ec$b), sum2(ec[[cfg$dp]]), sum2(ec$Intercept))
# frmtmb Laplace, both routes
a <- frm(bf(y | mi() ~ x), family = cfg$f, data = d)
b <- frm(bf(y ~ x), family = cfg$f, data = dc)
pr("%s frm joint : b %.4f dp(log) %.4f (exp %.4f)", case, a$estimates$beta[[2]], a$estimates$betad[[1]], exp(a$estimates$betad[[1]]))
pr("%s frm cc    : b %.4f dp(log) %.4f (exp %.4f)", case, b$estimates$beta[[2]], b$estimates$betad[[1]], exp(b$estimates$betad[[1]]))
# frm_sample on the joint, with brms's default priors
s <- tr(frm_sample(bf(y | mi() ~ x), family = cfg$f, data = d, chains = 4,
                   iter = 2000, seed = 11, refresh = 0, cores = 1))
if (is_err(s)) {
  pr("frm_sample: %s", msg(s, 300))
  q("no")
}
dm <- posterior::as_draws_df(s)
cn <- names(dm)
pr("%s draw names: %s", case, paste(head(cn, 12), collapse = " "))
bn <- grep("^b_x$", cn, value = TRUE)
dn <- cn[grepl(paste0("^", cfg$dp, "$|^b_", cfg$dp, "_Intercept$"), cn)]
pr("%s frm_sample: b %s dp %s", case, sum2(dm[[bn]]), if (length(dn)) sum2(dm[[dn[1]]]) else "?")
mis <- grep("^miss_", cn, value = TRUE)
pr("%s miss_1 draw range [%.3f, %.3f]; brms Ymi[1] range [%.3f, %.3f]", case,
   min(dm[[mis[1]]]), max(dm[[mis[1]]]), min(ej$Ymi[, 1]), max(ej$Ymi[, 1]))
vs <- cfg$inv(as.matrix(as.data.frame(dm)[, mis]))
if (identical(case, "student")) vs <- as.matrix(as.data.frame(dm)[, mis])
pr("%s value mean (frm_sample via chart) %s; brms Ymi mean %s",
   case, paste(sprintf("%.3f", colMeans(vs)[1:4]), collapse = " "),
   paste(sprintf("%.3f", colMeans(ej$Ymi)[1:4]), collapse = " "))
pr("%s value sd (frm_sample via chart) %s; brms Ymi sd %s",
   case, paste(sprintf("%.3f", apply(vs, 2, sd)[1:4]), collapse = " "),
   paste(sprintf("%.3f", apply(ej$Ymi, 2, sd)[1:4]), collapse = " "))
ks <- vapply(seq_len(ncol(vs)), function(j) {
  suppressWarnings(ks.test(vs[, j], ej$Ymi[, j])$statistic)
}, 0)
pr("%s KS(frm_sample value, brms Ymi): median %.3f max %.3f on %d values", case, median(ks), max(ks), length(ks))
pe <- tr(posterior_epred(s, ndraws = 50))
pr("%s posterior_epred: %s", case, msg(pe, 90))
pp <- tr(posterior_predict(s, ndraws = 50))
pr("%s posterior_predict: %s NA at missing %s", case, msg(pp, 90), if (is_err(pp)) "-" else sum(is.na(pp[, mrows])))
ll <- tr(log_lik(s, ndraws = 50))
pr("%s log_lik: %s", case, msg(ll, 120))
ps <- tr(posterior::summarise_draws(posterior::as_draws_df(s)))
ss <- tr(summary(s))
pr("%s summary: %s", case, msg(ss, 120))
