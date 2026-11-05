# Short check: the small items my stale scripts no longer reach: m4, m7, m10, an
# all-missing response, a Gamma KS on Ymi against brms. Seed 11700.
source("dev/review-refuse-lib.R")
suppressPackageStartupMessages(library(frmtmb.sample))
quiet <- function(e) suppressMessages(suppressWarnings(e))
set.seed(11700)
n <- 150; lev <- paste0("L", 1:8); x <- rnorm(n)
dm <- data.frame(g1 = sample(lev, n, TRUE), g2 = sample(lev, n, TRUE), x = x, y = rnorm(n))
A0 <- exp(-abs(outer(1:8, 1:8, "-")) / 3); dimnames(A0) <- list(lev, lev)
r <- tr(meas(list(bf = "bf(y ~ x + (1 + x || mm(g1, g2, cov = A)))", data = dm, data2 = list(A = A0))))
pr("m4 || over mm(cov): %s const %s grad %s frm[%s] l1[%s]", r$verdict, format(r$l1_const %||% NA, digits = 3), format(r$l1_grad %||% NA, digits = 3), substr(r$frm_error %||% "", 1, 80), substr(r$l1_error %||% "", 1, 80))
d4 <- data.frame(x = x, y = runif(n, -3, 3), s = 0.3); d4$y[3] <- pi
a <- tr(frm(bf(y | mi(s) ~ x), family = von_mises(), data = d4)); pr("m7 von_mises mi(s), y = pi: %s", msg(a, 130))
fam_args <- list(family = "tfam", dpars = "mu", links = list(mu = "log"),
  lpdf = function(y, dpars, aterms) RTMB::dgamma(y, 2, dpars[["mu"]] / 2, log = TRUE))
vy <- function(y, ...) if (any(y <= 0)) stop("tfam: response must be positive") else invisible(NULL)
d <- data.frame(x = x, y = rgamma(n, 2, 2 / exp(0.3 + 0.5 * x))); d$y[seq(5, n, 10)] <- NA
a <- tr(frm(bf(y | mi() ~ x), family = do.call(frmtmb_family, c(fam_args, list(accepts_aterms = "mi", valid_y = vy))), data = d))
pr("m10 custom positive, no chart: %s", msg(a, 130))
d2 <- d; d2$y <- NA_real_
a <- tr(frm(bf(y | mi() ~ x), family = stats::Gamma(link = "log"), data = d2)); pr("all missing, Gamma: %s", msg(a, 130))
a <- tr(frm(bf(y | mi() ~ x), data = d2)); pr("all missing, gaussian: %s", msg(a, 130))
# Ymi against brms, Gamma
fb <- brms::bf(y | mi() ~ x); fam <- brms::brmsfamily("gamma", "log")
pri <- brms::default_prior(fb, data = d, family = fam)
code <- brms::make_stancode(fb, data = d, family = fam, prior = pri)
sdat <- E$brms_standata(fb, data = d, family = fam, prior = pri)
sb <- quiet(rstan::sampling(E$brms_stan_model(code), data = sdat, chains = 4, iter = 2000, seed = 11, refresh = 0, cores = 1))
eb <- rstan::extract(sb)
s <- quiet(frm_sample(bf(y | mi() ~ x), family = stats::Gamma(link = "log"), data = d, chains = 4, iter = 2000, seed = 11, refresh = 0, cores = 1))
dd <- as.data.frame(posterior::as_draws_df(s))
nn <- paste0("Ymi[", seq_len(ncol(eb$Ymi)), "]")
ks <- vapply(seq_along(nn), function(j) suppressWarnings(ks.test(dd[[nn[j]]], eb$Ymi[, j])$statistic), 0)
pr("Ymi against brms, by index: KS median %.3f max %.3f on %d; range frm %.3f %.3f, brms %.3f %.3f", median(ks), max(ks), length(ks), min(dd[nn]), max(dd[nn]), min(eb$Ymi), max(eb$Ymi))
s2 <- quiet(frm_sample(frm(bf(y | mi() ~ x), family = stats::Gamma(link = "log"), data = d), chains = 4, iter = 2000, seed = 11, refresh = 0, cores = 1))
dd2 <- as.data.frame(posterior::as_draws_df(s2))
ks2 <- vapply(seq_along(nn), function(j) suppressWarnings(ks.test(dd2[[nn[j]]], eb$Ymi[, j])$statistic), 0)
pr("frm_sample(fit) by index: KS median %.3f max %.3f", median(ks2), max(ks2))
