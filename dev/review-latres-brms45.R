# Attack 4. (a) brms's posterior_epred() and log_lik() under cov = FALSE against
# its own Stan program, at one parameter point, read by hand. (b) cov = TRUE on
# new rows: fitted() at err = 0 against brms's epred mean with a large sderr on
# a log link, and predict().
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
suppressPackageStartupMessages({library(rstan); library(brms)})
hel <- new.env(parent = asNamespace("frmtmb"))
hel$`%||%` <- `%||%`
sys.source("dev/equiv-lib.R", envir = hel)
stan_model_of <- hel$brms_stan_model

fixed_fit <- function(form, d, fam, pars, prior) {
  code <- brms::stancode(form, d, family = fam, prior = prior)
  sdat <- brms::standata(form, d, family = fam, prior = prior)
  m <- stan_model_of(code)
  sf <- rstan::sampling(m, data = sdat, algorithm = "Fixed_param", chains = 1,
                        iter = 1, warmup = 0, init = list(pars), refresh = 0)
  scaf <- brm(form, d, family = fam, prior = prior, empty = TRUE)
  scaf$fit <- sf
  scaf <- brms::rename_pars(scaf)
  list(fit = scaf, sf = sf, sdat = sdat, code = code, m = m)
}
set.seed(167)
d <- expand.grid(t = 1:5, g = factor(1:30))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
d$y <- rpois(n, exp(0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)))
d$pos <- rgamma(n, 4, 4 / exp(0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)))
zz <- rnorm(n)
pars_of <- function(sh = NULL) {
  p <- list(b = array(0.4, 1), Intercept = 0.5 + 0.4 * mean(d$x), ar = array(0.6, 1),
            zerr = zz, sderr = 0.5)
  if (!is.null(sh)) p$shape <- sh
  p
}
pri <- c(set_prior("", class = "sderr"), set_prior("", class = "Intercept"))

cat("== (a) poisson cov = FALSE\n")
ff <- fixed_fit(y ~ x + ar(t, g), d, poisson(), pars_of(), pri)
up <- rstan::unconstrain_pars(ff$sf, pars_of())
lp <- as.numeric(rstan::log_prob(ff$sf, up, adjust_transform = FALSE))
ll <- log_lik(ff$fit)
ep <- posterior_epred(ff$fit)
cat("log_prob", lp, " sum(log_lik) + log N(zerr)", sum(ll) + sum(dnorm(zz, log = TRUE)), " gap", lp - sum(ll) - sum(dnorm(zz, log = TRUE)), "\n")
# the program's mean, by hand from its text
err <- 0.5 * zz
lag <- c(0, err[-n])
lag[d$t == 1] <- 0
eta_b <- 0.5 + 0.4 * d$x + err
eta_m <- eta_b + 0.6 * lag
cat("hand sum dpois with program's mean", sum(dpois(d$y, exp(eta_m), log = TRUE)) + sum(dnorm(zz, log = TRUE)), "\n")
cat("hand sum dpois without AR sum     ", sum(dpois(d$y, exp(eta_b), log = TRUE)) + sum(dnorm(zz, log = TRUE)), "\n")
cat("brms epred vs exp(eta without AR): max rel", max(abs(ep[1, ] / exp(eta_b) - 1)), "\n")
cat("brms epred vs exp(eta with AR)   : max rel", max(abs(ep[1, ] / exp(eta_m) - 1)), "\n")
cat("posterior_linpred(incl_autocor = FALSE) vs eta without err: max", max(abs(posterior_linpred(ff$fit, incl_autocor = FALSE)[1, ] - (0.5 + 0.4 * d$x))), "\n")
cat("predict() mean under brms uses mean:", mean(exp(eta_b)), "; the program's mean", mean(exp(eta_m)), "\n")
print(tail(strsplit(ff$code, "\n")[[1]][grep("Err|err", strsplit(ff$code, "\n")[[1]])], 8))

cat("\n== (a2) Gamma log cov = FALSE\n")
ffg <- fixed_fit(pos ~ x + ar(t, g), d, Gamma(link = "log"), pars_of(4),
                 c(pri, set_prior("", class = "shape")))
upg <- rstan::unconstrain_pars(ffg$sf, pars_of(4))
lpg <- as.numeric(rstan::log_prob(ffg$sf, upg, adjust_transform = FALSE))
llg <- log_lik(ffg$fit)
epg <- posterior_epred(ffg$fit)
cat("log_prob", lpg, " sum(log_lik)+logN", sum(llg) + sum(dnorm(zz, log = TRUE)), " gap", lpg - sum(llg) - sum(dnorm(zz, log = TRUE)), "\n")
mu_m <- exp(eta_b) + 0.6 * lag
cat("brms epred vs exp(eta_b): max rel", max(abs(epg[1, ] / exp(eta_b) - 1)), "; vs program's mean: max rel", max(abs(epg[1, ] / mu_m - 1)), "\n")

cat("\n== control: cov = TRUE on poisson\n")
fft <- fixed_fit(y ~ x + ar(t, g, cov = TRUE), d, poisson(), pars_of(), pri)
upt <- rstan::unconstrain_pars(fft$sf, pars_of())
lpt <- as.numeric(rstan::log_prob(fft$sf, upt, adjust_transform = FALSE))
cat("gap", lpt - sum(log_lik(fft$fit)) - sum(dnorm(zz, log = TRUE)), "\n")

cat("\n== frm() rows under cov = FALSE at the same point\n")
f <- frm(y ~ x + ar(t, g), d, family = poisson())
bk <- fi("frame_latres")(f$frame)[[1]]
f2 <- f
f2$estimates$beta[] <- c(0.5, 0.4)
f2$estimates$theta[] <- c(log(0.5), qlogis((0.6 + 1) / 2))
f2$estimates$b[bk$b_idx] <- zz
fm <- fitted(f2)[, 1]
cat("frm fitted vs program's mean: max rel", max(abs(fm / exp(eta_m) - 1)), "\n")
cat("frm fitted vs brms epred: max rel", max(abs(fm / ep[1, ] - 1)), " mean(frm/brms)", mean(fm / ep[1, ]), "\n")

cat("\n== (b) new rows: err = 0 against brms's mean, sderr large\n")
bgrid <- list()
for (s in c(0.25, 0.5, 1, 1.5, 2)) {
  pp <- list(b = array(0.4, 1), Intercept = 0.5 + 0.4 * mean(d$x), ar = array(0.5, 1), zerr = zz, sderr = s)
  ft <- fixed_fit(y ~ x + ar(t, g, cov = TRUE), d, poisson(), pp, pri)
  set.seed(1)
  nd <- data.frame(t = 1:5, g = factor(rep("new", 5)), x = 0)
  # many fixed-point copies: the brms mean over its draws of err
  draws <- replicate(4000, as.numeric(posterior_epred(ft$fit, newdata = nd, allow_new_levels = TRUE)))
  ar_m <- apply(draws, 1, mean)
  fr <- f <- frm(y ~ x + ar(t, g, cov = TRUE), d, family = poisson())
  fr$estimates$beta[] <- c(0.5, 0.4)
  fr$estimates$theta[] <- c(log(s), 0.5 / sqrt(1 - 0.25))
  fv <- fitted(fr, newdata = nd, allow_new_levels = TRUE)[, 1]
  set.seed(2)
  sims <- simulate(fr, nsim = 4000, newdata = nd, allow_new_levels = TRUE, propagate_error = FALSE)
  pm <- colMeans(t(as.matrix(sims)))
  msd <- s / sqrt(1 - 0.25)
  bgrid[[length(bgrid) + 1]] <- data.frame(sderr = s, marg_sd = msd, exp_half = exp(msd^2 / 2),
      brms_epred = mean(ar_m), frm_fitted = mean(fv), ratio = mean(ar_m) / mean(fv),
      frm_simulate_mean = mean(pm), truth = exp(0.5) * exp(msd^2 / 2))
  print(bgrid[[length(bgrid)]], digits = 4)
}
print(do.call(rbind, bgrid), digits = 4)
