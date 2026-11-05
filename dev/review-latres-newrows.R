# Attack 4 (b): new rows under cov = TRUE. brms's epred reads one draw of err per
# posterior draw; held parameters by Fixed_param, 4000 copies of one point, so its
# mean over draws is its expected value over the latent. fitted() is at err = 0.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
suppressPackageStartupMessages({library(rstan); library(brms)})
hel <- new.env(parent = asNamespace("frmtmb"))
hel$`%||%` <- `%||%`
sys.source("dev/equiv-lib.R", envir = hel)
set.seed(167)
d <- expand.grid(t = 1:5, g = factor(1:30))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
d$y <- rpois(n, exp(0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)))
d$bin <- rbinom(n, 1, plogis(0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)))
zz <- rnorm(n)
pri <- c(set_prior("", class = "sderr"), set_prior("", class = "Intercept"))
form <- y ~ x + ar(t, g, cov = TRUE)
code <- brms::stancode(form, d, family = poisson(), prior = pri)
sdat <- brms::standata(form, d, family = poisson(), prior = pri)
m <- hel$brms_stan_model(code)
codeb <- brms::stancode(bin ~ x + ar(t, g, cov = TRUE), d, family = bernoulli(), prior = pri)
sdatb <- brms::standata(bin ~ x + ar(t, g, cov = TRUE), d, family = bernoulli(), prior = pri)
mb <- hel$brms_stan_model(codeb)
fit_p <- frm(form, d, family = poisson())
fit_b <- frm(bin ~ x + ar(t, g, cov = TRUE), d, family = bernoulli())
nd <- data.frame(t = 1:5, g = factor(rep("new", 5)), x = 0)
out <- list()
for (s in c(0.25, 0.5, 1, 1.5, 2)) {
  pp <- list(b = array(0.4, 1), Intercept = 0.5 + 0.4 * mean(d$x), ar = array(0.5, 1), zerr = zz, sderr = s)
  for (fam in c("poisson", "bernoulli")) {
    sf <- rstan::sampling(if (fam == "poisson") m else mb, data = if (fam == "poisson") sdat else sdatb, algorithm = "Fixed_param", chains = 1,
                          iter = 4000, warmup = 0, init = list(pp), refresh = 0, seed = 1)
    scaf <- brm(if (fam == "poisson") form else bin ~ x + ar(t, g, cov = TRUE), d, family = if (fam == "poisson") brms::brmsfamily("poisson") else brms::brmsfamily("bernoulli"), prior = pri, empty = TRUE)
    scaf$fit <- sf
    bf <- brms::rename_pars(scaf)
    set.seed(3)
    ep <- brms::posterior_epred(bf, newdata = nd, allow_new_levels = TRUE)
    pr <- brms::posterior_predict(bf, newdata = nd, allow_new_levels = TRUE)
    fr <- if (fam == "poisson") fit_p else fit_b
    fr$estimates$beta[] <- c(0.5, 0.4)
    fr$estimates$theta[] <- c(log(s), 0.5 / sqrt(1 - 0.25))
    fv <- fitted(fr, newdata = nd, allow_new_levels = TRUE)[, 1]
    sim <- simulate(fr, nsim = 4000, newdata = nd, allow_new_levels = TRUE, seed = 2)
    msd <- s / sqrt(1 - 0.25)
    out[[length(out) + 1]] <- data.frame(fam = fam, sderr = s, marg_sd = msd,
      expected = if (fam == "poisson") exp(0.5 + msd^2 / 2) else integrate(function(u) plogis(0.5 + msd * u) * dnorm(u), -Inf, Inf)$value,
      brms_epred = mean(ep), frm_fitted = mean(fv), ratio_brms_frm = mean(ep) / mean(fv),
      brms_predict = mean(pr), frm_simulate = mean(as.matrix(sim)))
    print(out[[length(out)]], digits = 4)
  }
}
cat("\n")
print(do.call(rbind, out), digits = 4, row.names = FALSE)
# on the fitted fit: fitted(newdata = the same rows) against fitted() and the mean
f <- frm(form, d, family = poisson())
th <- f$estimates$theta
cat("estimate: sderr", exp(th[1]), "ar", th[2] / sqrt(1 + th[2]^2), "\n")
cat("mean fitted()", mean(fitted(f)[, 1]), " mean fitted(newdata = d)", mean(fitted(f, newdata = d)[, 1]), " mean y", mean(d$y), " mean predict(newdata = d)", mean(predict(f, newdata = d, ndraws = 4000)[, 1]), "\n")
print(body(fi("latres_refuse_newdata"))[[1]])
# new levels of a random effect, the precedent: what does fitted() do there?
fr <- frm(y ~ x + (1 | g), d, family = poisson())
fr$estimates$theta[] <- log(1)
cat("RE precedent, new level, sd = 1: fitted", fitted(fr, newdata = nd, allow_new_levels = TRUE)[1, 1], " exp(b0)", exp(fr$estimates$beta[1]), "\n")
