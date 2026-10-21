# Reviewer: frmtmb.sample paths for negbinomial2 that the lane's draws
# test does not reach: predictive_error() and pp_check(error_binned)
# must NOT refuse negbinomial2 (the polytomous guard absent), a user
# prior on sigma reaches sigma, loo() runs, and the sampled log density
# at a draw equals the density plus brms's gamma(0.4, 0.3) prior.
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R")
suppressMessages(library(frmtmb.sample))
options(mc.cores = 1)
set.seed(8)
n <- 150
dd <- data.frame(x = rnorm(n))
dd$y <- rnbinom(n, mu = exp(0.8 + 0.5 * dd$x), size = 2.5)
fit <- frm(y ~ x, family = brmsfamily("negbinomial2"), data = dd)
ds <- suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 1, iter = 400, warmup = 200, refresh = 0,
             seed = 3)))
probe("predictive_error nb2 dims", dim(predictive_error(ds)))
probe("pp_check error_binned nb2 class",
      class(pp_check(ds, type = "error_binned")))
probe("loo nb2", loo::loo(log_lik(ds))$estimates)
probe("prior_summary nb2", prior_summary(ds))
ds2 <- suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 1, iter = 400, warmup = 200, refresh = 0,
             seed = 3, prior = set_prior("gamma(50, 500)", class = "sigma"))))
probe("user prior on sigma: posterior mean of sigma, default vs gamma(50, 500)",
      c(default = mean(posterior::as_draws_df(ds)$sigma),
        tight = mean(posterior::as_draws_df(ds2)$sigma), prior_mean = 0.1))
probe("prior_summary with user sigma prior", prior_summary(ds2))
