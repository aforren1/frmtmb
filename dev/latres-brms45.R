# brms-45: under cov = FALSE on a family with no residual, brms 2.23.0's methods
# leave out the AR sum of the Stan program they summarize. At the same parameter
# values (Stan's Fixed_param), log_prob is compared with log_lik() summed over
# the rows plus the density of zerr, for cov = FALSE and, as the control, cov =
# TRUE. Seed 167. Prints every number dev/latres-findings.md quotes.
setwd("C:/Users/adf44/source/r/frmtmb-wt-latres")
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(rstan)
})
hel <- new.env(parent = asNamespace("frmtmb"))
hel$`%||%` <- function(a, b) if (is.null(a)) b else a
sys.source("dev/equiv-lib.R", envir = hel)
set.seed(167)
d <- expand.grid(t = 1:5, g = factor(1:30))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
e <- as.vector(replicate(30, arima.sim(list(ar = 0.5), 5, sd = 0.5)))
d$count <- rpois(n, exp(0.8 + 0.4 * d$x + e))
for (ft in c("count ~ x + ar(t, g)", "count ~ x + ar(t, g, cov = TRUE)")) {
  f <- as.formula(ft)
  fit <- suppressWarnings(frm(frmtmb::bf(f), family = poisson(), data = d))
  fb <- brms::bf(f)
  fam <- brms::brmsfamily("poisson")
  bfit <- suppressWarnings(hel$brms_fixed_fit(fb, fam, d, fit))
  prior <- hel$brms_flat_prior(fb, data = d, family = fam)
  sdat <- hel$brms_standata(fb, data = d, family = fam, prior = prior)
  code <- brms::make_stancode(fb, data = d, family = fam, prior = prior)
  sf <- suppressMessages(rstan::sampling(hel$brms_stan_model(code), data = sdat,
                                         chains = 0))
  pars <- hel$stan_pars_from_fit(fit, sdat, code, NULL)
  lp <- rstan::log_prob(sf, rstan::unconstrain_pars(sf, pars),
                        adjust_transform = FALSE)
  zp <- sum(dnorm(pars$zerr, log = TRUE))
  ll <- sum(brms::log_lik(bfit)[1L, ])
  cat(sprintf(paste0("%s: log_prob %.6f; sum(brms log_lik) + log N(zerr) %.6f; ",
                     "gap %.6g; sum(frmtmb rows) + log N(zerr) %.6f; gap %.3g\n"),
              ft, lp, ll + zp, lp - ll - zp,
              sum(hel$frm_row_loglik(fit)) + zp,
              lp - sum(hel$frm_row_loglik(fit)) - zp))
}
