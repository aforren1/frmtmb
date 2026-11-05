# me() beside a latent residual, L1 through the instrument with brms's me() on the search path. Seed 167.
setwd("C:/Users/adf44/source/r/frmtmb-wt-latres")
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib","C:/Users/adf44/source/r/rellib-r16","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win", FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/stan-cache")
suppressPackageStartupMessages({library(frmtmb); library(rstan)})
hel <- new.env(parent = asNamespace("frmtmb")); hel$`%||%` <- function(a, b) if (is.null(a)) b else a
sys.source("dev/equiv-lib.R", envir = hel)
set.seed(167); d <- expand.grid(t = 1:5, g = factor(1:30)); d <- d[order(d$g, d$t), ]; n <- nrow(d)
d$x <- rnorm(n); d$z <- rnorm(n); d$sx <- 0.2
e <- as.vector(replicate(30, arima.sim(list(ar = 0.5), 5, sd = 0.6)))
d$count <- rpois(n, exp(0.8 + 0.4 * d$x + e))
f <- count ~ x + me(z, sx) + ar(t, g, cov = TRUE)
fit <- suppressWarnings(frm(frmtmb::bf(f), family = poisson(), data = d))
out <- hel$with_brms_me(function() {
  fb <- brms::bf(f); fam <- brms::brmsfamily("poisson")
  prior <- hel$brms_flat_prior(fb, data = d, family = fam)
  hel$eq_l1(list(data = d), fb, fam, d, fit, prior, joint = TRUE)
})
cat(sprintf("me + ar: const %.4g grad %.3g logLik %.6f\n", out$at$const, out$at$grad, logLik(fit)))
