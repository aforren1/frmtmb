# Item 3 against brms 2.23.0, run: conditional_effects()'s se__ on a poisson
# and a lognormal mean, frmtmb's (READERS_ARM build) against brms's sampled
# fit, 4 chains of 1000 draws, seed 3. Data seed 7971, as test-predict-draw-
# mean.R builds it.
source("dev/readers-setup.R")
suppressPackageStartupMessages(library(brms))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
set.seed(7971)
n <- 300
d <- data.frame(x = rnorm(n))
eta <- 0.4 + 0.7 * d$x
d$cnt <- rpois(n, exp(eta))
d$pos <- exp(eta + rnorm(n, 0, exp(-0.5 + 0.3 * d$x)))
ic <- list(x = c(-1, 0, 1))
for (sp in list(list(frmtmb::bf(cnt ~ x), poisson(), brms::bf(cnt ~ x),
                     "poisson"),
                list(frmtmb::bf(pos ~ x, sigma ~ x), frmtmb::lognormal(),
                     brms::bf(pos ~ x, sigma ~ x), "lognormal"))) {
  f <- frm(sp[[1]], family = sp[[2]], data = d)
  b <- brm(sp[[3]], family = sp[[4]], data = d, chains = 4, iter = 2000,
           seed = 3, refresh = 0, silent = 2)
  cf <- frmtmb::conditional_effects(f, effects = "x",
                                    int_conditions = ic)[[1]]
  cb <- brms::conditional_effects(b, effects = "x", int_conditions = ic)[[1]]
  cat("\n==", sp[[4]], "\n")
  print(data.frame(x = cf$x, frmtmb_est = cf$estimate__, frmtmb_se = cf$se__,
                   brms_est = cb$estimate__, brms_se = cb$se__,
                   ratio = cf$se__ / cb$se__))
}
