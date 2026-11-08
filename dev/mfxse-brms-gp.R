# brms 2.23.0's posterior SD of the average slope in x on the exact gp()
# designs of dev/mfxse-gpslope.R (data seed 8150), with the slope taken as
# dev/mfxse-gpslope.R takes it, from posterior_epred() at x +- h/2. brms
# draws the field at unseen positions jointly, so the difference carries the
# field's conditional variance as well. 2 chains of 2000, default priors.
#   Rscript dev/mfxse-brms-gp.R
source("C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-lib.R")
suppressPackageStartupMessages(library(brms))
set.seed(8150)
m <- 90
d <- data.frame(x = runif(m, 0, 6), z = runif(m, 0, 6))
d$y <- sin(d$x) + 0.3 * d$z + rnorm(m, sd = 0.3)
h <- 1e-4
for (sp in c("y ~ z + gp(x)", "y ~ gp(x, gr = FALSE)")) {
  b <- quiet(brm(as.formula(sp), data = d, chains = 2, iter = 2000,
                 cores = 1, seed = 8150, refresh = 0, backend = "rstan"))
  # the field at the fitted positions plus the slope through x +- h/2 of the
  # same draws: brms predicts new gp positions jointly per draw
  nd <- rbind(transform(d, x = x + h / 2), transform(d, x = x - h / 2))
  ep <- posterior_epred(b, newdata = nd)
  s <- rowMeans((ep[, seq_len(m)] - ep[, m + seq_len(m)]) / h)
  cat(sprintf("BRMSGP %-24s mean %.5f sd %.5f rhat %.3f\n", sp, mean(s),
              sd(s), max(rhat(b), na.rm = TRUE)))
}
