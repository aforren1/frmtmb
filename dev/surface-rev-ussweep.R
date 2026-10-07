# Reviewer of lane surface, claim 9: does frm_sample(fit) on a
# correlated (1 + x | g) block, which the fit route samples CENTERED
# (flat priors), meet the +Inf the lane fixed for gp()? The us block's
# log density at log sd -1137.64 is +Inf on both builds
# (dev/surface-rev-otherblocks.R).
#   Rscript dev/surface-rev-ussweep.R <arm> <structure> <seed>
a <- commandArgs(TRUE)
source("dev/surface-rev-env.R")
rev_env(a[1])
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(5)
n_g <- 12
d <- data.frame(g = factor(rep(seq_len(n_g), each = 5)), x = rnorm(n_g * 5))
d$y <- 1 + 0.5 * d$x + rnorm(n_g, 0, 0.3)[d$g] + rnorm(n_g * 5)
f <- switch(a[2],
  us = bf(y ~ x + (1 + x | g)),
  gp = { d$xx <- round(runif(nrow(d), 0, 6), 1)
         d$y <- 0.5 + sin(d$xx) + rnorm(nrow(d), 0, 0.3)
         bf(y ~ gp(xx)) })
fit <- suppressWarnings(suppressMessages(frm(f, family = gaussian(),
                                             data = d)))
ds <- suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 1, iter = 600, refresh = 0,
             seed = as.integer(a[3]))))
sp <- rstan::get_sampler_params(ds$stanfit, inc_warmup = FALSE)[[1]]
dr <- as.matrix(ds)
cat(sprintf("%s %s seed %s: accept %.3f, stepsize %.3g, divergent %d, distinct draws of col 1 %d, min theta1 %.1f\n",
            a[1], a[2], a[3], mean(sp[, "accept_stat__"]),
            sp[1, "stepsize__"], sum(sp[, "divergent__"]),
            length(unique(dr[, 1])),
            min(rstan::extract(ds$stanfit, permuted = FALSE,
                               inc_warmup = TRUE)[, , grep("theta",
              dimnames(rstan::extract(ds$stanfit, permuted = FALSE))[[3]])[1]])))
