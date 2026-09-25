source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-helpers.R")
# Above the 0.05 slope threshold, where the default does NOT engage: is
# the plain fit short against autoscale = TRUE? Design of r2-m1-battery.R
# (uniform column, 25 groups x 12), seeds 101..108.
res <- NULL
for (fam in c("gaussian", "poisson", "bernoulli")) for (sx in c(0.02, 0.06, 0.08, 0.15))
for (ssd in c(0, 0.3)) for (seed in 101:108) {
  set.seed(seed)
  ng <- 25; per <- 12; n <- ng * per
  g <- factor(rep(seq_len(ng), each = per))
  xs <- runif(n, 0, 1); x <- xs * sx / sd(xs)
  eta <- 0.2 + 0.5 * xs + rnorm(ng, 0, 0.6)[g] + rnorm(ng, 0, ssd)[g] * xs
  y <- switch(fam, gaussian = eta + rnorm(n), poisson = rpois(n, exp(eta)),
              bernoulli = rbinom(n, 1, plogis(eta)))
  famo <- switch(fam, gaussian = gaussian(), poisson = poisson(), bernoulli = bernoulli())
  d <- data.frame(y, x, g)
  a <- fitw(y ~ x + (1 + x | g), data = d, family = famo)
  t <- fitw(y ~ x + (1 + x | g), data = d, family = famo,
            control = frmtmb_control(autoscale = TRUE))
  res <- rbind(res, data.frame(fam, sx, ssd, seed, engaged = !is.null(a$tpl),
    short = as.numeric(logLik(t$fit)) - as.numeric(logLik(a$fit)),
    codeD = a$fit$opt$convergence, nwD = length(a$warn)))
}
agg <- aggregate(short ~ fam + sx, res, function(v) c(max = max(v), n_gt_1e3 = sum(v > 1e-3)))
print(agg)
cat("engaged by sx:\n"); print(tapply(res$engaged, res$sx, mean))
print(res[res$short > 1e-3, ])
