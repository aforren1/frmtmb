# Review of lane nlstart: a prior on a coefficient nl_par_units() scales.
# The unit is read from the body alone, so a prior much tighter than the
# data puts its curvature times unit^2 into the optimizer's coordinates.
# Each row: body `a * xs` (xs = x * s) or `a + b * x` (y at scale s), a
# normal(m, sd) prior, against the linear spelling with the same prior.
# Usage: REV_ARM=after|before Rscript dev/review-nlstart-prior.R
source("dev/review-nlstart-common.R")
cat("ROWHEAD case s prior_sd seed ll_body conv_body nwarn ll_lin conv_lin",
    "gap\n")
for (seed in c(6201, 6202, 6203)) {
  set.seed(seed)
  n <- 200
  d <- data.frame(x = rnorm(n))
  d$y <- 0.5 * d$x + rnorm(n)
  for (s in c(1e-2, 1e-4, 1e-8, 1e-20)) {
    d$xs <- d$x * s
    for (psd in c(1e-3, 1, 1e3)) {
      pr_b <- set_prior(sprintf("normal(0, %g)", psd), nlpar = "a")
      pr_l <- set_prior(sprintf("normal(0, %g)", psd), class = "b")
      fb <- tt(frm(bf(y ~ a * xs, a ~ 1, nl = TRUE), data = d, prior = pr_b))
      fl <- tt(frm(y ~ 0 + xs, data = d, prior = pr_l))
      cat(sprintf("ROW ax %g %g %d %.8f %s %d %.8f %s %.3g\n", s, psd, seed,
                  llv(fb), if (iserr(fb)) NA else fb$opt$convergence,
                  length(attr(fb, "w")), llv(fl),
                  if (iserr(fl)) NA else fl$opt$convergence,
                  llv(fl) - llv(fb)))
    }
  }
  for (s in c(1e2, 1e4, 1e6)) {
    d$yl <- s * (1 + 0.05 * d$x) + rnorm(n, 0, s * 1e-3)
    for (psd in c(1e-3, 1, 1e3)) {
      pr_b <- c(set_prior(sprintf("normal(0, %g)", psd), nlpar = "b"),
                set_prior(sprintf("normal(%g, %g)", s, s), nlpar = "a"))
      pr_l <- c(set_prior(sprintf("normal(0, %g)", psd), class = "b"),
                set_prior(sprintf("normal(%g, %g)", s, s),
                          class = "Intercept"))
      fb <- tt(frm(bf(yl ~ a + b * x, a + b ~ 1, nl = TRUE), data = d,
                   prior = pr_b))
      fl <- tt(frm(yl ~ x, data = d, prior = pr_l))
      cat(sprintf("ROW apbx %g %g %d %.8f %s %d %.8f %s %.3g\n", s, psd, seed,
                  llv(fb), if (iserr(fb)) NA else fb$opt$convergence,
                  length(attr(fb, "w")), llv(fl),
                  if (iserr(fl)) NA else fl$opt$convergence,
                  llv(fl) - llv(fb)))
    }
  }
  flush.console()
}
cat("DONE\n")
