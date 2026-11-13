# Lane nlstop, item 4850: the rows of dev/review-nlstart-prior.R (63 rows:
# a * xs and a + b * x against the linear spellings, each under a normal
# prior), with the penalized log-likelihood the fits maximize, the code
# and the warnings of both spellings.
# Usage: NLSTOP_ARM=before|after Rscript dev/nlstop-prior.R
source("dev/nlstop-common.R")
cat("ROWHEAD case s prior_sd seed llp_body conv_body nw_body llp_lin",
    "conv_lin nw_lin wc_lin t_lin gap\n")
row <- function(case, s, psd, seed, fb, fl, tl) {
  cat(sprintf("ROW %s %g %g %d %.8f %s %d %.8f %s %d %s %.3f %.3g\n", case, s,
              psd, seed, llp(fb), if (iserr(fb)) NA else fb$opt$convergence,
              length(attr(fb, "w")), llp(fl),
              if (iserr(fl)) NA else fl$opt$convergence,
              length(attr(fl, "w")), wsum(fl), tl, llp(fl) - llp(fb)))
}
timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  v <- expr
  attr(v, "t") <- proc.time()[["elapsed"]] - t0
  v
}
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
      fl <- timed(tt(frm(y ~ 0 + xs, data = d, prior = pr_l)))
      row("ax", s, psd, seed, fb, fl, attr(fl, "t"))
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
      fl <- timed(tt(frm(yl ~ x, data = d, prior = pr_l)))
      row("apbx", s, psd, seed, fb, fl, attr(fl, "t"))
      if (s == 1e6 && psd == 1e3) {
        cat("  lin beta:", format(fl$estimates$beta, digits = 7),
            " body beta:", format(fb$estimates$beta, digits = 7), "\n")
        for (w in attr(fl, "w")) cat("  W lin:", substr(w, 1, 200), "\n")
        for (m in attr(fl, "m")) cat("  M lin:", substr(m, 1, 200), "\n")
      }
    }
  }
  flush.console()
}
cat("DONE\n")
