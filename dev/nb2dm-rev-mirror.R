# Reviewer: negbinomial vs negbinomial2 at the Poisson limit (seeds 4 and
# 6 of dev/nb2dm-rev-runoff.R): same likelihood, mirrored parameter.
# Where do the two fits part?
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R")
for (seed in c(4, 6, 1)) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = rnorm(n))
  d$y <- rpois(n, exp(1 + 0.3 * d$x))
  f1 <- suppressWarnings(frm(y ~ x, family = negbinomial(), data = d))
  f2 <- suppressWarnings(frm(y ~ x, family = brmsfamily("negbinomial2"),
                             data = d))
  cat(sprintf("seed %d: logLik nb %.12f nb2 %.12f\n", seed, logLik(f1),
              logLik(f2)))
  cat("  start nb :", format(f1$obj$par, digits = 10), "\n")
  cat("  start nb2:", format(f2$obj$par, digits = 10), "\n")
  cat("  opt nb   :", format(f1$opt$par, digits = 10), " iters",
      f1$opt$iterations, f1$opt$message, "\n")
  cat("  opt nb2  :", format(f2$opt$par, digits = 10), " iters",
      f2$opt$iterations, f2$opt$message, "\n")
  # the two objectives at mirrored points: identical functions?
  p <- f1$opt$par
  q <- p
  q[length(q)] <- -q[length(q)]
  cat(sprintf("  fn nb(p) %.15g  fn nb2(mirror p) %.15g\n", f1$obj$fn(p),
              f2$obj$fn(q)))
  g1 <- f1$obj$gr(p)
  g2 <- f2$obj$gr(q)
  cat("  gr nb(p):", format(g1, digits = 6), " gr nb2(mirror):",
      format(g2, digits = 6), "\n")
  cat("  lower/upper nb :", format(f1$opt$lower %||% NA), "/",
      format(f1$opt$upper %||% NA), "\n")
}
