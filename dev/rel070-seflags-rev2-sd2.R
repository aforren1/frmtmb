# Reviewer of lane seflags, claim 2: the barely-read sd of
# dev/seflags-rev-sd.R (y ~ a + 1e-4 * u, u ~ 0 + (1 | g), no group
# effect in the data). The SE the fit reports for theta_1, its Hessian
# row, and the profile of the objective in theta_1 (others re-optimized
# is too costly here; others held) from -12 to +12.
#   Rscript dev/seflags-rev2-sd2.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-release/dev/rel070-seflags-rev2-common.R")
for (seed in 1:3) for (G in c(3, 5, 10)) {
  set.seed(seed)
  n <- 20 * G
  d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(G), each = 20)))
  d$y <- 1 + 0.5 * d$x + rnorm(n, 0, 1)
  d$kk <- 1e-4
  r <- cap(frm(bf(y ~ a + kk * u, a ~ 1 + x, u ~ 0 + (1 | g), nl = TRUE),
               data = d))
  f <- r$v
  se <- se_outer(f)
  H <- ns$fit_outer_hessian(f)$H
  on <- ns$outer_par_names(f)
  j <- match("theta_1", on)
  p <- f$opt$par
  f0 <- f$obj$fn(p)
  prof <- vapply(c(-12, -6, -3, 3, 6, 12), function(h) {
    q <- p; q[j] <- h; f0 - f$obj$fn(q)
  }, 0)
  f$obj$fn(p)
  cat(sprintf(paste0("seed %d G %2d theta_1 %.3g | SE %-10s | H[j,j] %.3g |",
                     " lost %s | warn %d msg %d | gain at theta -12 -6 -3",
                     " 3 6 12: %s\n"),
              seed, G, p[j], format(se[j], digits = 4), H[j, j], lost_of(f),
              length(r$w), length(r$m),
              paste(format(prof, digits = 3), collapse = " ")))
}
