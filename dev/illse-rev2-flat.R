# Reviewer of lane illse, round 2, item 5: can headroom_flat_rows() hide
# a real non-convergence? (a) Synthetic rows: a row empty by rounding
# with a large gradient, and one at the gradient bound. (b) A real fit
# whose group sd is stuck near zero from the start while the data have a
# group effect (sd 0.5): what the convergence check and the SE check say,
# against the fit from the default start.
#   REV_ARM=lane|base Rscript dev/illse-rev2-flat.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
if (exists("headroom_flat_rows", ns)) {
  H <- diag(c(100, 50, 1e-12))
  for (g3 in c(1e-2, 1e-3, 1e-4, 4e-5)) {
    fl <- ns$headroom_flat_rows(H, c(0, 0, g3), rep(1, 3), 1e-3)
    cat(sprintf("synthetic: row max 1e-12, gradient %.0e -> left out %s; a Newton step at that curvature gains %.1e\n",
                g3, fl[3], g3^2 / (2 * 1e-12)))
  }
  H2 <- H; H2[3, 3] <- NaN
  cat("synthetic: non-finite row, gradient 1e-2 -> left out",
      ns$headroom_flat_rows(H2, c(0, 0, 1e-2), rep(1, 3), 1e-3)[3], "\n")
}
for (seed in 1:4) for (st in c(-6, -10, -15)) for (im in c(1000, 3)) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = stats::rnorm(n), g = factor(rep(1:20, length.out = n)))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(20, 0, 0.5)[d$g] + stats::rnorm(n, 0, 0.4)
  ref <- suppressWarnings(frm(y ~ x + (1 | g), data = d))
  r <- cap(frm(y ~ x + (1 | g), data = d, start = list(theta = st),
               control = frmtmb_control(optCtrl = list(iter.max = im, eval.max = 1000))))
  f <- r$v
  if (inherits(f, "rev_err")) { cat("ERROR", f, "\n"); next }
  l <- ns$sdr_of(f)$se_lost
  w <- paste(r$w, collapse = " / ")
  clause <- if (grepl("so the estimates are at the maximum", w, fixed = TRUE)) "AT_MAX" else
    if (grepl("short of the maximum", w, fixed = TRUE)) "SHORT" else
      if (grepl("could not be measured", w, fixed = TRUE)) "UNUSABLE" else "-"
  cat(sprintf("s%d start log sd %3d iter.max %4d | code %d | sd %.2e (ref %.2f) | short of ref %.3e | conv clause %-8s | theta_1 %s | warn %s\n",
              seed, st, im, f$opt$convergence, exp(f$opt$par[["theta"]]),
              exp(ref$opt$par[["theta"]]), f$opt$objective - ref$opt$objective, clause,
              if ("theta_1" %in% names(l)) l[["theta_1"]] else "kept", wtag(r$w)))
}
