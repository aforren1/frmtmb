# Reviewer of lane illse: raw polynomials with (1 | g) on data with no
# group effect (the sd runs to zero), seeds 1-20, degrees 5-8. Per fit:
# the sd, whether the route took the fit, what theta_1 got (lost reason
# or its SE), the warnings.
#   REV_ARM=lane|base Rscript dev/illse-rev-edge4.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
for (deg in 5:8) for (seed in 1:20) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = stats::runif(n, 1, 2), g = factor(rep(1:20, length.out = n)))
  d$y <- sin(3 * d$x) + stats::rnorm(n, 0, 0.3)
  xn <- paste0("x", seq_len(deg))
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "), "+ (1 | g)"))
  r <- cap(frm(fo, data = d))
  f <- r$v
  if (inherits(f, "rev_err")) { cat(deg, seed, "ERROR\n"); next }
  sd <- exp(f$opt$par[names(f$opt$par) == "theta"])
  l <- ns$sdr_of(f)$se_lost
  V <- suppressWarnings(vcov(f, full = TRUE))
  se_t <- sqrt(V[nrow(V), nrow(V)])
  cat(sprintf("deg %d s%2d | code %d | sd %.1e | route %s | theta_1 %s | n lost %d | warn %s\n",
              deg, seed, f$opt$convergence, sd, !is.null(ns$sdr_of(f)$cov_fixed_root),
              if ("theta_1" %in% names(l)) paste0("lost:", l[["theta_1"]]) else sprintf("SE %.2e", se_t),
              length(l), wtag(r$w)))
}
