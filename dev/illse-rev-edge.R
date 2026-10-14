# Reviewer of lane illse: the route's gate on an sd at its edge. A raw
# polynomial with (1 | g) whose group sd estimates at zero (the data have
# no group effect): what is lost, lane and base.
#   REV_ARM=lane|base Rscript dev/illse-rev-edge.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
for (deg in c(5, 7)) for (seed in 1:3) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = stats::runif(n, 1, 2), g = factor(rep(1:20, length.out = n)))
  d$y <- sin(3 * d$x) + stats::rnorm(n, 0, 0.3)
  for (k in seq_len(deg)) d[[paste0("x", k)]] <- d$x^k
  fo <- stats::as.formula(paste("y ~", paste(paste0("x", seq_len(deg)), collapse = " + "), "+ (1 | g)"))
  r <- cap(frm(fo, data = d))
  f <- r$v
  l <- ns$sdr_of(f)$se_lost
  cat(sprintf("deg %d s%d | code %d | sd %.2e | route %s | lost %s | warn %s\n", deg, seed,
              f$opt$convergence, exp(f$opt$par[names(f$opt$par) == "theta"]),
              !is.null(ns$sdr_of(f)$cov_fixed_root),
              if (length(l)) paste(paste0(names(l), ":", l), collapse = " ") else "-",
              if (length(r$w)) substr(paste(r$w, collapse = " / "), 1, 120) else "-"))
}
