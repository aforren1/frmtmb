# Review B: zero_inflated_poisson on data with NO zero inflation, over
# seeds. About half the seeds put the zi MLE inside its range (a sampling
# fluke of a few zeros), with a small likelihood gain over zi = 0 and a
# large standard error. Which are held, and does predict() agree with
# fitted() on the rest?
# Usage: Rscript ppdraw-review-b.R <lib or "base"> [n] [seeds]
args <- commandArgs(TRUE)
lib <- if (length(args) && args[1] != "base") args[1] else NULL
n <- if (length(args) >= 2) as.integer(args[2]) else 200L
nseed <- if (length(args) >= 3) as.integer(args[3]) else 30L
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
source("C:/Users/adf44/source/r/frmtmb-wt-ppdraw/dev/reviews/ppdraw-review-common.R")
cat("n", n, "seeds", nseed, "\n")
cat(sprintf("%4s %9s %9s %10s %10s %7s %8s %8s %7s %5s %s\n", "seed",
            "zi_est", "zi_sd", "minus2sd", "plus2sd", "LRgain", "pmean",
            "fmean", "z", "held", "share_zi_gt_0.5"))
for (seed in seq_len(nseed)) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n))
  d$y <- stats::rpois(n, exp(1 + 0.3 * d$x))
  fit <- q(frm(bf(y ~ x), family = zero_inflated_poisson(), data = d))
  f0 <- q(frm(bf(y ~ x), family = poisson(), data = d))
  gain <- -fit$opt$objective - (-f0$opt$objective)
  fs <- flat_sides(fit)
  j <- grep("^zi_", rownames(fs))
  r <- pz(fit, ndraws = 400L, seed = 1)
  # share of the draws' zi above 0.5, from the same draw law
  set.seed(1)
  dr <- q(frmtmb:::predict_par_drawer(fit, TRUE, 400L))
  mp <- frmtmb:::fit_draw_space(fit)$map
  jz <- grep("^zi_", mp$names)
  v0 <- frmtmb:::fit_outer_vector(fit, mp)
  zis <- vapply(seq_len(400), function(s) {
    v0[jz] + attr(dr, "delta")(s)[jz]
  }, 0)
  cat(sprintf("%4d %9.3f %9.3g %10.3g %10.3g %7.2g %8.4f %8.4f %7.2f %5s %.3f\n",
              seed, fs[j, "est"], fs[j, "sd"], fs[j, "minus2sd"],
              fs[j, "plus2sd"], gain, r$pmean, r$fmean, r$z,
              length(r$held) > 0, mean(stats::plogis(zis) > 0.5)))
}
