# The paired lane-against-base sweep for the residual hazard inside the
# declaring dpar: 12 seeds at each of four covariate scales.
d <- "C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-log/"
a <- read.csv(paste0(d, "alphascale-base.csv"))
b <- read.csv(paste0(d, "alphascale-fixed.csv"))
stopifnot(nrow(a) == nrow(b), identical(a$scale, b$scale),
          identical(a$seed, b$seed))
d_ll <- b$ll - a$ll
cat("paired fits:", nrow(a), "\n")
cat("lane BETTER than base:", sum(d_ll > 1e-6),
    "  worse:", sum(d_ll < -1e-6),
    "  tied:", sum(abs(d_ll) <= 1e-6), "\n")
cat("worst deficit:", format(min(d_ll), digits = 6),
    "  best gain:", format(max(d_ll), digits = 6), "\n")
cat("\nby scale (mean d logLik, lane minus base):\n")
for (s in unique(a$scale)) {
  k <- a$scale == s
  cat(sprintf("  %-8g n=%2d  mean %+10.5f  better %2d  worse %2d  worst %+9.5f\n",
              s, sum(k), mean(d_ll[k]), sum(d_ll[k] > 1e-6),
              sum(d_ll[k] < -1e-6), min(d_ll[k])))
}
cat("\nnonzero convergence codes: base", sum(a$conv != 0),
    " lane", sum(b$conv != 0), "\n")
