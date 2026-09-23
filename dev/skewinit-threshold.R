# Where the |alpha| distribution separates: the stationary point the
# optimizer stops at, against the smallest genuine optimum seen.
d <- read.csv("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-log/calib.csv")
d$best <- pmax(d$ll_sn, d$ll_refit)
d$short <- d$best - d$ll
stat <- d[d$short > 1e-6, ]
good <- d[d$short <= 1e-6, ]
cat("rows short > 1e-6 (not at the optimum):", nrow(stat), "\n")
cat("  max |alpha|:", format(max(abs(stat$alpha)), digits = 8), "\n")
cat("  sorted top 8:",
    paste(format(sort(abs(stat$alpha), decreasing = TRUE)[1:8], digits = 4),
          collapse = " "), "\n")
cat("rows short <= 1e-6 (at the optimum):", nrow(good), "\n")
cat("  min |alpha|:", format(min(abs(good$alpha)), digits = 8), "\n")
cat("  sorted bottom 8:",
    paste(format(sort(abs(good$alpha))[1:8], digits = 4), collapse = " "),
    "\n")
for (t in c(0.01, 0.02, 0.05, 0.1, 0.2)) {
  fire <- abs(d$alpha) < t
  cat(sprintf(
    "t=%-5g fires %3d/240  FN(miss short>1e-6) %2d  FP(fires, short<=1e-6) %2d",
    t, sum(fire), sum(!fire & d$short > 1e-6), sum(fire & d$short <= 1e-6)))
  cat(sprintf("  FP(fires, gain<=0.01) %2d\n",
              sum(fire & d$short <= 0.01)))
}
