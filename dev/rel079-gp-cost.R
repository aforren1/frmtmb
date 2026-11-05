# The cost line of the exact gp() on the merged build, one fit at a time, after
# release fix R1 of the ncgp review (six length scales in the warm start). The
# design is the lane's (NEWS, Costs): y ~ gp(x) on 60, 200 and 500 distinct positions,
# and on 1000 rows at 500 positions. The start's share is read by timing
# gp_marginal_max() on the frame.
#   Rscript dev/rel079-gp-cost.R <library>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), as.character(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
one <- function(n, m, seed) {
  set.seed(seed)
  x <- if (n == m) sort(runif(m, 0, 10)) else sample(sort(runif(m, 0, 10)), n, TRUE)
  d <- data.frame(x = x, y = sin(x) + rnorm(n, 0, 0.3))
  t0 <- proc.time()[[3]]
  fit <- suppressWarnings(suppressMessages(frm(y ~ gp(x), data = d)))
  tf <- proc.time()[[3]] - t0
  ts <- NA_real_
  if (exists("gp_marginal_max", ns)) {
    t0 <- proc.time()[[3]]
    try(ns$gp_marginal_max(fit$frame), silent = TRUE)
    ts <- proc.time()[[3]] - t0
  }
  cat(sprintf("n %5d positions %4d seed %d fit %.2f s, start %.2f s, conv %d, logLik %.4f\n",
              n, m, seed, tf, ts, fit$opt$convergence, as.numeric(logLik(fit))))
}
for (s in 1:2) {
  one(60, 60, s)
  one(200, 200, s)
}
one(500, 500, 1)
one(1000, 500, 1)
