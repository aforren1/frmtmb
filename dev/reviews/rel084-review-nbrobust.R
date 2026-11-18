# rel084 release review, nbrobust fix round 1 (records only): the record
# says a negbinomial fit started far out on the poisson plateau of
# overdispersed data stops at its start on both builds; 0.83.0 warned
# "Large maximum absolute gradient" by accident, 0.84.0 says only the SE
# check's false "flat". New data: seeds 51 and 52, n = 200, size 2 and 1.
#   REVLIB=r23|r22 Rscript dev/reviews/rel084-review-nbrobust.R
source("dev/reviews/rel084-review-common.R")
for (sd in c(51, 52)) {
  set.seed(sd)
  d <- data.frame(x = rnorm(200))
  d$y <- rnbinom(200, size = if (sd == 51) 2 else 1, mu = exp(1 + 0.4 * d$x))
  ref <- conds_of(frm(y ~ x, family = negbinomial(), data = d))
  llr <- as.numeric(logLik(ref$v))
  cat(sprintf("seed %d default: code %d lsh %.4f ll %.6f [%s]\n", sd,
              ref$v$opt$convergence, ref$v$opt$par[["betad"]], llr, cls_line(ref)))
  for (st in c(20, 25, 30, 35)) {
    r <- conds_of(frm(y ~ x, family = negbinomial(), data = d,
                      start = list(betad = st)))
    f <- r$v
    if (is.null(f)) { cat("  start", st, "ERROR", cls_line(r), "\n"); next }
    wm <- vapply(Filter(function(s) s$k == "W", r$said), `[[`, "", "msg")
    cat(sprintf("  start %d code %d lsh %.4f ll %.6f below %.3f grad %s flat %s [%s]\n",
                st, f$opt$convergence, f$opt$par[["betad"]],
                as.numeric(logLik(f)), llr - as.numeric(logLik(f)),
                any(grepl("Large maximum absolute gradient", wm)),
                any(grepl("do not determine", wm) | grepl("flat", wm)),
                cls_line(r)))
  }
}
cat("DONE\n")
