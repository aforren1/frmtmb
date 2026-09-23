# Lane wt-reunc: summarize dev/reunc-coverage.R's result files into the
# block dev/reunc-findings.md carries verbatim.
#
#   Rscript dev/reunc-coverage-sum.R dev/reunc-log/cov
a <- commandArgs(trailingOnly = TRUE)
dir <- if (length(a)) a[1] else "dev/reunc-log/cov"
files <- sort(list.files(dir, pattern = "^cov-.*[.]rds$", full.names = TRUE))
for (f in files) {
  x <- readRDS(f)
  r <- x$res
  ok <- r[!r$failed, ]
  reps <- sort(unique(ok$rep))
  cat(sprintf("design %s  build %s (frmtmb %s)  script dev/reunc-coverage.R\n",
              x$design, if (x$lane) "lane" else "base", x$version))
  cat(sprintf(paste0("replicates with results %d of %d (distinct seeds %d, ",
                     "failed %d)  points per replicate %d  ndraws %d  ",
                     "seed0 %d\n"),
              length(reps), x$nrep, length(unique(ok$seed)),
              length(unique(r$rep[r$failed])), x$m, x$ndraws, x$seed0))
  missing <- setdiff(seq_len(x$nrep), c(reps, r$rep[r$failed]))
  if (length(missing)) cat("GAPS in the replicate grid:", missing, "\n")
  cat(sprintf("%-12s %11s %7s %17s %7s %7s %7s %8s\n", "arm", "hits",
              "cover", "binomial 95% CI", "rep se", "z", "p", "width"))
  for (arm in unique(ok$arm)) {
    s <- ok[ok$arm == arm, ]
    k <- sum(s$hit)
    n <- nrow(s)
    ci <- stats::binom.test(k, n, 0.95)$conf.int
    pr <- tapply(s$hit, s$rep, mean)
    se <- stats::sd(pr) / sqrt(length(pr))
    z <- (mean(pr) - 0.95) / se
    cat(sprintf("%-12s %5d/%5d %7.4f  (%.4f, %.4f) %7.4f %7.2f %7.4f %8.4f\n",
                arm, k, n, k / n, ci[1], ci[2], se, z,
                2 * stats::pnorm(-abs(z)), stats::median(s$width)))
  }
  # the conditional reading, which the interval does NOT claim: coverage
  # by the size of the group's realized effect, terciles of |u|
  q <- stats::quantile(ok$absu[ok$arm == "oracle_pred"], c(1 / 3, 2 / 3))
  cat("coverage by |u| tercile (the claim is the average, not these):\n")
  for (arm in intersect(c("pred_ml", "pred_reml", "oracle_pred", "fit_ml",
                          "oracle_mean"), unique(ok$arm))) {
    s <- ok[ok$arm == arm, ]
    tc <- cut(s$absu, c(-Inf, q, Inf), labels = c("low", "mid", "high"))
    v <- tapply(s$hit, tc, mean)
    cat(sprintf("  %-12s low %.4f  mid %.4f  high %.4f\n", arm, v[1], v[2],
                v[3]))
  }
  cat("\n")
}
