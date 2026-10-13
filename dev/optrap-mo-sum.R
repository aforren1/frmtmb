# Lane optrap: summary of dev/optrap-mo-study.R runs, per arm.
#   Rscript dev/optrap-mo-sum.R
for (arm in c("r7", "optrap")) {
  f <- Sys.glob(sprintf("dev/optrap-log/mo-%s-*.tsv", arm))
  X <- do.call(rbind, lapply(f, read.delim))
  X <- X[order(X$seed), ]
  cat(sprintf(paste0("%-6s seeds %d | at the exact max (gap < 1e-6) %d | ",
                     "max gap %.4g (seed %s) | gap > 1e-6: %s | code != 0 %d ",
                     "| warned %d | evals %d\n"),
              arm, nrow(X), sum(X$gap < 1e-6), max(X$gap),
              paste(X$seed[X$gap == max(X$gap)], collapse = ","),
              paste(X$seed[X$gap >= 1e-6], collapse = ","),
              sum(X$code != 0), sum(X$nwarn > 0), sum(X$nfev)))
}
