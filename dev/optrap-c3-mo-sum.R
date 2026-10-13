# Reviewer of lane optrap, item 5: summary of dev/optrap-mo-study.R rerun
# on the lane's final build and rellib-r7 (dev/optrap-rev-log/mo/).
# Lane optrap, punch round 1: arms p2 (punch round 2) and c3 (the consolidation fixes) added.
for (rng in list(1:200, 201:400)) for (arm in c("r7", "optrap", "p2", "c3")) {
  f <- Sys.glob(sprintf(c("dev/optrap-rev-log/mo/mo-%s-*.tsv",
                         "dev/optrap-rev-log/mo/%s-*.tsv",
                         "dev/optrap-log/p2/mo/%s-*.tsv",
                         "dev/optrap-log/c3/mo/%s-*.tsv"), arm))
  X <- do.call(rbind, lapply(f, read.delim))
  X <- X[X$seed %in% rng, ]
  X <- X[order(X$seed), ]
  cat(sprintf(paste0("%-6s seeds %d-%d n %d | at the exact max (gap < 1e-6)",
                     " %d | max gap %.4g | gap > 1e-6: %s | code != 0 %d ",
                     "| warned %d | evals %d\n"),
              arm, min(rng), max(rng), nrow(X), sum(X$gap < 1e-6),
              max(X$gap), paste(X$seed[X$gap >= 1e-6], collapse = ","),
              sum(X$code != 0), sum(X$nwarn > 0), sum(X$nfev)))
}
