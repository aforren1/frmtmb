# Reviewer of lane optrap: dev/optrap-rev-mo-resid.R's 400 seeds.
X <- do.call(rbind, lapply(Sys.glob("dev/optrap-rev-log/mo-resid/*.tsv"),
                           read.delim))
for (rng in list(1:200, 201:400)) {
  x <- X[X$seed %in% rng, ]
  cat(sprintf("resid-unit seeds %d-%d n %d | at the exact max %d | max gap %.4g | gap > 1e-6: %s | code != 0 %d | warned %d | evals %d\n",
              min(rng), max(rng), nrow(x), sum(x$gap < 1e-6), max(x$gap),
              paste(x$seed[x$gap >= 1e-6], collapse = ","), sum(x$code != 0),
              sum(x$nwarn > 0), sum(x$nfev)))
}
