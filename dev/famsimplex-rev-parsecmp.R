a <- readRDS("famsimplex-rev-parse-lane.rds"); b <- readRDS("famsimplex-rev-parse-base.rds")
for (k in names(b)) if (is.null(b[[k]]$err)) {
  cat(k, " logLik identical:", identical(a[[k]]$ll, b[[k]]$ll),
      " coef identical:", identical(a[[k]]$coef, b[[k]]$coef),
      " se identical:", identical(a[[k]]$se, b[[k]]$se),
      " max rel coef diff:", signif(max(abs(a[[k]]$coef - b[[k]]$coef) / abs(b[[k]]$se)), 3), "\n")
}
