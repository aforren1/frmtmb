# Lane optrap, punch round 1: dev/optrap-rev-scale-sum.R with the lane's
# punch-round build (tag p1) beside rellib-r7 and round 0 (tag lane). The
# reference of a seed is the best unscaled log likelihood any arm (both
# builds, every scale, the external fit) reached; "short" is more than
# 1e-4 below it.
#   OPTRAP_P1TAG=p1g Rscript dev/optrap-p1-scale-sum.R (default p1)
fs <- Sys.glob("dev/optrap-rev-log/scale-*-*.tsv")
X <- do.call(rbind, lapply(fs, function(f) {
  x <- read.delim(f, stringsAsFactors = FALSE, quote = "")
  x$build <- sub(".*-([a-z0-9]+)[.]tsv$", "\\1", basename(f))
  x
}))
X <- X[X$build %in% c("lane", "r7", Sys.getenv("OPTRAP_P1TAG", "p1")), ]
key <- paste(X$case, X$seed)
ref <- tapply(pmax(X$llu, X$ext, na.rm = TRUE), key, max, na.rm = TRUE)
X$gap <- X$llu - ref[key]
X$extgap <- X$ext - ref[key]
cat(sprintf("%-12s %-6s | %-24s | %-24s | ext short\n", "case", "scale",
            "r7: short(min) code!=0 warned err",
            "lane: short(min) code!=0 warned err | p1: same"))
for (cs in unique(X$case)) for (sc in unique(X$scale[X$case == cs])) {
  out <- character()
  for (b in c("r7", "lane", Sys.getenv("OPTRAP_P1TAG", "p1"))) {
    x <- X[X$case == cs & X$scale == sc & X$build == b, ]
    sh <- x$gap < -1e-4 & !is.na(x$gap)
    out <- c(out, sprintf("%2d/%2d (%9.3g) %2d %2d %2d", sum(sh), nrow(x),
                          if (any(sh)) min(x$gap[sh]) else
                            max(abs(x$gap), na.rm = TRUE),
                          sum(x$code != 0, na.rm = TRUE),
                          sum(x$nwarn > 0, na.rm = TRUE),
                          sum(is.na(x$code))))
  }
  x <- X[X$case == cs & X$scale == sc & X$build == "r7", ]
  eg <- x$extgap
  cat(sprintf("%-12s %-6g | %s | %s | %s | %s\n", cs, sc, out[1], out[2],
              out[3],
              if (all(is.na(eg))) "-" else
                sprintf("%d (%.3g)", sum(eg < -1e-4, na.rm = TRUE),
                        min(eg, na.rm = TRUE))))
}
cat("\nunits seen on the lane:\n")
L <- X[X$build == "lane", ]
print(table(L$case, L$units != "1"))
cat("\nwarnings:\n")
W <- X[nzchar(X$warn) & !is.na(X$warn), ]
print(table(paste(W$case, W$build, substr(W$warn, 1, 60))))
