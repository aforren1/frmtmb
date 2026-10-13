# Lane optrap, punch round 1 (m3): the review's mixture scale case on
# rellib-r7, the punch build (p1, mixtures without a unit) and the
# measured attempt with a unit on each component's mu (p1mix). Short is
# more than 1e-4 below the best any arm reached for the seed.
X <- do.call(rbind, lapply(c("r7", "p1", "p1mix"), function(b) {
  x <- read.delim(sprintf("dev/optrap-rev-log/scale-mixture-%s.tsv", b),
                  quote = "")
  x$build <- b
  x
}))
ref <- tapply(X$llu, X$seed, max, na.rm = TRUE)
X$gap <- X$llu - ref[as.character(X$seed)]
for (sc in unique(X$scale)) {
  cat(sprintf("%-6g", sc))
  for (b in c("r7", "p1", "p1mix")) {
    x <- X[X$scale == sc & X$build == b, ]
    sh <- x$gap < -1e-4
    cat(sprintf(" | %-5s short %2d (%8.3g) code!=0 %2d warned %2d", b,
                sum(sh), min(x$gap), sum(x$code != 0), sum(x$nwarn > 0)))
  }
  cat("\n")
}
