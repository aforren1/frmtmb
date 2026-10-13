# Lane optrap, punch round 1 (m6): huber() at scale 1 per seed, rellib-r7
# against round 0 (lane) and the punch build (p1), from
# dev/optrap-rev-scale.R's logs.
X <- lapply(c("r7", "lane", "p1"), function(b) {
  x <- read.delim(sprintf("dev/optrap-rev-log/scale-huber-%s.tsv", b),
                  quote = "")
  x$build <- b
  x
})
X <- do.call(rbind, X)
ref <- tapply(pmax(X$llu, X$ext, na.rm = TRUE), X$seed, max, na.rm = TRUE)
X$gap <- X$llu - ref[as.character(X$seed)]
for (sc in unique(X$scale)) {
  cat("== scale", sc, "\n")
  w <- reshape(X[X$scale == sc, c("seed", "build", "gap", "code", "units")],
               idvar = "seed", timevar = "build", direction = "wide")
  print(w[order(w$seed), ], digits = 3, row.names = FALSE)
}
