# Lane optrap, punch round 1 (m6): huber() fits short of the best any arm
# reached, split into flagged (code != 0 or a warning) and silent, per
# build and scale; and whether a silent short of p1 is short on r7 too.
X <- do.call(rbind, lapply(c("r7", "lane", "p1i"), function(b) {
  x <- read.delim(sprintf("dev/optrap-rev-log/scale-huber-%s.tsv", b),
                  quote = "")
  x$build <- b
  x
}))
ref <- tapply(X$llu, X$seed, max, na.rm = TRUE)
X$gap <- X$llu - ref[as.character(X$seed)]
X$short <- X$gap < -1e-4
X$flag <- X$code != 0 | X$nwarn > 0
cat("short / silent short (max silent gap) per build and scale\n")
for (sc in unique(X$scale)) {
  cat(sprintf("%-6g", sc))
  for (b in c("r7", "lane", "p1i")) {
    x <- X[X$scale == sc & X$build == b, ]
    s <- x$short & !x$flag
    cat(sprintf(" | %s %2d / %2d (%8.3g)", b, sum(x$short), sum(s),
                if (any(s)) min(x$gap[s]) else 0))
  }
  cat("\n")
}
p1 <- X[X$build == "p1i" & X$short & !X$flag, c("seed", "scale", "gap")]
r7 <- X[X$build == "r7" & X$scale == 1, c("seed", "gap")]
cat("\np1i silent shorts, with r7's gap at scale 1 for the same seed\n")
p1$r7_scale1 <- r7$gap[match(p1$seed, r7$seed)]
print(p1[order(p1$seed), ], digits = 3, row.names = FALSE)
