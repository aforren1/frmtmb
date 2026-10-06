# Reviewer, punch 2b: per fit, dev/nanse-rev4-ridge.R's lane against the
# round-2 merge build (r2) and base: finite ridge SEs, lost counts.
#   Rscript dev/nanse-rev4-cmp.R
rd <- function(a) {
  x <- readLines(sprintf("C:/Users/adf44/source/r/frmtmb-wt-nanse/dev/nanse-rev4-log/ridge-%s.txt", a))
  i <- grep("^ design", x)[1]
  j <- grep("^fits with a finite", x) - 2
  read.table(text = x[i:j], header = TRUE)
}
L <- rd("lane"); R <- rd("r2"); B <- rd("base")
k <- c("design", "sdg", "seed")
m <- merge(merge(L, R, by = k, suffixes = c(".lane", ".r2")), B[, c(k, "ridge_finite", "lost")], by = k)
cat("fits:", nrow(m), "\n")
cat("ridge_finite differs lane vs r2:", sum(m$ridge_finite.lane != m$ridge_finite.r2), "\n")
cat("ridge_finite differs lane vs base:", sum(m$ridge_finite.lane != m$ridge_finite), "\n")
cat("lost differs lane vs r2:", sum(m$lost.lane != m$lost.r2), "\n")
print(m[m$lost.lane != m$lost.r2 | m$ridge_finite.lane != m$ridge_finite.r2,
        c(k, "ridge_finite.lane", "ridge_finite.r2", "ridge_finite", "lost.lane", "lost.r2", "min_finite_se.lane")], row.names = FALSE)
cat("\nfinite ridge SEs on the lane, smallest per design:\n")
print(aggregate(min_finite_se.lane ~ design + sdg, m, min))
