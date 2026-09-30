# Punch 2: is the worst sweep point (a = 0.528, b = 8.93e6, 1.8 sd above
# m) short of steps? log_ibeta_half() at N = 50, 100, 200, 400.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
f <- frmtmb:::log_ibeta_half
for (p in list(c(0.528, 8.93e6, 1.8), c(3.23, 9.34e6, 1.38), c(21.8, 6.2e6, 0.812))) {
  a <- p[1]; b <- p[2]; s <- a + b
  m <- (a + 1) / (s + 2); sd <- sqrt(a * b / (s * s * (s + 1)))
  x <- m + p[3] * sd
  ref <- pbeta(x, a, b, log.p = TRUE)
  cat(sprintf("a %g b %g k %g:", a, b, p[3]))
  for (N in c(50L, 100L, 200L, 400L)) {
    cat(sprintf("  N %d err %.2e", N, abs(f(x, a, b, N) - ref) / max(1, abs(ref))))
  }
  cat("\n")
}
