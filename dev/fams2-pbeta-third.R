# Which third-order partials of log(RTMB::pbeta(q, a, b)) are non-finite
# at the ordinary points dev/fams2-pbeta-grid.R flagged? Arguments in
# (q, log a, log b), taped directly.
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(RTMB)
f <- function(p) log(pbeta(p[1], exp(p[2]), exp(p[3])))
x0 <- c(0.2, 0, 1)
F <- MakeTape(f, x0)
J <- F$jacfun(); H <- J$jacfun(); T3 <- H$jacfun()
pts <- rbind(c(0.268941421, 1, 4.481689), c(0.091122961, 0.04978707, 1.648721),
             c(0.268941421, 0.22313016, 2.718282), c(0.1, 3, 20),
             c(0.1, 0.5, 0.5), c(0.3, 2, 2))
for (i in seq_len(nrow(pts))) {
  x <- c(pts[i, 1], log(pts[i, 2]), log(pts[i, 3]))
  t3 <- T3(x)   # 9 x 3: rows are (i, j) of the Hessian, columns the third index
  cat(sprintf("q=%.4f a=%.4f b=%.4f value %.6f\n", pts[i, 1], pts[i, 2],
              pts[i, 3], F(x)))
  lab <- c("q", "a", "b")
  nf <- which(!is.finite(t3), arr.ind = TRUE)
  if (!nrow(nf)) { cat("  all finite\n"); next }
  for (r in seq_len(nrow(nf))) {
    ij <- nf[r, 1]; k <- nf[r, 2]
    ii <- (ij - 1) %% 3 + 1; jj <- (ij - 1) %/% 3 + 1
    cat("  d3/d", lab[ii], "d", lab[jj], "d", lab[k], " = ", t3[ij, k], "\n",
        sep = "")
  }
}
