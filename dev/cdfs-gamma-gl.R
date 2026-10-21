# Prototype: the incomplete gamma's tail on the far side of x from the
# mode a - 1 by Gauss-Legendre quadrature of the integrand divided by
# its peak (Numerical Recipes 3rd ed., gammpapprox(), with its ranges),
# for the middle band x in [0.7 (a + 1), a + 1 + 3 sqrt(a)]. How many
# nodes, and from which shape? Output: dev/cdfs-out/gamma-gl.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
E <- new.env(parent = asNamespace("frmtmb"))
sys.source("R/families-cdf.R", envir = E)
gl_unit <- get("gl_unit", E)
lgamma_binet <- get("lgamma_binet", E)
tail_gl <- function(x, a, g) {
  a1 <- a - 1
  s <- sqrt(a1)
  # log of a1^a1 e^-a1 / Gamma(a), from Stirling with Binet's remainder
  lpre <- a1 * log1p(-1 / a) + 1 - 0.5 * log(2 * pi * a) - lgamma_binet(a)
  up <- x >= a1
  lim <- if (up) max(a1 + 11.5 * s, x + 6 * s) else
    max(0, min(a1 - 7.5 * s, x - 5 * s))
  t <- x + (lim - x) * g$x
  d <- (t - a1) / a1
  f <- exp(a1 * (log1p(d) - d))
  ans <- abs(lim - x) * sum(g$w * f)
  c(up = up, l = lpre + log(ans))
}
rel <- function(v, r) abs(v - r) / pmax(1, abs(r))
for (n in c(18, 24, 32, 40)) {
  g <- gl_unit(n)
  out <- character()
  for (a in c(30, 60, 100, 150, 400, 1e3, 1e4, 1e5, 1e6)) {
    xs <- c(0.7 * (a + 1), 0.8 * (a + 1), a - 3 * sqrt(a), a - 1.2,
            a - 0.7, a, a + sqrt(a), a + 1 + 2 * sqrt(a),
            a + 1 + 3 * sqrt(a))
    worst <- 0
    for (x in xs) {
      r <- tail_gl(x, a, g)
      ref <- if (r[["up"]] == 1) {
        pgamma(x, a, lower.tail = FALSE, log.p = TRUE)
      } else {
        pgamma(x, a, log.p = TRUE)
      }
      worst <- max(worst, rel(r[["l"]], ref))
    }
    out <- c(out, sprintf("%g:%.0e", a, worst))
  }
  cat(sprintf("n %2d  %s\n", n, paste(out, collapse = " ")))
}
