# Lane optima, item 2: where is RTMB's pnorm(log.p = TRUE) derivative
# wrong or not finite in the lower tail? A log grid of x in [1, 1e300]
# and a fine grid about 4.463108e9, d/dx log Phi(-x) against its
# asymptote -x (relative error), counting non-finite values.
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
chk <- function(x, lab) {
  tp <- MakeTape(function(e) pnorm(-e, log.p = TRUE), x)
  v <- tp(x)
  g <- diag(tp$jacobian(x))
  # exact: -phi(x)/Phi(-x), which tends to -x (1 + 1/x^2)^-1 ...
  ref <- -x / (1 - 1 / x^2 + 3 / x^4)
  ok <- x > 30
  cat(sprintf("%-26s n %d | value non-finite %d | derivative non-finite %d | max rel err (x > 30) %.3g at x = %.4g\n",
              lab, length(x), sum(!is.finite(v)), sum(!is.finite(g)),
              max(abs(g[ok] / ref[ok] - 1), na.rm = TRUE),
              x[ok][which.max(abs(g[ok] / ref[ok] - 1))]))
  invisible(list(x = x, g = g))
}
r1 <- chk(10^seq(0, 300, length.out = 3001), "log grid 1 .. 1e300")
r2 <- chk(4.463108491e9 + seq(-5, 5, length.out = 2001), "about 4.463108e9")
r3 <- chk(10^seq(5, 7, length.out = 2001), "log grid 1e5 .. 1e7")
bad <- r2$x[!is.finite(r2$g)]
if (length(bad)) cat("non-finite at", format(head(bad), digits = 17), "\n")
