# Lane famend, fix round F1: where the GEV's tape Hessian is not finite
# on outside rows after the ratio clamp, NaN or Inf, and the true size
# of the second derivative there (exp(-h) C^2 x^2).
.libPaths(c("C:/Users/adf44/source/r/wt-famend-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
fam <- gen_extreme_value()
for (xi in c(0.06, 1)) {
  xs <- seq(-1 / xi - 0.01, -1 / xi - 60, by = -0.01)
  for (x in xs) {
    y <- 0.3 + exp(-0.2) * x
    p <- c(0.3, -0.2, log1p(xi))
    tp <- RTMB::MakeTape(function(p) {
      fam$lpdf(y, list(mu = p[1], sigma = exp(p[2]), xi = expm1(p[3]),
                       .eta_sigma = p[2], .eta_xi = p[3]), list())
    }, p)
    H <- tp$jacfun()$jacobian(p)
    if (any(!is.finite(H))) {
      cat(sprintf("xi %g x %.2f h %.1f: NaN %d Inf %d | finite max %.3g\n",
                  xi, x, x * 34.539, sum(is.nan(H)), sum(is.infinite(H)),
                  max(abs(H[is.finite(H)]))))
    }
  }
}
