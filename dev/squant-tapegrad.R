# Lane squant, fix round 1: the worst relative error of the taped
# d lF / d mu and d lS / d mu of quantile_elf()'s cdf_log against -f / F
# and f / S, per (p, smoothing, side), on the grid of the test "the taped
# log CDF pair has the derivatives -f / F and f / S".
# Run: Rscript dev/squant-tapegrad.R <lib> > dev/squant-out/tapegrad-<arm>.txt
.libPaths(c(commandArgs(TRUE)[1], "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "BLAS", extSoftVersion()[["BLAS"]], "\n")
t <- c(-700, -690, -689, -500, -300, -200, -150, -100, -40, -10, -1, 0,
       1, 10, 40, 100, 150, 200, 300, 500, 689, 690, 700)
for (p in c(0.1, 0.5, 0.9)) for (lam in c(0.25, 1)) {
  fam <- quantile_elf(smoothing = lam)
  q <- lam * t
  dp <- function(m) list(mu = m + 0 * q, sigma = 1, quantile = p)
  for (side in c("lF", "lS")) {
    tp <- RTMB::MakeTape(function(m) fam$cdf_log(q, dp(m), list())[[side]], 0)
    g <- as.numeric(tp$jacobian(0))
    v <- fam$cdf_log(q, dp(0), list())[[side]]
    an <- (if (side == "lF") -1 else 1) * exp(fam$lpdf(q, dp(0), list()) - v)
    err <- abs(g / an - 1)
    cat(sprintf("p %.1f lam %.2f %s: worst %.2g at t = %g | %s\n", p, lam,
                side, max(err), t[which.max(err)],
                paste(sprintf("%g:%.1g", t, err), collapse = " ")))
  }
}
