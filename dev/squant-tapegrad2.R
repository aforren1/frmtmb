# Lane squant, fix round 1: at the points where dev/squant-tapegrad.R's
# relative error is not small after the fix, the size of the derivative
# itself against the tail's slope max(a, b) / (lambda sigma), and the
# tape's value. Run: Rscript dev/squant-tapegrad2.R > dev/squant-out/tapegrad2.txt
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
t <- c(-700, -690, -689, -500, -300, -200, -150, -100, -40, -10, -1, 0,
       1, 10, 40, 100, 150, 200, 300, 500, 689, 690, 700)
for (p in c(0.1, 0.5, 0.9)) for (lam in c(0.25, 1)) {
  fam <- quantile_elf(smoothing = lam)
  q <- lam * t
  dp <- function(m) list(mu = m + 0 * q, sigma = 1, quantile = p)
  slope <- lam * max(p, 1 - p) / lam
  for (side in c("lF", "lS")) {
    tp <- RTMB::MakeTape(function(m) fam$cdf_log(q, dp(m), list())[[side]], 0)
    g <- as.numeric(tp$jacobian(0))
    v <- fam$cdf_log(q, dp(0), list())[[side]]
    an <- (if (side == "lF") -1 else 1) * exp(fam$lpdf(q, dp(0), list()) - v)
    bad <- abs(g / an - 1) > 1e-8
    if (any(bad)) {
      cat(sprintf("p %.1f lam %.2f %s: t %s | |an| / slope %s | tape %s\n",
                  p, lam, side, paste(t[bad], collapse = " "),
                  paste(sprintf("%.1e", abs(an[bad]) / slope), collapse = " "),
                  paste(sprintf("%.1e", g[bad]), collapse = " ")))
    }
    cat(sprintf("p %.1f lam %.2f %s: max |tape - an| / slope %.2g\n", p, lam,
                side, max(abs(g - an)) / slope))
  }
}
