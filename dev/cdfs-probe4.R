# log_ibeta_half() restricted to x <= 0.6, by shape pair, and the two
# sides it would be blended from. Output: dev/cdfs-out/probe4.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
lih <- asNamespace("frmtmb")$log_ibeta_half
rel <- function(a, b) abs(a - b) / pmax(1, abs(b))
sh <- c(0.01, 0.05, 0.2, 0.5, 1, 2, 7.3, 40, 300, 5000, 1e5, 1e6)
x <- c(1e-300, 1e-100, 1e-30, 1e-12, 1e-6, 1e-3, seq(0.01, 0.6, by = 0.01))
worst <- 0
for (a in sh) for (b in sh) {
  v <- lih(x, a, b)
  r <- suppressWarnings(stats::pbeta(x, a, b, log.p = TRUE))
  ok <- is.finite(r)
  e <- rel(v[ok], r[ok])
  if (max(e) > 1e-12 || any(!is.finite(v[ok]))) {
    cat(sprintf("a %7g b %7g  worst %.2e at x = %.6g (ref %.6g) nonfinite %d\n",
                a, b, max(e), x[ok][which.max(e)], r[ok][which.max(e)],
                sum(!is.finite(v[ok]))))
  }
  worst <- max(worst, e)
}
cat("overall worst on x <= 0.6:", worst, "\n")
