# log_ibeta_half() above x = 0.6 for whole shapes, against
# Rmpfr::pbetaI (exact for whole shapes). Output: dev/cdfs-out/dbg4.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(Rmpfr)})
lih <- asNamespace("frmtmb")$log_ibeta_half
cat("lih at x = 1 - 2^-53:", lih(1 - 2^-53, 3, 0.5), " at x = 1:", lih(1, 3, 0.5), "\n")
for (ab in list(c(50000, 25), c(3000, 25), c(300, 2), c(40, 1), c(500, 300), c(20, 20), c(5, 1), c(2, 3))) {
  a <- ab[1]; b <- ab[2]; s <- a + b; m <- a / s; sd <- sqrt(a * b / (s^2 * (s + 1)))
  x <- sort(unique(c(0.61, 0.7, 0.8, 0.9, 0.95, m - c(50, 10, 3, 1) * sd, m, m + sd)))
  x <- x[x > 0.6 & x < 1]
  for (xx in x) {
    r <- as.numeric(log(pbetaI(mpfr(xx, 300), a, b, precBits = 300)))
    v <- lih(xx, a, b); st <- stats::pbeta(xx, a, b, log.p = TRUE)
    cat(sprintf("a %6g b %4g x %.6f k %8.1f  ref %.10g  lih rel %.1e  stats rel %.1e\n", a, b, xx, (xx - m) / sd, r,
                abs(v - r) / max(1, abs(r)), abs(st - r) / max(1, abs(r))))
  }
}
