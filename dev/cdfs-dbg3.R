# log_ibeta_half() on the student-t tail with x above 0.6 (a = nu / 2,
# b = 1/2), against stats::pbeta(log.p = TRUE) (TOMS 708), by the value
# of I_x. Output: dev/cdfs-out/dbg3.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
lih <- asNamespace("frmtmb")$log_ibeta_half
rel <- function(v, r) abs(v - r) / pmax(1, abs(r))
res <- NULL
for (nu in c(1.01, 1.7, 3, 8, 30, 200, 2e3, 2e4, 4e5)) {
  z <- c(seq(0.02, 3, by = 0.02), 3.5, 4, 5, 7, 10, 20, 30, 40, 60, 100,
         300) * c(1, rep(1, 165))[1]
  z <- c(z, sqrt(nu) * c(0.1, 0.3, 0.5, 0.7, 0.8))
  x <- nu / (nu + z^2)
  keep <- x > 0.6
  x <- x[keep]
  z <- z[keep]
  r <- stats::pbeta(x, nu / 2, 0.5, log.p = TRUE)
  v <- lih(x, nu / 2, 0.5)
  res <- rbind(res, data.frame(nu = nu, z = z, x = x, lI = r, e = rel(v, r)))
}
res <- res[is.finite(res$lI), ]
cat("points", nrow(res), "\n")
for (cut in c(-690, -100, -10, log(0.01), log(0.1), log(0.3), log(0.6),
              log(0.9), 0)) {
  sub <- res[res$lI <= cut, ]
  cat(sprintf("I_x <= %-12.4g n %4d worst %.2e\n", exp(cut), nrow(sub),
              if (nrow(sub)) max(sub$e) else NA))
}
print(head(res[order(-res$e), ], 12))
