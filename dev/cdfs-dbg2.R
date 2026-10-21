# log_ibeta_half() above x = 0.6: where does it hold? Error against
# Rmpfr::pbetaI in 200 bits, by the position of x relative to the mean
.libPaths(c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(Rmpfr)})
lih <- asNamespace("frmtmb")$log_ibeta_half
rel <- function(v, r) { e <- abs(v - r) / pmax(1, abs(r)); e[v == r] <- 0; e }
res <- NULL
for (a in c(0.5, 3, 30, 300, 3000, 5e4, 5e5)) for (b in c(0.05, 0.5, 2, 25, 300)) {
  s <- a + b; m <- a / s; sd <- sqrt(a * b / (s^2 * (s + 1)))
  x <- c(0.61, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999, m - c(40, 20, 10, 5, 2, 1, 0) * sd, m + c(1, 3) * sd)
  x <- x[x > 0.6 & x < 1]
  if (!length(x)) next
  v <- lih(x, a, b)
  r <- stats::pbeta(x, a, b, log.p = TRUE)
  res <- rbind(res, data.frame(a = a, b = b, x = x, k = (x - m) / sd, F = r, e = rel(v, r)))
}
res <- res[is.finite(res$F), ]
cat("points", nrow(res), "\n")
for (cut in c(-30, -10, -5, -2, 0, 2, 10)) {
  sub <- res[res$k <= cut, ]
  cat(sprintf("x at or below mean %+g sd: n %3d worst %.2e\n", cut, nrow(sub), if (nrow(sub)) max(sub$e) else NA))
}
print(head(res[order(-res$e), ], 15))
