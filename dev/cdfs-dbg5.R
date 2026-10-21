# Is log_ibeta_half(x, a, b) exact wherever I_x(a, b) <= 0.7, on the
# whole unit interval? Reference: stats::pbeta(log.p = TRUE) where it
# raises no warning, and Rmpfr::pbetaI (whole shapes, 300 bits) where it
# does. Output: dev/cdfs-out/dbg5.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(Rmpfr)})
lih <- asNamespace("frmtmb")$log_ibeta_half
rel <- function(v, r) abs(v - r) / pmax(1, abs(r))
set.seed(5)
sh <- c(0.01, 0.05, 0.3, 0.5, 1, 1.7, 3, 7.3, 25, 120, 900, 7000, 6e4, 5e5)
x <- c(1e-200, 1e-50, 1e-10, 1e-4, 0.01, 0.05, seq(0.1, 0.9, by = 0.05),
       0.95, 0.99, 0.999, 1 - 1e-5, 1 - 1e-8, 1 - 1e-12)
res <- NULL
for (a in sh) for (b in sh) {
  s <- a + b
  m <- a / s
  sd <- sqrt(a * b / (s^2 * (s + 1)))
  xx <- sort(unique(pmin(pmax(c(x, m + c(-30, -6, -2, -1, 0, 1, 2) * sd),
                             1e-300), 1 - 1e-15)))
  for (xi in xx) {
    r <- withCallingHandlers(
      tryCatch(stats::pbeta(xi, a, b, log.p = TRUE), warning = function(w) NA),
      warning = function(w) invokeRestart("muffleWarning"))
    if (is.na(r)) next
    if (!is.finite(r) || r > log(0.7)) next
    res <- rbind(res, data.frame(a = a, b = b, x = xi, lI = r,
                                 e = rel(lih(xi, a, b), r)))
  }
}
cat("stats::pbeta reference, I <= 0.7: points", nrow(res), " worst",
    format(max(res$e), digits = 3), "\n")
print(head(res[order(-res$e), ], 10))
# whole shapes where stats::pbeta underflows or warns
res2 <- NULL
for (ab in list(c(1, 1e5), c(1e5, 1), c(25, 5e4), c(5e4, 25), c(3, 900),
                c(900, 3), c(7000, 7000), c(60000, 500))) {
  a <- ab[1]
  b <- ab[2]
  s <- a + b
  m <- a / s
  sd <- sqrt(a * b / (s^2 * (s + 1)))
  for (xi in sort(unique(pmin(pmax(c(m - c(3000, 300, 40, 5) * sd, 0.05, 0.5,
                                     0.95), 1e-300), 1 - 1e-15)))) {
    r <- as.numeric(log(pbetaI(mpfr(xi, 300), a, b, precBits = 300)))
    if (!is.finite(r) || r > log(0.7)) next
    res2 <- rbind(res2, data.frame(a = a, b = b, x = xi, lI = r,
                                   e = rel(lih(xi, a, b), r)))
  }
}
cat("Rmpfr::pbetaI reference, I <= 0.7: points", nrow(res2), " worst",
    format(max(res2$e), digits = 3), "\n")
print(head(res2[order(-res2$e), ], 6))
