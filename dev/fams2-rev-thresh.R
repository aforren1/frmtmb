# Reviewer: two thresholds in the shape sum s = a + b.
#  1. where RTMB::dbeta(log = TRUE)'s gradient first goes NaN (the
#     interior of Beta(), zero_inflated_beta(), zoib and xbeta)
#  2. where log_ibeta_half()'s 50-step error at the switch point first
#     exceeds 1e-10 and 1e-6 (worst over x within +-1e-2 relative of m
#     and r = a / s in 0.02..0.48)
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
lib <- frmtmb:::log_ibeta_half
D <- MakeTape(function(p) RTMB::dbeta(p[1], p[2], p[3], log = TRUE),
              c(0.3, 1, 1))
DJ <- D$jacfun()
ss <- round(10^seq(3, 4, by = 0.02))
cat("== 1. RTMB::dbeta gradient finite? (x = mean and x = mean +- 2 sd) ==\n")
first_bad <- NA
for (s in ss) {
  ok <- TRUE
  for (r in c(0.05, 0.3, 0.5)) {
    a <- r * s; b <- (1 - r) * s; sd <- sqrt(r * (1 - r) / (s + 1))
    for (x in r + c(-2, 0, 2) * sd) {
      if (!all(is.finite(DJ(c(x, a, b))))) ok <- FALSE
    }
  }
  if (!ok && is.na(first_bad)) first_bad <- s
}
cat("first shape sum with a NaN gradient:", first_bad, "\n")
for (s in c(first_bad - 1, first_bad, 2000, 2500, 3000)) {
  cat(sprintf("  s %d: grad at (0.3, 0.3 s, 0.7 s) %s\n", s,
              paste(format(DJ(c(0.3, 0.3 * s, 0.7 * s)), digits = 4),
                    collapse = " ")))
}
cat("\n== 2. log_ibeta_half() error at the switch, by shape sum ==\n")
for (s in c(1e3, 2e3, 3e3, 5e3, 7e3, 1e4, 2e4, 5e4, 1e5)) {
  worst <- 0
  for (r in seq(0.02, 0.48, by = 0.02)) {
    a <- r * s; b <- (1 - r) * s; m <- (a + 1) / (s + 2)
    for (dx in seq(-1e-2, 1e-2, length.out = 41)) {
      x <- m * (1 + dx)
      if (x >= 0.5) next
      ref <- pbeta(x, a, b, log.p = TRUE)
      e <- abs(lib(x, a, b) - ref) / max(1, abs(ref))
      worst <- max(worst, e)
    }
  }
  cat(sprintf("  s %.0e: worst error %.2e\n", s, worst))
}
