.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-lib", "C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
ns <- asNamespace("frmtmb")
"c" <- RTMB::ADoverload("c")
for (q in c(-30, -4, -1, -0.2, 0.3, 1.5, 6, 40)) {
  f <- function(p) { lp <- ns$student_cdf_log((q - p[1]) / exp(p[2]), 1 + exp(p[3])); c(lp$lF, lp$lS) }
  p0 <- c(0.3, 0.2, log(3))
  F <- MakeTape(f, p0)
  g <- F$jacobian(p0); r <- numDeriv::jacobian(f, p0)
  r2 <- numDeriv::jacobian(function(p) c(stats::pt((q - p[1]) / exp(p[2]), 1 + exp(p[3]), log.p = TRUE), stats::pt((q - p[1]) / exp(p[2]), 1 + exp(p[3]), lower.tail = FALSE, log.p = TRUE)), p0)
  cat("q", q, "\n"); print(rbind(ad = as.vector(g), fd_self = as.vector(r), fd_ref = as.vector(r2)))
}
