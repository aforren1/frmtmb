.libPaths(c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
E <- new.env(parent = asNamespace("frmtmb")); sys.source("R/families-cdf.R", envir = E)
rel <- function(v, r) { e <- abs(v - r) / pmax(1, abs(r)); e[v == r] <- 0; e }
sh <- c(0.01, 0.2, 1, 2, 7.3, 40, 300, 5000, 1e5)
x <- c(1e-300, 1e-30, 1e-8, 1e-3, seq(0.01, 0.99, by = 0.01), 1 - 1e-3, 1 - 1e-8)
for (a in sh) for (b in sh) {
  r <- E$log_ibeta_two(x, 1 - x, a, b)
  rF <- suppressWarnings(stats::pbeta(x, a, b, log.p = TRUE))
  rS <- suppressWarnings(stats::pbeta(x, a, b, lower.tail = FALSE, log.p = TRUE))
  eF <- rel(r$lF, rF); eS <- rel(r$lS, rS)
  eF[!is.finite(rF)] <- 0; eS[!is.finite(rS)] <- 0
  if (max(eF) > 1e-12) cat(sprintf("lF a %g b %g worst %.2e at x=%.10g v=%.6g ref=%.6g\n", a, b, max(eF), x[which.max(eF)], r$lF[which.max(eF)], rF[which.max(eF)]))
  if (max(eS) > 1e-12) cat(sprintf("lS a %g b %g worst %.2e at x=%.10g v=%.6g ref=%.6g\n", a, b, max(eS), x[which.max(eS)], r$lS[which.max(eS)], rS[which.max(eS)]))
}
z <- c(-1e6, -1e3, -50, -10, -3, -1, -0.3, -1e-4, 0, 1e-4, 0.3, 1, 3, 10, 50, 1e3, 1e6)
r <- E$student_cdf_log(z, 1e6); print(cbind(z, r$lF, stats::pt(z, 1e6, log.p=TRUE)))
print(E$skew_normal_cdf_log(c(-38,-30,-20,-10,-5,5,10,20,30,38), sqrt(2/pi)/sqrt(2), sqrt(1-1/pi), 1)$lF - 2*pnorm(c(-38,-30,-20,-10,-5,5,10,20,30,38), log.p=TRUE))
