# Who is right at shapes (1e5, 7), x = 0.99: log_ibeta_half() or
# stats::pbeta()? 300-bit Rmpfr::pbetaI as the referee.
.libPaths(c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(Rmpfr)})
lih <- asNamespace("frmtmb")$log_ibeta_half
for (x in c(0.99, 0.995)) {
  r <- as.numeric(log(pbetaI(mpfr(x, 300), 1e5, 7, precBits = 300)))
  cat(sprintf("x %g ref %.12g lih %.12g stats %.12g\n", x, r, lih(x, 1e5, 7), stats::pbeta(x, 1e5, 7, log.p = TRUE)))
}
