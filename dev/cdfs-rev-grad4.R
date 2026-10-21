# Reviewer: (1) ex-gaussian AD derivatives against central differences
# of the 320-bit closed form (so the reference has no cancellation);
# (2) skew-normal AD derivatives over a fine q sweep at alpha = -10 and
# +10, against Richardson differences on doubles.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(Rmpfr)
})
ns <- asNamespace("frmtmb")
exg_mp <- function(q, mu, lsig, lbet, prec = 400) {
  sg <- exp(mpfr(lsig, prec))
  bt <- exp(mpfr(lbet, prec))
  z <- (mpfr(q, prec) - mu + bt) / sg
  u <- sg / bt
  ex <- exp(u^2 / 2 - u * z) * pnorm(z - u)
  c(log(pnorm(z) - ex), log(pnorm(-z) + ex))
}
cat("== ex-gaussian at q = 0: AD vs 400-bit central differences (h = 1e-8)\n")
for (cs in list(c(-2, 0, -10), c(2, 0, -10), c(-2, 3, -10), c(-2, 6, -10),
                c(-2, 10, 0), c(2, 0, -5), c(-2, 0, -7), c(0, 0, -3), c(-2, 2, -2))) {
  tp <- RTMB::MakeTape(function(p) {
    r <- ns$exgaussian_cdf_log(0, p[1], exp(p[2]), exp(p[3]))
    c(r$lF, r$lS)
  }, cs)
  ad <- tp$jacobian(cs)
  h <- 1e-8
  fd <- sapply(1:3, function(j) {
    e <- replace(numeric(3), j, h)
    a <- exg_mp(0, cs[1] + e[1], cs[2] + e[2], cs[3] + e[3])
    b <- exg_mp(0, cs[1] - e[1], cs[2] - e[2], cs[3] - e[3])
    as.numeric((a - b) / (2 * h))
  })
  e <- max(abs(ad - fd) / pmax(1e-300, abs(fd)))
  cat(sprintf("mu %g log sigma %g log beta %g (u = %.3g): worst rel err %.2e\n", cs[1], cs[2],
              cs[3], exp(cs[2] - cs[3]), e))
  cat("   AD lF:", signif(ad[1, ], 6), "  ref:", signif(fd[1, ], 6), "\n")
  cat("   AD lS:", signif(ad[2, ], 6), "  ref:", signif(fd[2, ], 6), "\n")
}
cat("== skew-normal: AD vs Richardson over q, mu 0 sigma 1\n")
for (al in c(-10, -3, 3, 10)) {
  qs <- seq(-6, 6, by = 0.1)
  worst <- 0; wq <- NA
  for (q in qs) {
    f <- function(p) {
      r <- ns$skew_normal_cdf_log(q, p[1], exp(p[2]), p[3])
      c(r$lF, r$lS)
    }
    p <- c(0, 0, al)
    ad <- RTMB::MakeTape(f, p)$jacobian(p)
    fd <- numDeriv::jacobian(f, p, method = "Richardson")
    e <- max(abs(ad - fd) / pmax(1, abs(fd)))
    if (e > worst) { worst <- e; wq <- q; wad <- ad; wfd <- fd }
  }
  cat(sprintf("alpha %g: worst rel err %.2e at q = %g\n", al, worst, wq))
  if (worst > 1e-6) { print(signif(wad, 7)); print(signif(wfd, 7)) }
}
