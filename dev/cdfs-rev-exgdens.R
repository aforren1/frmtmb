# Reviewer: is the ex-gaussian AD error at large sigma / beta also in
# the density (RTMBdist::dexgauss, pre-existing), or only in the CDF?
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(Rmpfr)})
ns <- asNamespace("frmtmb")
f <- ns$as_frmtmb_family(exgaussian())
dmp <- function(y, mu, ls, lb, prec = 400) {
  sg <- exp(mpfr(ls, prec)); bt <- exp(mpfr(lb, prec))
  lam <- 1 / bt; m <- mu - bt
  # log density of N(m, sg) + Exp(lam)
  as.numeric(log(lam) + lam * (m - y) + lam^2 * sg^2 / 2 + log(pnorm((y - m) / sg - lam * sg)))
}
for (cs in list(c(-2, 0, -10), c(2, 0, -7), c(-2, 0, -3))) {
  g <- function(p) f$lpdf(0, list(mu = p[1], sigma = exp(p[2]), beta = exp(p[3])), list())
  ad <- RTMB::MakeTape(g, cs)$jacobian(cs)
  h <- 1e-8
  fd <- sapply(1:3, function(j) {
    e <- replace(numeric(3), j, h)
    (dmp(0, cs[1] + e[1], cs[2] + e[2], cs[3] + e[3]) - dmp(0, cs[1] - e[1], cs[2] - e[2], cs[3] - e[3])) / (2 * h)
  })
  cat(sprintf("lpdf mu %g log sigma %g log beta %g (u = %.3g): AD %s | ref %s\n", cs[1], cs[2], cs[3],
              exp(cs[2] - cs[3]), paste(signif(ad, 6), collapse = ","), paste(signif(fd, 6), collapse = ",")))
}
cat("== the same against Richardson differences of f$lpdf on doubles, and of brms::dexgaussian\n")
for (cs in list(c(-2, 0, -10), c(2, 0, -7), c(-2, 0, -3), c(-2, 0, 0))) {
  g <- function(p) f$lpdf(0, list(mu = p[1], sigma = exp(p[2]), beta = exp(p[3])), list())
  gb <- function(p) brms::dexgaussian(0, p[1], exp(p[2]), exp(p[3]), log = TRUE)
  ad <- RTMB::MakeTape(g, cs)$jacobian(cs)
  cat(sprintf("u = %-8.3g AD %s | FD(lpdf) %s | FD(brms) %s | mpfr value %.10g lpdf %.10g brms %.10g\n",
              exp(cs[2] - cs[3]), paste(signif(ad, 6), collapse = ","),
              paste(signif(numDeriv::grad(g, cs), 6), collapse = ","),
              paste(signif(numDeriv::grad(gb, cs), 6), collapse = ","),
              dmp(0, cs[1], cs[2], cs[3]), g(cs), gb(cs)))
}
