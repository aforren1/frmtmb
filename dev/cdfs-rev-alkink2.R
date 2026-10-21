# Reviewer: asym_laplace's CDF slot on the tape: d/d mu of F(q) and
# lS(q) against central differences, over u = (q - mu) / sigma.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
f <- ns$as_frmtmb_family(asym_laplace())
for (q in c(-3, -1.5, -0.2, 0.2, 2, 4)) {
  g <- function(p) c(f$lcdf(q, list(mu = p[1], sigma = exp(p[2]), quantile = (1 / (1 + exp(-p[3])))), list()),
                     f$lccdf(q, list(mu = p[1], sigma = exp(p[2]), quantile = (1 / (1 + exp(-p[3])))), list()))
  p <- c(0.1, 0, 0)
  ad <- RTMB::MakeTape(g, p)$jacobian(p)
  fd <- numDeriv::jacobian(g, p)
  cat(sprintf("q %4g: F  AD %s | FD %s\n        lS AD %s | FD %s\n", q,
              paste(signif(ad[1, ], 6), collapse = ","), paste(signif(fd[1, ], 6), collapse = ","),
              paste(signif(ad[2, ], 6), collapse = ","), paste(signif(fd[2, ], 6), collapse = ",")))
}
cat("== with the quantile's eta attached as the frame attaches it? the pair directly\n")
for (q in c(-1.5, 2)) {
  g <- function(p) {
    pq <- (1 / (1 + exp(-p[3])))
    r <- ns$asym_laplace_cdf_log(q, p[1], exp(p[2]), pq, list(l = log(pq), l1m = log1p(-pq)))
    c(r$lF, r$lS)
  }
  p <- c(0.1, 0, 0)
  cat("q", q, " AD", signif(RTMB::MakeTape(g, p)$jacobian(p), 6), "\n      FD", signif(numDeriv::jacobian(g, p), 6), "\n")
}
