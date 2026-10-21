# Fix round 1 debugging: the von Mises interval on the far side, the
# pieces log_window() reads.
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
k <- 50
mu <- -2
y <- c(0.6, 1.0, 1.154)
y2 <- y + 0.4
p1 <- ns$von_mises_cdf_log(y, mu, k)
p2 <- ns$von_mises_cdf_log(y2, mu, k)
print(rbind(lF1 = p1$lF, lS1 = p1$lS, lF2 = p2$lF, lS2 = p2$lS))
f <- ns$as_frmtmb_family(von_mises())
print(ns$row_lpdf(f, y, y, list(mu = rep(mu, 3), kappa = rep(k, 3)),
                  list(cens = rep(2, 3), cens_y2 = y2), NULL))
dens <- function(t) exp(k * (cos(t - mu) - 1)) /
  (2 * pi * besselI(k, 0, expon.scaled = TRUE))
print(log(mapply(function(a, b) integrate(dens, a, b, rel.tol = 1e-12)$value,
                 y, y2)))
