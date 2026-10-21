# Fix round 1 debugging: which piece of the ex-gaussian survivor carries
# the 1e-4 error in d/dmu at sigma / beta = 1e3.
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(frmtmb)})
ns <- asNamespace("frmtmb")
jac <- function(f, x) as.numeric(MakeTape(f, x)$jacobian(x))
for (z in c(-3, -0.5, 0.5, 3)) {
  cat(sprintf("z %4.1f  d log S: AD %.10f exact %.10f | d log F: AD %.10f exact %.10f\n", z,
              jac(function(p) RTMB::pnorm(p, lower.tail = FALSE, log.p = TRUE), z),
              -dnorm(z) / pnorm(z, lower.tail = FALSE),
              jac(function(p) RTMB::pnorm(p, log.p = TRUE), z),
              dnorm(z) / pnorm(z)))
}
h <- 1e-4
for (t in c(5, 6.5, 8, 999.5)) {
  fd <- (ns$log_mills(t + h) - ns$log_mills(t - h)) / (2 * h)
  cat(sprintf("log_mills t %g: AD %.10g FD %.10g\n", t,
              jac(function(p) ns$log_mills(p), t), fd))
}
f <- function(p) ns$exgaussian_cdf_log(-1.5, p, 1, 1e-3)$lS
cat("exgauss lS d/dmu AD", jac(f, -2), " FD",
    (f(-2 + 1e-6) - f(-2 - 1e-6)) / 2e-6, "\n")
g <- function(p) ns$exgauss_lx(p, 1000)
cat("exgauss_lx d/dz AD", jac(g, 0.501), " FD",
    (g(0.501 + 1e-6) - g(0.501 - 1e-6)) / 2e-6, "\n")
