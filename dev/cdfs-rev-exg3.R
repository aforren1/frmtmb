# Reviewer, fix round 1: exgaussian()'s new closed-form density and its
# CDF pair against 400-bit Rmpfr, sigma / beta from 1e-3 to 1e4, y over
# both tails; and the AD gradient of the density against 400-bit
# central differences.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(Rmpfr)})
ns <- asNamespace("frmtmb")
f <- ns$as_frmtmb_family(exgaussian())
ref <- function(y, mu, sg, bt, prec = 400) {
  y <- mpfr(y, prec); sg <- mpfr(sg, prec); bt <- mpfr(bt, prec)
  z <- (y - mu + bt) / sg
  u <- sg / bt
  # log density: -log(beta) + u^2/2 - u z + log Phi(z - u)
  ld <- -log(bt) + u^2 / 2 - u * z + log(pnorm(z - u))
  ex <- exp(u^2 / 2 - u * z) * pnorm(z - u)
  c(as.numeric(ld), as.numeric(log(pnorm(z) - ex)), as.numeric(log(pnorm(-z) + ex)))
}
worst <- c(d = 0, F = 0, S = 0)
for (u in c(1e-3, 1e-2, 0.1, 1, 5, 19, 20, 21, 50, 100, 1e3, 1e4)) {
  sg <- 1; bt <- sg / u; mu <- 0
  sdv <- sqrt(sg^2 + bt^2)
  ys <- mu + sdv * c(-30, -8, -3, -1, 0, 1, 3, 8, 30)
  e <- matrix(0, length(ys), 3)
  for (i in seq_along(ys)) {
    r <- ref(ys[i], mu, sg, bt)
    g <- c(f$lpdf(ys[i], list(mu = mu, sigma = sg, beta = bt), list()),
           unlist(ns$exgaussian_cdf_log(ys[i], mu, sg, bt))[c("lF", "lS")])
    e[i, ] <- abs(g - r) / pmax(1, abs(r))
    e[i, !is.finite(r)] <- NA
  }
  w <- apply(e, 2, max, na.rm = TRUE)
  worst <- pmax(worst, w)
  cat(sprintf("sigma / beta %-6g worst rel log err: density %.2e  lF %.2e  lS %.2e\n", u, w[1], w[2], w[3]))
}
cat(sprintf("ALL: density %.2e  lF %.2e  lS %.2e\n", worst[1], worst[2], worst[3]))
cat("== density gradient d/d(mu, log sigma, log beta) vs 400-bit differences\n")
for (cs in list(c(0.5, 0, -7), c(-1, 0, -3), c(2, 1, 0), c(0, 0, 3), c(0, 0, -9.2))) {
  g <- function(p) f$lpdf(0, list(mu = p[1], sigma = exp(p[2]), beta = exp(p[3])), list())
  ad <- as.numeric(RTMB::MakeTape(g, cs)$jacobian(cs))
  h <- 1e-8
  fd <- sapply(1:3, function(j) {
    e <- replace(numeric(3), j, h)
    (ref(0, cs[1] + e[1], exp(cs[2] + e[2]), exp(cs[3] + e[3]))[1] -
       ref(0, cs[1] - e[1], exp(cs[2] - e[2]), exp(cs[3] - e[3]))[1]) / (2 * h)
  })
  cat(sprintf("mu %g log sigma %g log beta %g (u %.3g): AD %s | ref %s | err %.1e\n", cs[1], cs[2], cs[3],
              exp(cs[2] - cs[3]), paste(signif(ad, 7), collapse = ","), paste(signif(fd, 7), collapse = ","),
              max(abs(ad - fd) / pmax(1, abs(fd)))))
}
