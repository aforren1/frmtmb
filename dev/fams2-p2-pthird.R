# Punch 2, n1: log_pbeta_ad(), the pbeta branch of log_ibeta_half(),
# over the region it is read in (shapes 150 to 1e7, x within 6.5 sd of
# the mean; seed 11, 20000 points, as dev/fams2-p1-pthird.R), then at
# the mean exactly, at the tie of its second form, (a + 1) / (a + b + 1),
# and at each end of its blend. Derivatives up to third, in (x, a, b);
# value against stats::pbeta(log.p = TRUE), relative, floored at one.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
f <- frmtmb:::log_pbeta_ad
F <- MakeTape(function(p) f(p[1], p[2], p[3]), c(0.3, 300, 700))
J <- F$jacfun(); H <- J$jacfun(); T3 <- H$jacfun()
ok <- function(p) all(is.finite(J(p))) && all(is.finite(H(p))) &&
  all(is.finite(T3(p)))
err <- function(p) {
  r <- pbeta(p[1], p[2], p[3], log.p = TRUE)
  abs(F(p) - r) / max(1, abs(r))
}
set.seed(11)
M <- 20000
a <- exp(runif(M, log(150), log(1e7)))
b <- exp(runif(M, log(150), log(1e7)))
s <- a + b
sd <- sqrt(a * b / (s * s * (s + 1)))
x <- a / s + runif(M, -6.5, 6.5) * sd
bad <- 0; e <- numeric(M)
for (i in seq_len(M)) {
  p <- c(x[i], a[i], b[i])
  if (!ok(p)) bad <- bad + 1
  e[i] <- err(p)
}
cat("random: points", M, " non-finite derivatives", bad,
    " value error max", signif(max(e), 3), "\n")
cat("special points (non-finite derivatives, value error):\n")
for (ab in list(c(300, 700), c(150, 150), c(151, 1e7), c(2e4, 8e4),
                c(1e7, 1e7), c(268.941421, 731.058579), c(1e7, 150))) {
  A <- ab[1]; B <- ab[2]; S <- A + B; gap <- B / (S * (S + 1))
  pts <- c(mean = A / S, tie2 = (A + 1) / (S + 1),
           lo = A / S + 0.25 * gap, hi = A / S + 0.75 * gap)
  for (nm in names(pts)) {
    p <- c(pts[[nm]], A, B)
    cat(sprintf("  a %-9g b %-9g %-4s finite %s  err %.1e\n", A, B, nm,
                ok(p), err(p)))
  }
}
