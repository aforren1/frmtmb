# Reviewer: AD gradient of the skew-normal and ex-gaussian log CDF
# pairs against Richardson differences of the same functions on
# doubles, over z and the shape; and RTMB::pnorm(log.p = TRUE)'s own
# derivative in the far lower tail against the Mills ratio.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
chk <- function(lab, fun, p) {
  tp <- RTMB::MakeTape(function(pp) fun(pp), p)
  ad <- as.numeric(tp$jacobian(p))
  fd <- numDeriv::jacobian(function(pp) fun(pp), p, method = "Richardson")
  e <- max(abs(ad - fd) / pmax(1, abs(fd)))
  cat(sprintf("%-52s AD %s | FD %s | err %.1e%s\n", lab,
              paste(signif(ad, 7), collapse = ","), paste(signif(fd, 7), collapse = ","),
              e, if (e > 1e-5) "  <-- WRONG" else ""))
}
cat("== skew_normal lF and lS: d/d(mu, log sigma, alpha) at q = 0\n")
for (al in c(-50, -10, -5, -2, 0.5, 2, 5, 10, 50)) {
  for (mu in c(-3, -1, 0, 1, 3)) {
    f1 <- function(p) ns$skew_normal_cdf_log(0, p[1], exp(p[2]), p[3])$lF
    f2 <- function(p) ns$skew_normal_cdf_log(0, p[1], exp(p[2]), p[3])$lS
    chk(sprintf("skew_normal lF alpha %g mu %g", al, mu), f1, c(mu, 0, al))
    chk(sprintf("skew_normal lS alpha %g mu %g", al, mu), f2, c(mu, 0, al))
  }
}
cat("== owen_g_log pieces at alpha = -10: d/dz, d/da\n")
for (z in c(-2, -1, -0.5, -0.2, 0.2, 0.5, 1, 2)) {
  chk(sprintf("owen_g_log z %g a 10", z), function(p) ns$owen_g_log(p[1], p[2]), c(z, 10))
  chk(sprintf("owen_g_log z %g a -10", z), function(p) ns$owen_g_log(p[1], p[2]), c(z, -10))
  chk(sprintf("owen_t_ad z %g a -10", z), function(p) ns$owen_t_ad(p[1], p[2]), c(z, -10))
  chk(sprintf("owen_t_ad z %g a 10", z), function(p) ns$owen_t_ad(p[1], p[2]), c(z, 10))
}
cat("== RTMB::pnorm(w, log.p = TRUE) derivative vs the Mills ratio\n")
for (w in c(-5, -10, -20, -37, -38, -40, -100, -1e3, -2e4)) {
  tp <- RTMB::MakeTape(function(p) RTMB::pnorm(p, log.p = TRUE), w)
  ad <- as.numeric(tp$jacobian(w))
  ex <- exp(dnorm(w, log = TRUE) - pnorm(w, log.p = TRUE))
  cat(sprintf("w %-7g value %.10g (stats %.10g) d/dw AD %.10g exact %.10g\n", w,
              tp(w), pnorm(w, log.p = TRUE), ad, ex))
}
cat("== exgaussian lF, lS: d/d(mu, log sigma, log beta) at q = 0\n")
for (ls in c(-3, 0, 3, 6, 10)) for (lb in c(-10, -3, 0, 3)) {
  f1 <- function(p) ns$exgaussian_cdf_log(0, p[1], exp(p[2]), exp(p[3]))$lF
  f2 <- function(p) ns$exgaussian_cdf_log(0, p[1], exp(p[2]), exp(p[3]))$lS
  for (mu in c(-2, 2)) {
    chk(sprintf("exgauss lF log sigma %g log beta %g mu %g", ls, lb, mu), f1, c(mu, ls, lb))
    chk(sprintf("exgauss lS log sigma %g log beta %g mu %g", ls, lb, mu), f2, c(mu, ls, lb))
  }
}
