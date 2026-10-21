# Reviewer, fix round 1: poisson_cdf_log()'s three forms (q = 0, sums to
# q = 30, incomplete gamma above) and nbinom_cdf_log()'s q = 0 closed
# form, at and around the switch points, values against stats' log
# tails and AD derivatives against the analytic d/d log mu.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
qs <- c(0, 1, 2, 29, 30, 31, 32)
wv <- 0; wg <- 0; where <- ""
for (mu in c(1e-10, 1e-3, 0.5, 5, 15, 25, 30, 31, 45, 60, 100, 1e3)) {
  for (q in qs) {
    p <- ns$poisson_cdf_log(q, mu)
    rF <- ppois(q, mu, log.p = TRUE); rS <- ppois(q, mu, lower.tail = FALSE, log.p = TRUE)
    e <- max(abs(c(p$lF - rF, p$lS - rS)) / pmax(1, abs(c(rF, rS))))
    if (!is.finite(e)) e <- Inf
    if (e > wv) { wv <- e; where <- sprintf("q %g mu %g: lF %.15g (%.15g) lS %.15g (%.15g)", q, mu, p$lF, rF, p$lS, rS) }
    tp <- RTMB::MakeTape(function(lm) { r <- ns$poisson_cdf_log(q, exp(lm)); c(r$lF, r$lS) }, log(mu))
    ad <- as.numeric(tp$jacobian(log(mu)))
    ld <- log(mu) + dpois(q, mu, log = TRUE)
    an <- c(-exp(ld - rF), exp(ld - rS))
    eg <- max(abs(ad - an) / pmax(1, abs(an)))
    if (!is.finite(eg)) eg <- Inf
    if (eg > wg) { wg <- eg; wherg <- sprintf("q %g mu %g: AD %s analytic %s", q, mu,
                                              paste(signif(ad, 10), collapse = ","), paste(signif(an, 10), collapse = ","))}
  }
}
cat(sprintf("poisson values: worst rel log err %.2e at %s\n", wv, where))
cat(sprintf("poisson d/d log mu: worst rel err %.2e at %s\n", wg, wherg))
cat("== the survivor's blend region (1 <= q <= 30, S near 1e-3) and the 40-term sum's reach\n")
for (cs in list(c(10, 3), c(30, 15), c(30, 25), c(20, 1e-3), c(30, 29), c(5, 1))) {
  q <- cs[1]; mu <- cs[2]
  p <- ns$poisson_cdf_log(q, mu)
  cat(sprintf("q %2g mu %-6g lS %.15g ref %.15g rel %.1e\n", q, mu, p$lS,
              ppois(q, mu, lower.tail = FALSE, log.p = TRUE),
              abs(p$lS / ppois(q, mu, lower.tail = FALSE, log.p = TRUE) - 1)))
}
cat("== nbinom at q = 0 and 1\n")
w0 <- 0
for (sh in c(1e-4, 0.5, 2, 1e3, 1e8)) for (mu in c(1e-8, 0.3, 5, 1e4)) for (q in c(0, 1)) {
  p <- ns$nbinom_cdf_log(q, log(mu), log(sh))
  rF <- pnbinom(q, size = sh, mu = mu, log.p = TRUE)
  rS <- pnbinom(q, size = sh, mu = mu, lower.tail = FALSE, log.p = TRUE)
  e <- max(abs(c(p$lF - rF, p$lS - rS)) / pmax(1, abs(c(rF, rS))))
  if (e > w0) { w0 <- e; wh <- sprintf("q %g shape %g mu %g", q, sh, mu) }
}
cat(sprintf("nbinom q in {0, 1}: worst rel log err %.2e at %s\n", w0, wh))
