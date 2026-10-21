# Reviewer, fix round 1: row_lpdf()'s log window on rows in an upper
# tail and at a point mass, against references that do not share its
# code (Rmpfr, mass sums, stats' log tails), and its AD gradient against
# Richardson differences of the same function on doubles.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(Rmpfr)})
ns <- asNamespace("frmtmb")
F <- function(f) ns$as_frmtmb_family(f)
rl <- function(fam, y, dp, av) ns$row_lpdf(F(fam), y, y, dp, av, NULL)
lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
out <- function(lab, got, ref) {
  e <- abs(got - ref) / pmax(1, abs(ref))
  cat(sprintf("%-58s got %-16.10g ref %-16.10g rel %.1e%s\n", lab, got, ref, e,
              if (!is.finite(got) || e > 1e-10) "  <--" else ""))
}
mp <- function(x) mpfr(x, 300)
cat("== gaussian (mu 0, sigma 1): interval and truncation windows, 300-bit pnorm\n")
for (w in list(c(8, 9), c(20, 20.5), c(37, 38), c(-9, -8), c(-0.1, 0.1))) {
  ref <- as.numeric(if (w[1] > 0) log(pnorm(-mp(w[1])) - pnorm(-mp(w[2]))) else log(pnorm(mp(w[2])) - pnorm(mp(w[1]))))
  out(sprintf("gaussian interval (%g, %g]", w[1], w[2]),
      rl(gaussian(), w[1], list(mu = 0, sigma = 1), list(cens = 2, cens_y2 = w[2])), ref)
  # a truncated row at the window's middle: lpdf - log mass
  ym <- mean(w)
  ref2 <- dnorm(ym, log = TRUE) - ref
  out(sprintf("gaussian trunc(lb = %g, ub = %g) row", w[1], w[2]),
      rl(gaussian(), ym, list(mu = 0, sigma = 1), list(trunc_lb = w[1], trunc_ub = w[2])), ref2)
}
for (lb in c(10, 30, 38)) {
  ref <- dnorm(lb + 0.01, log = TRUE) - pnorm(lb, lower.tail = FALSE, log.p = TRUE)
  out(sprintf("gaussian trunc(lb = %g) alone", lb),
      rl(gaussian(), lb + 0.01, list(mu = 0, sigma = 1), list(trunc_lb = lb)), ref)
}
cat("== left-censored row below 1e-308 (gaussian, student)\n")
out("gaussian left at -40", rl(gaussian(), -40, list(mu = 0, sigma = 1), list(cens = -1)),
    pnorm(-40, log.p = TRUE))
out("student nu 5 left at -1e60", rl(student(), -1e60, list(mu = 0, sigma = 1, nu = 5),
                                      list(cens = -1)), pt(-1e60, 5, log.p = TRUE))
cat("== Gamma (mu 2, shape 3): upper-tail interval, Rmpfr igamma\n")
for (w in list(c(20, 21), c(40, 40.5), c(0.001, 0.002))) {
  a <- 3; s <- 1.5   # rate shape / mu
  Q <- function(x) igamma(mp(a), mp(x) * s) / gamma(mp(a))
  ref <- as.numeric(log(Q(w[1]) - Q(w[2])))
  out(sprintf("Gamma interval (%g, %g]", w[1], w[2]),
      rl(Gamma(link = "log"), w[1], list(mu = 2, shape = 3), list(cens = 2, cens_y2 = w[2])), ref)
}
cat("== poisson and negbinomial: intervals and windows by mass sums\n")
for (w in list(c(40, 45), c(80, 81), c(2, 9), c(0, 3))) {
  ref <- lse(dpois((w[1] + 1):w[2], 5, log = TRUE))
  out(sprintf("poisson mu 5 interval (%g, %g]", w[1], w[2]),
      rl(poisson(), w[1], list(mu = 5), list(cens = 2, cens_y2 = w[2])), ref)
  refn <- lse(dnbinom((w[1] + 1):w[2], mu = 5, size = 2, log = TRUE))
  out(sprintf("negbinomial mu 5 shape 2 interval (%g, %g]", w[1], w[2]),
      rl(negbinomial(), w[1], list(mu = 5, shape = 2), list(cens = 2, cens_y2 = w[2])), refn)
}
for (lb in c(1, 30, 60)) {
  ref <- dpois(lb, 5, log = TRUE) - ppois(lb - 1, 5, lower.tail = FALSE, log.p = TRUE)
  out(sprintf("poisson mu 5 trunc(lb = %g), row at lb", lb),
      rl(poisson(), lb, list(mu = 5), list(trunc_lb = lb)), ref)
}
cat("== point masses\n")
hu <- 0.3; zi <- 0.3
for (y2 in c(1e-6, 0.5, 30)) {
  ref <- log1p(-hu) + pgamma(y2, 3, 1.5, log.p = TRUE)
  out(sprintf("hurdle_gamma interval (0, %g]", y2),
      rl(hurdle_gamma(), 0, list(mu = 2, shape = 3, hu = hu), list(cens = 2, cens_y2 = y2)), ref)
}
out("hurdle_gamma left at 0 (the mass)", rl(hurdle_gamma(), 0, list(mu = 2, shape = 3, hu = 1e-12),
                                            list(cens = -1)), log(1e-12))
out("hurdle_gamma right at 0", rl(hurdle_gamma(), 0, list(mu = 2, shape = 3, hu = hu),
                                  list(cens = 1)), log1p(-hu))
for (w in list(c(0, 3), c(0, 40), c(-1, 0), c(-1, 2))) {
  k <- (w[1] + 1):w[2]
  lm <- ifelse(k == 0, log(zi + (1 - zi) * exp(-5)), log1p(-zi) + dpois(k, 5, log = TRUE))
  out(sprintf("zero_inflated_poisson interval (%g, %g]", w[1], w[2]),
      rl(zero_inflated_poisson(), w[1], list(mu = 5, zi = zi), list(cens = 2, cens_y2 = w[2])),
      lse(lm))
}
out("zero_inflated_poisson trunc(lb = 1), row 1", rl(zero_inflated_poisson(), 1,
    list(mu = 5, zi = zi), list(trunc_lb = 1)),
    log1p(-zi) + dpois(1, 5, log = TRUE) - log((1 - zi) * (1 - exp(-5))))
k <- 1:3
out("hurdle_poisson interval (0, 3]", rl(hurdle_poisson(), 0, list(mu = 5, hu = hu),
    list(cens = 2, cens_y2 = 3)),
    log1p(-hu) + lse(dpois(k, 5, log = TRUE)) - log1p(-exp(-5)))
out("zero_inflated_beta interval (0, 0.3]", rl(zero_inflated_beta(), 0,
    list(mu = 0.4, phi = 5, zi = zi), list(cens = 2, cens_y2 = 0.3)),
    log1p(-zi) + pbeta(0.3, 2, 3, log.p = TRUE))
cat("== von Mises far-side arc and a short arc at -pi (integrate, max divided out)\n")
vmref <- function(a, b, mu, k) {
  ld <- function(t) k * (cos(t - mu) - 1) - log(2 * pi * besselI(k, 0, TRUE))
  m <- max(ld(c(a, b)))
  m + log(integrate(function(t) exp(ld(t) - m), a, b, rel.tol = 1e-13)$value)
}
for (cs in list(c(2.5, 3, 0, 30), c(-3.14, -3.1, 0, 30), c(-pi, -pi + 1e-9, 0.5, 2),
                c(2.9, 3.1, 0, 300))) {
  out(sprintf("von_mises interval (%g, %g] mu %g kappa %g", cs[1], cs[2], cs[3], cs[4]),
      rl(von_mises(), cs[1], list(mu = cs[3], kappa = cs[4]), list(cens = 2, cens_y2 = cs[2])),
      vmref(cs[1], cs[2], cs[3], cs[4]))
}
cat("== AD vs Richardson of row_lpdf on doubles, interval rows in a tail\n")
gchk <- function(lab, fam, y, av, mk) {
  f <- function(p) sum(ns$row_lpdf(F(fam), y, y, mk(p), av, NULL))
  p0 <- attr(mk, "p0")
  ad <- as.numeric(RTMB::MakeTape(f, p0)$jacobian(p0))
  fd <- numDeriv::grad(f, p0)
  e <- max(abs(ad - fd) / pmax(1, abs(fd)))
  cat(sprintf("%-46s fn %-12.6g AD %s | FD %s | err %.1e%s\n", lab, f(p0),
              paste(signif(ad, 7), collapse = ","), paste(signif(fd, 7), collapse = ","), e,
              if (e > 1e-6) "  <--" else ""))
}
yy <- c(8, 20, 3, -9, 0.2); av <- list(cens = c(2, 2, 1, -1, 0), cens_y2 = c(9, 20.5, 3, -9, 0.2))
mkg <- function(p) list(mu = p[1], sigma = exp(p[2])); attr(mkg, "p0") <- c(0.1, 0.05)
gchk("gaussian tail intervals", gaussian(), yy, av, mkg)
mks <- function(p) list(mu = p[1], sigma = exp(p[2]), alpha = p[3]); attr(mks, "p0") <- c(0, 0, -5)
gchk("skew_normal tail intervals", skew_normal(), c(3, 4, 2, -3, 0), list(cens = c(2, 2, 1, -1, 0),
     cens_y2 = c(3.5, 4.5, 2, -3, 0)), mks)
mkh <- function(p) list(mu = exp(p[1]), shape = exp(p[2]), hu = 1 / (1 + exp(-p[3])))
attr(mkh, "p0") <- c(0.7, 4, -1)
gchk("hurdle_gamma (0, y2] at shape 55", hurdle_gamma(), c(0, 0, 0.5, 3),
     list(cens = c(2, 2, -1, 1), cens_y2 = c(0.5, 1, 0.5, 3)), mkh)
mkz <- function(p) list(mu = exp(p[1]), zi = 1 / (1 + exp(-p[2]))); attr(mkz, "p0") <- c(1.5, -1)
gchk("zero_inflated_poisson intervals and edges", zero_inflated_poisson(), c(0, 0, 30, 2, 0),
     list(cens = c(2, -1, 2, 1, 0), cens_y2 = c(3, 0, 35, 2, 0)), mkz)
mkv <- function(p) list(mu = 2 * atan(p[1]), kappa = exp(p[2])); attr(mkv, "p0") <- c(0, 4)
gchk("von_mises far-side intervals", von_mises(), c(2.5, -3.1, 1, -1),
     list(cens = c(2, 2, 1, -1), cens_y2 = c(3, -3, 1, -1)), mkv)
mkp <- function(p) list(mu = exp(p[1])); attr(mkp, "p0") <- log(5)
gchk("poisson trunc(lb = 30) rows", poisson(), c(30, 31, 35), list(trunc_lb = 30), mkp)
gchk("poisson intervals (40, 45], trunc(ub = 50)", poisson(), c(40, 3, 1),
     list(cens = c(2, 1, -1), cens_y2 = c(45, 3, 1), trunc_ub = 50), mkp)
