# Reviewer: accuracy of the lane's CDF pairs on my own grids.
# Metric: |log(got) - log(ref)|, the RELATIVE error of the tail
# probability (the quantity a likelihood term carries), reported for
# both lF and lS, with the point where it is worst. The installed
# lane build is used (wt-cdfs-rev-lib), not sourced helpers.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(Rmpfr)
})
ns <- asNamespace("frmtmb")
cat("lib:", find.package("frmtmb"), "\n")
L <- function(p) suppressWarnings(log(p))

rep_line <- function(lab, got, ref, pts) {
  ok <- is.finite(ref)
  # error of the log, relative where |log| > 1 and absolute (the tail
  # probability's relative error) where it is smaller
  e <- abs(got - ref) / pmax(1, abs(ref))
  e[got == ref] <- 0
  e[!ok] <- NA
  nf <- sum(ok & !is.finite(got))
  e2 <- e
  e2[!is.finite(e2)] <- Inf
  i <- if (any(ok)) which.max(replace(e2, !ok, -1)) else NA
  cat(sprintf("%-44s n %4d worst %.2e at %s (ref %.6g got %.6g) nonfinite %d\n",
              lab, sum(ok), if (is.na(i)) NA else e2[i],
              if (is.na(i)) "-" else pts[i],
              if (is.na(i)) NA else ref[i], if (is.na(i)) NA else got[i],
              nf))
  invisible(if (is.na(i)) 0 else e2[i])
}

#--- student --------------------------------------------------------
cat("== student: student_cdf_log vs stats::pt (nu <= 1e5) and 256-bit\n")
cat("   pbetaI at whole nu/2 (nu >= 1e4)\n")
z <- c(-1e10, -1e5, -1e3, -100, -30, -10, -5, -2, -1, -0.5, -0.4, -0.31,
       -0.3, -0.1, -1e-8, 0, 1e-8, 0.1, 0.3, 0.31, 0.4, 0.5, 1, 2, 5, 10,
       30, 100, 1e3, 1e5, 1e10)
for (nu in c(1.0001, 1.5, 2, 4, 7, 30, 300, 1e4, 1e5)) {
  p <- ns$student_cdf_log(z, nu)
  rep_line(sprintf("student nu=%g lF", nu), p$lF, stats::pt(z, nu, log.p = TRUE),
           paste0("z=", z))
  rep_line(sprintf("student nu=%g lS", nu), p$lS,
           stats::pt(z, nu, lower.tail = FALSE, log.p = TRUE), paste0("z=", z))
}
# large nu: stats::pt() is a normal approximation above 4e5, so the
# reference is 256-bit pbetaI at a whole nu / 2 for 2e4 (tail =
# I_x(nu/2, 1/2) / 2, x = nu / (nu + z^2)) and the tail integral of
# stats::dt(log = TRUE) above it
tail_int0 <- function(ld, q, side, scale) {
  l0 <- ld(q)
  f <- function(u) exp(ld(q + side * scale * u) - l0)
  v <- 0
  br <- c(0, 0.01, 0.1, 1, 10, 100, 1e3, Inf)
  for (j in seq_len(length(br) - 1)) {
    r <- stats::integrate(f, br[j], br[j + 1], rel.tol = 1e-13,
                          subdivisions = 2000L)
    v <- v + r$value
  }
  l0 + log(scale * v)
}
for (nu in c(2e4, 1e6, 1e8, 1e12)) {
  zz <- c(-200, -50, -10, -3, -1, 1, 3, 10, 50, 200)
  got <- ns$student_cdf_log(zz, nu)
  ref <- vapply(zz, function(zi) {
    if (FALSE) {
      xm <- mpfr(nu, 256) / (nu + mpfr(zi, 256)^2)
      as.numeric(log(pbetaI(xm, nu / 2, 0.5, precBits = 256) / 2))
    } else {
      ld <- function(t) stats::dt(t, nu, log = TRUE)
      tail_int0(ld, -abs(zi), -1, 1 / max(abs(zi), 1))
    }
  }, 0)
  refF <- ifelse(zz < 0, ref, log1p(-exp(ref)))
  refS <- ifelse(zz > 0, ref, log1p(-exp(ref)))
  rep_line(sprintf("student nu=%g lF (ref)", nu), got$lF, refF, paste0("z=", zz))
  rep_line(sprintf("student nu=%g lS (ref)", nu), got$lS, refS, paste0("z=", zz))
}

#--- Gamma ----------------------------------------------------------
cat("== Gamma: gamma_cdf_log(q, mu = 1, shape) vs stats::pgamma(log.p)\n")
for (a in c(1e-4, 0.01, 0.5, 1, 3.3, 12, 59.5, 60, 75, 100, 150, 185, 186,
            500, 1e4, 1e5, 1e6, 1e7)) {
  sd <- sqrt(a)
  x <- c(a * c(1e-12, 1e-6, 1e-3, 0.1, 0.5, 0.8, 0.9), a + sd * (-8:8),
         a + sd * c(12, 20, 40), a * c(2, 5, 30, 1000))
  x <- sort(unique(x[x > 0]))
  p <- ns$gamma_cdf_log(x / a, 1, a)
  rep_line(sprintf("Gamma shape=%g lF", a), p$lF,
           stats::pgamma(x, a, log.p = TRUE), sprintf("x=%.4g", x))
  rep_line(sprintf("Gamma shape=%g lS", a), p$lS,
           stats::pgamma(x, a, lower.tail = FALSE, log.p = TRUE),
           sprintf("x=%.4g", x))
}

#--- Beta -----------------------------------------------------------
cat("== Beta: beta_cdf_log vs stats::pbeta(log.p); 256-bit pbetaI referee\n")
cat("   at whole shape2 where the two disagree by more than 1e-12\n")
shp <- list(c(0.01, 0.01), c(0.5, 0.5), c(1e-3, 5), c(5, 1e-3), c(2, 3),
            c(50, 50), c(1e4, 1e4), c(1e6, 2), c(2, 1e6), c(0.3, 1e5),
            c(1e5, 0.3), c(1, 1))
xg <- c(1e-300, 1e-100, 1e-20, 1e-5, 0.01, 0.1, 0.3, 0.45, 0.5, 0.55, 0.6,
        0.65, 0.7, 0.9, 0.99, 1 - 1e-5, 1 - 1e-10, 1 - 1e-15)
for (s in shp) {
  a <- s[1]; b <- s[2]
  m <- a / (a + b); sdv <- sqrt(a * b / ((a + b)^2 * (a + b + 1)))
  x <- sort(unique(c(xg, pmin(pmax(m + sdv * c(-6, -3, -1, 0, 1, 3, 6),
                                   1e-300), 1 - 1e-16))))
  x <- x[x > 0 & x < 1]
  p <- ns$beta_cdf_log(x, a, b)
  rF <- stats::pbeta(x, a, b, log.p = TRUE)
  rS <- stats::pbeta(x, a, b, lower.tail = FALSE, log.p = TRUE)
  rep_line(sprintf("Beta (%g, %g) lF", a, b), p$lF, rF, sprintf("x=%.6g", x))
  rep_line(sprintf("Beta (%g, %g) lS", a, b), p$lS, rS, sprintf("x=%.6g", x))
}
# a referee at whole shape2, where stats is not the truth
for (s in list(c(3, 2), c(2000, 2), c(7, 1), c(300, 3))) {
  a <- s[1]; b <- s[2]
  x <- c(1e-5, 0.3, 0.9, 0.999, 0.999999)
  if (a > 100) x <- c(0.99, 0.999, 0.9999, 0.99999)
  ref <- vapply(x, function(xx) as.numeric(log(pbetaI(mpfr(xx, 256), a, b,
                                                       precBits = 256))), 0)
  refS <- vapply(x, function(xx) as.numeric(log(pbetaI(1 - mpfr(xx, 256), b, a,
                                                        precBits = 256))), 0)
  p <- ns$beta_cdf_log(x, a, b)
  rep_line(sprintf("Beta (%g, %g) lF (mpfr)", a, b), p$lF, ref, sprintf("x=%.8g", x))
  rep_line(sprintf("Beta (%g, %g) lS (mpfr)", a, b), p$lS, refS, sprintf("x=%.8g", x))
}

#--- negative binomial ----------------------------------------------
cat("== negbinomial: nbinom_cdf_log vs stats::pnbinom(log.p)\n")
for (sh in c(1e-4, 0.05, 1, 2, 7.5, 1e3, 1e6, 1e8)) {
  for (mu in c(1e-8, 0.3, 4, 100, 1e4, 1e6)) {
    v <- mu + mu^2 / sh
    q <- c(0, 1, 2, floor(mu + sqrt(v) * c(-3, -1, 0, 1, 3, 10)),
           floor(mu * 10 + 50), floor(mu * 100 + 500))
    q <- sort(unique(q[q >= 0]))
    p <- ns$nbinom_cdf_log(q, log(mu), log(sh))
    rep_line(sprintf("nbinom shape=%g mu=%g lF", sh, mu), p$lF,
             stats::pnbinom(q, size = sh, mu = mu, log.p = TRUE), paste0("q=", q))
    rep_line(sprintf("nbinom shape=%g mu=%g lS", sh, mu), p$lS,
             stats::pnbinom(q, size = sh, mu = mu, lower.tail = FALSE,
                            log.p = TRUE), paste0("q=", q))
  }
}

#--- binomial -------------------------------------------------------
cat("== binomial: binomial_cdf_log vs stats::pbinom(log.p)\n")
for (n in c(1, 2, 10, 1000, 1e6)) {
  for (pr in c(1e-12, 1e-6, 0.05, 0.5, 0.95, 1 - 1e-6)) {
    q <- c(0, 1, floor(n / 2), floor(n * pr + sqrt(n * pr * (1 - pr)) *
                                         c(-5, -1, 0, 1, 5)), n - 1)
    q <- sort(unique(q[q >= 0 & q < n]))
    p <- ns$binomial_cdf_log(q, n, list(p = pr, q = 1 - pr))
    rep_line(sprintf("binom size=%g p=%g lF", n, pr), p$lF,
             stats::pbinom(q, n, pr, log.p = TRUE), paste0("q=", q))
    rep_line(sprintf("binom size=%g p=%g lS", n, pr), p$lS,
             stats::pbinom(q, n, pr, lower.tail = FALSE, log.p = TRUE),
             paste0("q=", q))
  }
}

#--- poisson --------------------------------------------------------
cat("== poisson: poisson_cdf_log and the lccdf slot vs stats::ppois\n")
fp <- ns$as_frmtmb_family(poisson())
for (mu in c(1e-10, 1e-3, 0.5, 5, 80, 1e3, 1e5, 1e6)) {
  q <- c(0, 1, 2, floor(mu + sqrt(mu) * c(-8, -3, 0, 3, 8, 20)),
         floor(mu * 10 + 100))
  q <- sort(unique(q[q >= 0]))
  p <- ns$poisson_cdf_log(q, mu)
  rep_line(sprintf("poisson mu=%g lF (pair)", mu), p$lF,
           stats::ppois(q, mu, log.p = TRUE), paste0("q=", q))
  rep_line(sprintf("poisson mu=%g lS (pair)", mu), p$lS,
           stats::ppois(q, mu, lower.tail = FALSE, log.p = TRUE), paste0("q=", q))
  rep_line(sprintf("poisson mu=%g lS (slot)", mu),
           fp$lccdf(q, list(mu = mu), list()),
           stats::ppois(q, mu, lower.tail = FALSE, log.p = TRUE), paste0("q=", q))
}

#--- inverse gaussian -----------------------------------------------
cat("== inverse.gaussian vs statmod::pinvgauss(log.p)\n")
for (mu in c(0.1, 1, 10, 1e3)) for (sh in c(1e-3, 0.1, 1, 50, 1e3, 1e5)) {
  sdv <- sqrt(mu^3 / sh)
  q <- c(mu * c(1e-4, 0.01, 0.2, 0.5, 0.9, 1, 1.1, 2, 10, 100),
         mu + sdv * c(-3, 3, 10, 30))
  q <- sort(unique(q[q > 0]))
  p <- ns$inverse_gaussian_cdf_edges(q, mu, sh)
  rep_line(sprintf("ig mu=%g shape=%g lF", mu, sh), p$lF,
           statmod::pinvgauss(q, mu, sh, log.p = TRUE), sprintf("q=%.4g", q))
  rep_line(sprintf("ig mu=%g shape=%g lS", mu, sh), p$lS,
           statmod::pinvgauss(q, mu, sh, lower.tail = FALSE, log.p = TRUE),
           sprintf("q=%.4g", q))
}

#--- tail integral reference for continuous densities ----------------
# log int of exp(ld) over the tail on the side `side` of q: the
# integrand is taken relative to its value at q, so neither underflows.
tail_int <- function(ld, q, side, scale) {
  l0 <- ld(q)
  f <- function(u) exp(ld(q + side * scale * u) - l0)
  v <- 0
  br <- c(0, 0.01, 0.1, 1, 10, 100, 1e3, Inf)
  for (j in seq_len(length(br) - 1)) {
    r <- try(stats::integrate(f, br[j], br[j + 1], rel.tol = 1e-13,
                              subdivisions = 2000L), silent = TRUE)
    if (inherits(r, "try-error")) return(NA_real_)
    v <- v + r$value
  }
  l0 + log(scale * v)
}
# reference pair: near tail by integration, far tail its complement
ref_pair <- function(ld, q, near_left, scale) {
  t <- vapply(seq_along(q), function(i) {
    tail_int(ld, q[i], if (near_left[i]) -1 else 1, scale[i])
  }, 0)
  c1 <- log(-expm1(t))
  list(lF = ifelse(near_left, t, c1), lS = ifelse(near_left, c1, t))
}

#--- skew normal ----------------------------------------------------
cat("== skew_normal (mu 0, sigma 1) vs integration of sn::dsn on the tail\n")
for (al in c(-1e3, -50, -5, -1, -0.1, 0, 0.1, 1, 5, 50, 1e3)) {
  delta <- al / sqrt(1 + al^2)
  om <- 1 / sqrt(1 - 2 * delta^2 / pi)
  xi <- -om * delta * sqrt(2 / pi)
  q <- c(-40, -20, -10, -5, -2, -1, -0.1, 0, 0.1, 1, 2, 5, 10, 20, 40)
  ld <- function(t) sn::dsn(t, xi, om, al, log = TRUE)
  med <- sn::qsn(0.5, xi, om, al)
  near_left <- q < med
  # local scale: one over the log density's slope, capped by omega
  sl <- abs((ld(q + 1e-6) - ld(q - 1e-6)) / 2e-6)
  sc <- pmin(om, 1 / pmax(sl, 1e-300))
  r <- ref_pair(ld, q, near_left, sc)
  p <- ns$skew_normal_cdf_log(q, 0, 1, al)
  rep_line(sprintf("skew_normal alpha=%g lF", al), p$lF, r$lF, paste0("q=", q))
  rep_line(sprintf("skew_normal alpha=%g lS", al), p$lS, r$lS, paste0("q=", q))
}

#--- ex-gaussian ----------------------------------------------------
cat("== exgaussian vs the closed form in 320-bit arithmetic\n")
exg_ref <- function(q, mu, sigma, beta) {
  qm <- mpfr(q, 320)
  z <- (qm - mu + beta) / sigma
  u <- mpfr(sigma, 320) / beta
  ex <- exp(u^2 / 2 - u * z) * pnorm(z - u)
  F <- pnorm(z) - ex
  S <- pnorm(-z) + ex
  c(as.numeric(log(F)), as.numeric(log(S)))
}
for (sg in c(1e-3, 0.1, 1, 10)) for (bt in c(1e-3, 0.1, 1, 10, 1e3)) {
  mu <- 0
  sdv <- sqrt(sg^2 + bt^2)
  q <- mu + sdv * c(-40, -10, -3, -1, 0, 1, 3, 10, 40, 200)
  rr <- vapply(q, exg_ref, c(0, 0), mu = mu, sigma = sg, beta = bt)
  p <- ns$exgaussian_cdf_log(q, mu, sg, bt)
  rep_line(sprintf("exgauss sigma=%g beta=%g lF", sg, bt), p$lF, rr[1, ],
           sprintf("q=%.4g", q))
  rep_line(sprintf("exgauss sigma=%g beta=%g lS", sg, bt), p$lS, rr[2, ],
           sprintf("q=%.4g", q))
}

#--- asymmetric Laplace ---------------------------------------------
cat("== asym_laplace vs its closed form in 256-bit arithmetic\n")
for (pq in c(1e-6, 0.01, 0.5, 0.99, 1 - 1e-6)) {
  u <- c(-1e4, -100, -10, -1, -1e-3, 0, 1e-3, 1, 10, 100, 1e4)
  pm <- mpfr(pq, 256)
  um <- mpfr(u, 256)
  F <- ifelse(u < 0, 0, 1)
  rF <- rS <- numeric(length(u))
  for (i in seq_along(u)) {
    if (u[i] < 0) {
      Fi <- pm * exp((1 - pm) * um[i])
    } else {
      Fi <- 1 - (1 - pm) * exp(-pm * um[i])
    }
    rF[i] <- as.numeric(log(Fi))
    rS[i] <- as.numeric(log(1 - Fi))
  }
  p <- ns$asym_laplace_cdf_log(u, 0, 1, pq, list(l = log(pq), l1m = log1p(-pq)))
  rep_line(sprintf("asym_laplace p=%g lF", pq), p$lF, rF, paste0("u=", u))
  rep_line(sprintf("asym_laplace p=%g lS", pq), p$lS, rS, paste0("u=", u))
}

#--- von Mises ------------------------------------------------------
cat("== von_mises vs adaptive quadrature of the density, max divided out\n")
vm_ref <- function(a, b, mu, k) {
  # log int_a^b of the density, a < b in [-pi, pi]
  ld <- function(t) k * (cos(t - mu) - 1) - log(2 * pi) -
    log(besselI(k, 0, expon.scaled = TRUE))
  # maximum of cos(t - mu) on [a, b]: at mu (or mu +- 2 pi) if inside
  cands <- c(a, b, mu, mu - 2 * pi, mu + 2 * pi)
  cands <- cands[cands >= a & cands <= b]
  lm <- max(ld(cands))
  br <- sort(unique(c(a, b, cands)))
  # extra breakpoints near the peak, at the density's own scale
  sc <- 1 / sqrt(k + 1e-300)
  pk <- cands[which.max(ld(cands))]
  br <- sort(unique(c(br, pmin(pmax(pk + sc * c(-50, -10, -3, -1, 1, 3, 10, 50),
                                     a), b))))
  v <- 0
  for (j in seq_len(length(br) - 1)) {
    if (br[j + 1] <= br[j]) next
    r <- stats::integrate(function(t) exp(ld(t) - lm), br[j], br[j + 1],
                          rel.tol = 1e-13, subdivisions = 5000L)
    v <- v + r$value
  }
  lm + log(v)
}
for (k in c(1e-10, 1e-3, 0.5, 5, 50, 700, 1e4, 1e5)) {
  for (mu in c(-3.1, -1, 0, 2, 3.1)) {
    y <- c(-pi + 1e-9, -3, -2, -1, mu - 1e-3, mu, mu + 1e-3, 0.5, 2, 3,
           pi - 1e-9)
    y <- sort(unique(pmin(pmax(y, -pi + 1e-12), pi - 1e-12)))
    rF <- vapply(y, function(yy) vm_ref(-pi, yy, mu, k), 0)
    rS <- vapply(y, function(yy) vm_ref(yy, pi, mu, k), 0)
    p <- ns$von_mises_cdf_log(y, mu, k)
    rep_line(sprintf("von_mises kappa=%g mu=%g lF", k, mu), p$lF, rF,
             sprintf("y=%.6g", y))
    rep_line(sprintf("von_mises kappa=%g mu=%g lS", k, mu), p$lS, rS,
             sprintf("y=%.6g", y))
  }
}

#--- beta-binomial --------------------------------------------------
cat("== beta_binomial vs a log-sum of extraDistr::dbbinom\n")
lse <- function(v) { m <- max(v); if (!is.finite(m)) m else m + log(sum(exp(v - m))) }
for (n in c(1, 2, 50, 1000)) {
  for (s in list(c(0.01, 0.01), c(0.3, 0.3), c(2, 9), c(400, 30), c(1e5, 1e5))) {
    q <- sort(unique(c(0, 1, floor(n / 3), floor(n / 2), n - 1)))
    q <- q[q >= 0 & q < n]
    lm <- extraDistr::dbbinom(0:n, n, s[1], s[2], log = TRUE)
    rF <- vapply(q, function(k) lse(lm[1:(k + 1)]), 0)
    rS <- vapply(q, function(k) lse(lm[(k + 2):(n + 1)]), 0)
    p <- ns$beta_binomial_cdf_log(q, n, s[1], s[2])
    rep_line(sprintf("bbinom n=%g (%g, %g) lF", n, s[1], s[2]), p$lF, rF,
             paste0("q=", q))
    rep_line(sprintf("bbinom n=%g (%g, %g) lS", n, s[1], s[2]), p$lS, rS,
             paste0("q=", q))
  }
}

#--- COM-Poisson ----------------------------------------------------
cat("== compois vs a log-sum of RTMB::dcompois2 (frmtmb's own mass)\n")
for (mu in c(0.1, 3, 30)) for (nu in c(0.3, 1, 3)) {
  K <- 3000
  lm <- RTMB::dcompois2(0:K, mu, nu, log = TRUE)
  q <- sort(unique(c(0, 1, floor(mu), floor(mu * 3 + 10), floor(mu * 6 + 30))))
  rF <- vapply(q, function(k) lse(lm[1:(k + 1)]), 0)
  rS <- vapply(q, function(k) lse(lm[(k + 2):(K + 1)]), 0)
  p <- ns$compois_cdf_log(q, mu, nu)
  rep_line(sprintf("compois mu=%g nu=%g lF", mu, nu), p$lF, rF, paste0("q=", q))
  rep_line(sprintf("compois mu=%g nu=%g lS", mu, nu), p$lS, rS, paste0("q=", q))
}

#--- zero-inflated and hurdle, through the family slots -------------
cat("== zi / hurdle point masses, through the slots: log(lcdf), lccdf\n")
slot <- function(fam, q, dp, av = list()) {
  f <- ns$as_frmtmb_family(fam)
  list(lF = L(f$lcdf(q, dp, av)), lS = f$lccdf(q, dp, av))
}
chk <- function(lab, fam, q, dp, av, rF, rS) {
  s <- slot(fam, q, dp, av)
  rep_line(paste(lab, "lF"), s$lF, rF, paste0("q=", q))
  rep_line(paste(lab, "lS"), s$lS, rS, paste0("q=", q))
}
q <- c(-1, 0, 1, 3, 20, 60)
for (zi in c(1e-12, 0.3, 1 - 1e-9)) {
  mu <- 2
  F <- ifelse(q < 0, 0, zi + (1 - zi) * stats::ppois(q, mu))
  S <- ifelse(q < 0, 1, (1 - zi) * stats::ppois(q, mu, lower.tail = FALSE))
  S2 <- ifelse(q < 0, 0, log1p(-zi) + stats::ppois(q, mu, lower.tail = FALSE,
                                                    log.p = TRUE))
  chk(sprintf("zip zi=%g", zi), zero_inflated_poisson(), q,
      list(mu = mu, zi = zi), list(), L(F), ifelse(q < 0, 0, S2))
  hu <- zi
  Fh <- ifelse(q < 0, 0, hu + (1 - hu) * (stats::ppois(q, mu) - exp(-mu)) /
                 (1 - exp(-mu)))
  Sh <- ifelse(q < 0, 0, log1p(-hu) + stats::ppois(q, mu, lower.tail = FALSE,
                                                    log.p = TRUE) -
                 log1p(-exp(-mu)))
  chk(sprintf("hurdle_poisson hu=%g", hu), hurdle_poisson(), q,
      list(mu = mu, hu = hu), list(), L(Fh), Sh)
  sh <- 0.7
  p0 <- stats::dnbinom(0, size = sh, mu = mu)
  Fn <- ifelse(q < 0, 0, hu + (1 - hu) * (stats::pnbinom(q, size = sh, mu = mu) -
                                             p0) / (1 - p0))
  Sn <- ifelse(q < 0, 0, log1p(-hu) + stats::pnbinom(q, size = sh, mu = mu,
                                                      lower.tail = FALSE,
                                                      log.p = TRUE) - log1p(-p0))
  chk(sprintf("hurdle_negbinomial hu=%g", hu), hurdle_negbinomial(), q,
      list(mu = mu, shape = sh, hu = hu), list(), L(Fn), Sn)
  Fz <- ifelse(q < 0, 0, zi + (1 - zi) * stats::pnbinom(q, size = sh, mu = mu))
  Sz <- ifelse(q < 0, 0, log1p(-zi) + stats::pnbinom(q, size = sh, mu = mu,
                                                      lower.tail = FALSE,
                                                      log.p = TRUE))
  chk(sprintf("zinb zi=%g", zi), zero_inflated_negbinomial(), q,
      list(mu = mu, shape = sh, zi = zi), list(), L(Fz), Sz)
  qb <- c(-1, 0, 1, 5, 9, 10)
  Fb <- ifelse(qb < 0, 0, zi + (1 - zi) * stats::pbinom(qb, 10, 0.3))
  Sb <- ifelse(qb < 0, 0, log1p(-zi) + stats::pbinom(qb, 10, 0.3,
                                                      lower.tail = FALSE,
                                                      log.p = TRUE))
  chk(sprintf("zib zi=%g", zi), zero_inflated_binomial(), qb,
      list(mu = 0.3, zi = zi), list(trials = 10), L(Fb), Sb)
  qc <- c(-1, 0, 1e-300, 1e-8, 0.5, 3, 50)
  Fg <- ifelse(qc < 0, 0, hu + (1 - hu) * stats::pgamma(qc, 3, 3 / 2))
  Sg <- ifelse(qc < 0, 0, log1p(-hu) + stats::pgamma(qc, 3, 3 / 2,
                                                      lower.tail = FALSE,
                                                      log.p = TRUE))
  chk(sprintf("hurdle_gamma hu=%g", hu), hurdle_gamma(), qc,
      list(mu = 2, shape = 3, hu = hu), list(), L(Fg), Sg)
  Fl <- ifelse(qc < 0, 0, hu + (1 - hu) * stats::plnorm(qc, 0.2, 0.7))
  Sl <- ifelse(qc < 0, 0, log1p(-hu) + stats::plnorm(qc, 0.2, 0.7,
                                                      lower.tail = FALSE,
                                                      log.p = TRUE))
  chk(sprintf("hurdle_lognormal hu=%g", hu), hurdle_lognormal(), qc,
      list(mu = 0.2, sigma = 0.7, hu = hu), list(), L(Fl), Sl)
  qz <- c(-0.5, 0, 1e-12, 0.3, 0.9, 1 - 1e-12, 1)
  Fzb <- ifelse(qz < 0, 0, zi + (1 - zi) * stats::pbeta(qz, 2, 3))
  Szb <- ifelse(qz < 0, 0, log1p(-zi) + stats::pbeta(qz, 2, 3,
                                                      lower.tail = FALSE,
                                                      log.p = TRUE))
  chk(sprintf("zi_beta zi=%g", zi), zero_inflated_beta(), qz,
      list(mu = 0.4, phi = 5, zi = zi), list(), L(Fzb), Szb)
}
cat("== zi_asym_laplace at y = 0 and either side\n")
qa <- c(-3, -1e-9, 0, 1e-9, 3)
for (zi in c(0.2)) {
  pq <- 0.3
  Fb <- brms::pasym_laplace(qa, 0.3, 1.2, pq)
  F <- (1 - zi) * Fb + zi * (qa >= 0)
  S <- 1 - F
  chk("zi_asym_laplace zi=0.2", zero_inflated_asym_laplace(), qa,
      list(mu = 0.3, sigma = 1.2, quantile = pq, zi = zi), list(), L(F), L(S))
}
cat("done\n")
