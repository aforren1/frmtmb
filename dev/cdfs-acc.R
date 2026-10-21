# Accuracy of the new log CDF pairs against independent references,
# both tails, on doubles. The helpers are sourced from the worktree's
# R/families-cdf.R into an environment under the base build's
# namespace, so no install is needed. Output: dev/cdfs-out/acc.txt
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(RTMB)
  library(Rmpfr)
})
wt <- "C:/Users/adf44/source/r/frmtmb-wt-cdfs"
E <- new.env(parent = asNamespace("frmtmb"))
sys.source(file.path(wt, "R/families-cdf.R"), envir = E)
# relative error on the log scale, floored at one, as fams2 measured it
rel <- function(v, r) {
  e <- abs(v - r) / pmax(1, abs(r))
  e[v == r] <- 0
  e
}
report <- function(lab, v, r) {
  ok <- is.finite(r)
  e <- rel(v[ok], r[ok])
  nf <- sum(!is.finite(v[ok]))
  cat(sprintf("%-46s n %5d  worst %.2e  median %.1e  nonfinite %d\n",
              lab, sum(ok), max(e), stats::median(e), nf))
  invisible(max(e))
}
set.seed(20261020)

cat("== log_ibeta_two against stats::pbeta(log.p = TRUE)\n")
sh <- c(0.01, 0.2, 1, 2, 7.3, 40, 300, 5000, 1e5)
x <- c(1e-300, 1e-30, 1e-8, 1e-3, seq(0.01, 0.99, by = 0.01), 1 - 1e-3,
       1 - 1e-8)
wF <- wS <- 0
for (a in sh) for (b in sh) {
  r <- E$log_ibeta_two(x, 1 - x, a, b)
  rF <- suppressWarnings(stats::pbeta(x, a, b, log.p = TRUE))
  rS <- suppressWarnings(stats::pbeta(x, a, b, lower.tail = FALSE,
                                      log.p = TRUE))
  okF <- is.finite(rF)
  okS <- is.finite(rS)
  wF <- max(wF, rel(r$lF[okF], rF[okF]))
  wS <- max(wS, rel(r$lS[okS], rS[okS]))
}
cat(sprintf("shapes %d x %d, x %d points: worst lF %.2e  worst lS %.2e\n",
            length(sh), length(sh), length(x), wF, wS))
# stats::pbeta() is the one off there: at whole shapes beside the worst
# point, 300-bit Rmpfr::pbetaI() referees
for (xx in c(0.99, 0.995)) {
  rr <- as.numeric(log(pbetaI(mpfr(xx, 300), 1e5, 7, precBits = 300)))
  cat(sprintf("x %g shapes (1e5, 7): referee %.12g  log_ibeta_two %.12g  stats %.12g\n",
              xx, rr, E$log_ibeta_two(xx, 1 - xx, 1e5, 7)$lF,
              stats::pbeta(xx, 1e5, 7, log.p = TRUE)))
}

cat("== log_igamma_two against stats::pgamma(log.p = TRUE)\n")
wP <- wQ <- 0
worst_at <- NULL
for (a in c(0.01, 0.3, 1, 3, 12, 20, 30, 300, 3000, 1e4)) {
  xs <- a * c(1e-6, 1e-3, 0.1, 0.3, 0.5, 0.6, 0.65, 0.7, 0.8, 0.9, 1, 1.1,
              1.5, 2, 5, 20, 100) + c(rep(0, 17))
  xs <- sort(unique(c(xs, a + 1 + sqrt(a) * c(-1, 0, 1, 2, 2.5, 3, 4, 8,
                                              30))))
  xs <- xs[xs > 0]
  r <- E$log_igamma_two(xs, a)
  rP <- stats::pgamma(xs, a, log.p = TRUE)
  rQ <- stats::pgamma(xs, a, lower.tail = FALSE, log.p = TRUE)
  okP <- is.finite(rP) & rP > -1e5
  okQ <- is.finite(rQ) & rQ > -1e5
  eP <- rel(r$lP[okP], rP[okP])
  eQ <- rel(r$lQ[okQ], rQ[okQ])
  if (max(eP) > wP) worst_at$P <- c(a = a, x = xs[okP][which.max(eP)])
  if (max(eQ) > wQ) worst_at$Q <- c(a = a, x = xs[okQ][which.max(eQ)])
  wP <- max(wP, eP)
  wQ <- max(wQ, eQ)
}
cat(sprintf("worst lP %.2e at a = %g x = %g; worst lQ %.2e at a = %g x = %g\n",
            wP, worst_at$P[1], worst_at$P[2], wQ, worst_at$Q[1],
            worst_at$Q[2]))

cat("== student-t against stats::pt(log.p = TRUE)\n")
z <- c(-1e6, -1e3, -50, -10, -3, -1, -0.3, -1e-4, 0, 1e-4, 0.3, 1, 3, 10,
       50, 1e3, 1e6)
# stats::pt() is exact through pbeta() to nu = 4e5 and a normal
# approximation above it, so the grid stops there
for (nu in c(1.01, 1.5, 2, 3.7, 4, 10, 30, 1e3, 4e5)) {
  r <- E$student_cdf_log(z, nu)
  report(sprintf("t nu = %g  lF", nu), r$lF, stats::pt(z, nu, log.p = TRUE))
  report(sprintf("t nu = %g  lS", nu), r$lS,
         stats::pt(z, nu, lower.tail = FALSE, log.p = TRUE))
}

cat("== skew-normal: alpha = +-1 is exact, Phi(z)^2 and 1 - Phi(-z)^2\n")
z <- c(-38, -30, -20, -10, -5, -2, -1, -0.1, 0, 0.1, 1, 2, 5, 10, 20, 30,
       38)
sn_q <- function(z, alpha) {
  # mu, sigma such that the standardized argument is z with xi = 0,
  # omega = 1
  delta <- alpha / sqrt(1 + alpha^2)
  om <- 1
  sig <- om * sqrt(1 - 2 * delta^2 / pi)
  mu <- om * delta * sqrt(2 / pi)
  list(q = z, mu = mu, sigma = sig)
}
p1 <- sn_q(z, 1)
r <- E$skew_normal_cdf_log(p1$q, p1$mu, p1$sigma, 1)
report("skew_normal alpha = 1  lF", r$lF, 2 * pnorm(z, log.p = TRUE))
report("skew_normal alpha = 1  lS", r$lS,
       log1p(-exp(2 * pnorm(z, log.p = TRUE))) +
         0 * z)
# the exact upper tail at alpha = 1: S = Phi(-z) (1 + Phi(z))
report("skew_normal alpha = 1  lS (exact form)", r$lS,
       pnorm(-z, log.p = TRUE) + log1p(pnorm(z)))
pm <- sn_q(z, -1)
r <- E$skew_normal_cdf_log(pm$q, pm$mu, pm$sigma, -1)
report("skew_normal alpha = -1 lS", r$lS, 2 * pnorm(-z, log.p = TRUE))
report("skew_normal alpha = -1 lF", r$lF,
       pnorm(z, log.p = TRUE) + log1p(pnorm(-z)))
cat("== skew-normal against sn::psn (body) and angle quadrature (tails)\n")
worst <- 0
for (alpha in c(-50, -5, -2, -0.3, 0, 0.3, 2, 5, 50)) {
  zz <- c(-4, -2, -1, -0.5, 0, 0.5, 1, 2, 4)
  pa <- sn_q(zz, alpha)
  r <- E$skew_normal_cdf_log(pa$q, pa$mu, pa$sigma, alpha)
  ref <- log(sn::psn(zz, 0, 1, alpha))
  refS <- log(1 - sn::psn(zz, 0, 1, alpha))
  ok <- ref > log(1e-6) & refS > log(1e-6)
  worst <- max(worst, rel(r$lF[ok], ref[ok]), rel(r$lS[ok], refS[ok]))
}
cat(sprintf("sn::psn body (both tails above 1e-6): worst %.2e\n", worst))
# The tail on the side of z is G(z, a) = (1 / pi) int_atan(a)^(pi / 2)
# exp(-z^2 / (2 cos^2 t)) dt with a = alpha for z < 0 and -alpha for
# z > 0; adaptive quadrature on the angle, its peak divided out (an
# Rmpfr Romberg rule on the density in z was tried first and was itself
# wrong by 3.8e-5 at z = -6, alpha = 20: the same number for two
# different implementations)
ref_G <- function(z, a) {
  ap <- max(a, 0)
  f <- function(t) exp(-z^2 / 2 * (tan(t)^2 - ap^2))
  tp <- atan(ap)
  i1 <- stats::integrate(f, tp, pi / 2, rel.tol = 1e-13, abs.tol = 0,
                         subdivisions = 5000L)$value
  i2 <- if (a < 0) stats::integrate(f, atan(a), 0, rel.tol = 1e-13,
                                    abs.tol = 0, subdivisions = 5000L)$value
  else 0
  log(i1 + i2) - log(pi) - z^2 * (1 + ap^2) / 2
}
worst <- c(lF = 0, lS = 0)
noref <- 0
cases <- expand.grid(z = c(-38, -30, -12, -6, -3, -1, -0.2),
                     alpha = c(-50, -5, -0.5, 0.5, 3, 20, 300))
for (i in seq_len(nrow(cases))) {
  zz <- cases$z[i]
  al <- cases$alpha[i]
  pa <- sn_q(zz, al)
  v <- E$skew_normal_cdf_log(pa$q, pa$mu, pa$sigma, al)$lF
  rr <- tryCatch(ref_G(zz, al), error = function(e) NA_real_)
  if (!is.finite(rr)) noref <- noref + 1 else worst["lF"] <- max(worst["lF"], rel(v, rr))
  pb <- sn_q(-zz, al)
  v <- E$skew_normal_cdf_log(pb$q, pb$mu, pb$sigma, al)$lS
  rr <- tryCatch(ref_G(-zz, -al), error = function(e) NA_real_)
  if (!is.finite(rr)) noref <- noref + 1 else worst["lS"] <- max(worst["lS"], rel(v, rr))
}
cat(sprintf("angle-quadrature tails, %d points per side (z to +-38, alpha to +-300), %d without a reference: worst lF %.2e  lS %.2e\n",
            nrow(cases), noref, worst["lF"], worst["lS"]))

cat("== ex-gaussian against the closed form in 256-bit arithmetic\n")
mp_exg <- function(q, mu, sigma, beta) {
  q <- mpfr(q, 256)
  u <- mpfr(sigma, 256) / beta
  z <- (q - mu + beta) / sigma
  x <- exp(u * u / 2 - u * z) * pnorm(z - u)
  list(lF = as.numeric(log(pnorm(z) - x)),
       lS = as.numeric(log(pnorm(-z) + x)))
}
worst <- c(lF = 0, lS = 0)
for (sigma in c(0.05, 1, 20)) for (beta in c(0.01, 1, 30)) {
  q <- c(-200, -50, -10, -3, 0, 1, 3, 10, 50, 300)
  v <- E$exgaussian_cdf_log(q, 0, sigma, beta)
  rr <- mp_exg(q, 0, sigma, beta)
  ok <- is.finite(rr$lF) & rr$lF > -1e4
  worst["lF"] <- max(worst["lF"], rel(v$lF[ok], rr$lF[ok]))
  ok <- is.finite(rr$lS) & rr$lS > -1e4
  worst["lS"] <- max(worst["lS"], rel(v$lS[ok], rr$lS[ok]))
}
cat(sprintf("exgaussian 9 (sigma, beta) x 10 q: worst lF %.2e  lS %.2e\n",
            worst["lF"], worst["lS"]))

cat("== asymmetric Laplace: brms::pasym_laplace() lower tail, and the\n")
cat("   closed form in 200-bit arithmetic for both (brms forms the\n")
cat("   upper tail as 1 - F)\n")
mp_al <- function(q, mu, sigma, p) {
  u <- (mpfr(q, 200) - mu) / sigma
  F <- ifelse(as.numeric(u) < 0, p * exp((1 - p) * u),
              1 - (1 - p) * exp(-p * u))
  F <- mpfr(F, 200)
  lF <- lS <- numeric(length(q))
  for (i in seq_along(q)) {
    ui <- u[i]
    if (as.numeric(ui) < 0) {
      lF[i] <- as.numeric(log(mpfr(p, 200)) + (1 - p) * ui)
      lS[i] <- as.numeric(log1p(-p * exp((1 - p) * ui)))
    } else {
      lF[i] <- as.numeric(log1p(-(1 - p) * exp(-p * ui)))
      lS[i] <- as.numeric(log(mpfr(1 - p, 200)) - p * ui)
    }
  }
  list(lF = lF, lS = lS)
}
q <- c(-100, -10, -1, -1e-3, 0, 1e-3, 1, 10, 100)
for (p in c(0.01, 0.3, 0.5, 0.9)) {
  gq <- list(l = log(p), l1m = log1p(-p))
  v <- E$asym_laplace_cdf_log(q, 0, 2, p, gq)
  rr <- mp_al(q, 0, 2, p)
  report(sprintf("asym_laplace p = %g lF vs brms", p), v$lF,
         brms::pasym_laplace(q, 0, 2, p, log.p = TRUE))
  report(sprintf("asym_laplace p = %g lF", p), v$lF, rr$lF)
  report(sprintf("asym_laplace p = %g lS", p), v$lS, rr$lS)
}

cat("== von Mises against adaptive quadrature of the density, its\n")
cat("   largest value on the interval divided out so that the far side\n")
cat("   of a concentrated distribution does not underflow\n")
mp_vm <- function(lo, hi, mu, k) {
  if (hi <= lo) return(-Inf)
  i0e <- besselI(k, 0, expon.scaled = TRUE)
  cand <- c(lo, hi, mu + 2 * pi * (-2:2))
  cand <- cand[cand >= lo & cand <= hi]
  cmax <- max(cos(cand - mu))
  f <- function(t) exp(k * (cos(t - mu) - cmax))
  # split at the peaks inside, so the adaptive rule sees them
  br <- sort(unique(c(lo, hi, cand)))
  tot <- 0
  for (j in seq_len(length(br) - 1L)) {
    tot <- tot + stats::integrate(f, br[j], br[j + 1L], rel.tol = 1e-13,
                                  abs.tol = 0, subdivisions = 2000L)$value
  }
  k * (cmax - 1) + log(tot) - log(2 * pi * i0e)
}
worst <- c(lF = 0, lS = 0)
wh <- list()
for (k in c(0, 0.01, 1, 5, 9, 9.7, 10.5, 30, 300, 3000)) {
  for (mu in c(-3, -1, 0, 0.5, 2.9)) {
    yy <- c(-pi + 1e-3, -3, -2, -1, 0, 1, 2, 3, pi - 1e-3)
    v <- E$von_mises_cdf_log(yy, mu, k)
    for (i in seq_along(yy)) {
      rF <- mp_vm(-pi, yy[i], mu, k)
      rS <- mp_vm(yy[i], pi, mu, k)
      if (is.finite(rF) && rF > -1e4) {
        e <- rel(v$lF[i], rF)
        if (is.na(e)) e <- Inf
        if (e > worst["lF"]) wh$lF <- c(k = k, mu = mu, y = yy[i], ref = rF)
        worst["lF"] <- max(worst["lF"], e)
      }
      if (is.finite(rS) && rS > -1e4) {
        e <- rel(v$lS[i], rS)
        if (is.na(e)) e <- Inf
        if (e > worst["lS"]) wh$lS <- c(k = k, mu = mu, y = yy[i], ref = rS)
        worst["lS"] <- max(worst["lS"], e)
      }
    }
  }
}
cat(sprintf("von_mises: worst lF %.2e at %s\n", worst["lF"],
            paste(names(wh$lF), signif(wh$lF, 4), sep = "=", collapse = " ")))
cat(sprintf("von_mises: worst lS %.2e at %s\n", worst["lS"],
            paste(names(wh$lS), signif(wh$lS, 4), sep = "=", collapse = " ")))

cat("== inverse gaussian against statmod::pinvgauss(log.p = TRUE)\n")
for (mu in c(0.1, 1, 10)) for (shape in c(0.1, 1, 50, 1e3)) {
  q <- mu * c(1e-3, 0.1, 0.5, 1, 2, 5, 20, 100)
  v <- E$inverse_gaussian_cdf_log(q, mu, shape)
  report(sprintf("inv.gaussian mu %g shape %g lF", mu, shape), v$lF,
         statmod::pinvgauss(q, mu, shape, log.p = TRUE))
  report(sprintf("inv.gaussian mu %g shape %g lS", mu, shape), v$lS,
         statmod::pinvgauss(q, mu, shape, lower.tail = FALSE, log.p = TRUE))
}

cat("== negative binomial against stats::pnbinom(log.p = TRUE)\n")
q <- c(-1, 0, 1, 3, 10, 30, 100, 1000)
for (mu in c(1e-6, 0.3, 4, 100, 1e4)) for (shape in c(0.05, 1, 2, 7.5, 1e3)) {
  v <- E$nbinom_cdf_log(q, log(mu), log(shape))
  report(sprintf("nbinom mu %g shape %g lF", mu, shape), v$lF,
         stats::pnbinom(q, size = shape, mu = mu, log.p = TRUE))
  report(sprintf("nbinom mu %g shape %g lS", mu, shape), v$lS,
         stats::pnbinom(q, size = shape, mu = mu, lower.tail = FALSE,
                        log.p = TRUE))
}

cat("== binomial against stats::pbinom(log.p = TRUE)\n")
for (size in c(1, 10, 1000)) for (p in c(1e-8, 0.05, 0.5, 0.97)) {
  q <- c(-1, 0, 1, round(size / 3), size - 1, size, size + 3)
  v <- E$binomial_cdf_log(q, size, list(p = p, q = 1 - p))
  report(sprintf("binomial size %g p %g lF", size, p), v$lF,
         stats::pbinom(q, size, p, log.p = TRUE))
  report(sprintf("binomial size %g p %g lS", size, p), v$lS,
         stats::pbinom(q, size, p, lower.tail = FALSE, log.p = TRUE))
}

cat("== poisson against stats::ppois(log.p = TRUE)\n")
for (mu in c(1e-6, 0.5, 5, 80, 3000)) {
  q <- c(-1, 0, 1, 3, 10, 60, 200, 700, 5000)
  v <- E$poisson_cdf_log(q, mu)
  report(sprintf("poisson mu %g lF", mu), v$lF, stats::ppois(q, mu,
                                                            log.p = TRUE))
  report(sprintf("poisson mu %g lS", mu), v$lS,
         stats::ppois(q, mu, lower.tail = FALSE, log.p = TRUE))
}

cat("== beta-binomial against a 200-bit sum of its mass function\n")
mp_bb <- function(q, size, a, b) {
  j <- mpfr(0:size, 200)
  lp <- lgamma(j + a) + lgamma(size - j + b) - lgamma(mpfr(size + a + b,
                                                         200)) +
    lgamma(mpfr(a + b, 200)) - lgamma(mpfr(a, 200)) - lgamma(mpfr(b, 200)) +
    lgamma(mpfr(size + 1, 200)) - lgamma(j + 1) - lgamma(size - j + 1)
  p <- exp(lp)
  sapply(q, function(qq) {
    lo <- if (qq < 0) mpfr(0, 200) else sum(p[seq_len(min(qq, size) + 1)])
    hi <- if (qq >= size) mpfr(0, 200) else sum(p[(max(qq, -1) + 2):(size + 1)])
    c(as.numeric(log(lo)), as.numeric(log(hi)))
  })
}
for (size in c(1, 7, 60)) for (ab in list(c(0.3, 0.3), c(2, 9), c(400, 30))) {
  q <- c(-1, 0, 1, floor(size / 2), size - 1, size)
  v <- E$beta_binomial_cdf_log(q, size, ab[1], ab[2])
  rr <- mp_bb(q, size, ab[1], ab[2])
  report(sprintf("beta_binomial size %g a %g b %g lF", size, ab[1], ab[2]),
         v$lF, rr[1, ])
  report(sprintf("beta_binomial size %g a %g b %g lS", size, ab[1], ab[2]),
         v$lS, rr[2, ])
}
