# Gradients of the new log CDF pairs on the tape: the AD gradient
# against Richardson finite differences (numDeriv) of an independent
# reference where one exists (stats, statmod), else of the same
# function on doubles; and the second and third derivatives' finiteness
# at random points, which the Laplace approximation needs. Seed 7.
# Output: dev/cdfs-out/grad.txt
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
ns <- asNamespace("frmtmb")
set.seed(7)
gerr <- function(g, r) max(abs(g - r) / pmax(1e-3, abs(r)))
check <- function(lab, fad, fref, p0, pts = NULL, third = TRUE) {
  F <- MakeTape(fad, p0)
  g <- F$jacobian(p0)
  r <- numDeriv::jacobian(fref, p0)
  e <- gerr(g, r)
  bad2 <- bad3 <- 0
  if (!is.null(pts)) {
    H <- F$jacfun()
    T3 <- if (third) H$jacfun()
    for (i in seq_len(nrow(pts))) {
      if (!all(is.finite(H$jacobian(pts[i, ])))) bad2 <- bad2 + 1
      if (third && !all(is.finite(T3$jacobian(pts[i, ])))) bad3 <- bad3 + 1
    }
  }
  cat(sprintf("%-44s grad rel err %.1e  nonfinite 2nd %d/%d  3rd %d/%d\n",
              lab, e, bad2, NROW(pts), bad3, NROW(pts)))
}
rpts <- function(n, lo, hi) {
  sapply(seq_along(lo), function(j) runif(n, lo[j], hi[j]))
}
"c" <- RTMB::ADoverload("c")
both <- function(lp) c(lp$lF, lp$lS)
NP <- 40

# student-t: p = (mu, log sigma, log(nu - 1)), q a grid
q <- c(-30, -4, -1, -0.2, 0.3, 1.5, 6, 40)
st <- function(p, f) {
  s <- exp(p[2])
  nu <- 1 + exp(p[3])
  f(q, p[1], s, nu)
}
check("student lF, lS vs stats::pt",
      function(p) st(p, function(q, m, s, nu) {
        both(ns$student_cdf_log((q - m) / s, nu))
      }),
      function(p) st(p, function(q, m, s, nu) {
        c(stats::pt((q - m) / s, nu, log.p = TRUE),
          stats::pt((q - m) / s, nu, lower.tail = FALSE, log.p = TRUE))
      }), c(0.3, 0.2, log(3)), rpts(NP, c(-3, -1, -2), c(3, 1, 4)))
# nu = 4 exactly, the start, where nu / 2 is whole
check("student at nu = 4 (start) vs stats::pt",
      function(p) st(p, function(q, m, s, nu) {
        both(ns$student_cdf_log((q - m) / s, nu))
      }),
      function(p) st(p, function(q, m, s, nu) {
        c(stats::pt((q - m) / s, nu, log.p = TRUE),
          stats::pt((q - m) / s, nu, lower.tail = FALSE, log.p = TRUE))
      }), c(0.3, 0.2, log(3)))

# Gamma: p = (log mu, log shape)
q <- c(0.01, 0.3, 1, 2.5, 6, 20)
check("Gamma lF, lS vs stats::pgamma",
      function(p) both(ns$gamma_cdf_log(q, exp(p[1]), exp(p[2]))),
      function(p) {
        sh <- exp(p[2])
        x <- q * sh / exp(p[1])
        c(stats::pgamma(x, sh, log.p = TRUE),
          stats::pgamma(x, sh, lower.tail = FALSE, log.p = TRUE))
      }, c(log(2), log(3)), rpts(NP, c(-2, -2), c(3, 6)))

# Beta: p = (logit mu, log phi)
q <- c(1e-6, 0.05, 0.3, 0.5, 0.71, 0.95, 1 - 1e-6)
check("Beta lF, lS vs stats::pbeta",
      function(p) {
        m <- plogis(p[1])
        ph <- exp(p[2])
        both(ns$beta_cdf_log(q, m * ph, (1 - m) * ph))
      },
      function(p) {
        m <- plogis(p[1])
        ph <- exp(p[2])
        c(stats::pbeta(q, m * ph, (1 - m) * ph, log.p = TRUE),
          stats::pbeta(q, m * ph, (1 - m) * ph, lower.tail = FALSE,
                       log.p = TRUE))
      }, c(0.4, log(6)), rpts(NP, c(-3, -1), c(3, 7)))

# negative binomial: p = (log mu, log shape), shape 2 exactly whole
q <- c(-1, 0, 1, 3, 9, 30)
nbr <- function(p) {
  c(stats::pnbinom(q, size = exp(p[2]), mu = exp(p[1]), log.p = TRUE),
    stats::pnbinom(q, size = exp(p[2]), mu = exp(p[1]), lower.tail = FALSE,
                   log.p = TRUE))
}
nba <- function(p) both(ns$nbinom_cdf_log(q, p[1], p[2]))
fin <- function(f) function(p) {
  v <- f(p)
  v[!is.finite(v)] <- 0
  v
}
check("negbinomial vs stats::pnbinom", nba, fin(nbr),
      c(log(4), log(1.7)), rpts(NP, c(-3, -2), c(5, 6)))
check("negbinomial at shape = 2 (whole) vs stats", nba, fin(nbr),
      c(log(4), log(2)))

# binomial: p = logit prob, 12 trials
q <- c(-1, 0, 2, 6, 11, 12)
check("binomial vs stats::pbinom",
      function(p) {
        pr <- plogis(p[1])
        both(ns$binomial_cdf_log(q, 12, list(p = pr, q = 1 - pr)))
      },
      fin(function(p) {
        c(stats::pbinom(q, 12, plogis(p[1]), log.p = TRUE),
          stats::pbinom(q, 12, plogis(p[1]), lower.tail = FALSE, log.p = TRUE))
      }), 0.3, rpts(NP, -4, 4))

# poisson: p = log mu
q <- c(-1, 0, 2, 7, 25, 80)
check("poisson vs stats::ppois",
      function(p) both(ns$poisson_cdf_log(q, exp(p[1]))),
      fin(function(p) {
        c(stats::ppois(q, exp(p[1]), log.p = TRUE),
          stats::ppois(q, exp(p[1]), lower.tail = FALSE, log.p = TRUE))
      }), log(5), rpts(NP, -3, 5))

# inverse gaussian: p = (log mu, log shape)
q <- c(0.05, 0.4, 1, 3, 12)
check("inverse.gaussian vs statmod::pinvgauss",
      function(p) both(ns$inverse_gaussian_cdf_log(q, exp(p[1]), exp(p[2]))),
      function(p) {
        c(statmod::pinvgauss(q, exp(p[1]), exp(p[2]), log.p = TRUE),
          statmod::pinvgauss(q, exp(p[1]), exp(p[2]), lower.tail = FALSE,
                             log.p = TRUE))
      }, c(0, log(2)), rpts(NP, c(-1, -2), c(1.5, 5)))

# the rest against their own double evaluation
self <- function(lab, f, p0, lo, hi, third = TRUE) {
  check(lab, f, fin(f), p0, rpts(NP, lo, hi), third)
}
q <- c(-5, -1.3, -0.2, 0, 0.4, 2, 7)
self("skew_normal (mu, log sigma, alpha)",
     function(p) both(ns$skew_normal_cdf_log(q, p[1], exp(p[2]), p[3])),
     c(0.2, 0, 2.5), c(-2, -1, -8), c(2, 1, 8))
self("skew_normal at alpha = 0 (the start's stationary point)",
     function(p) both(ns$skew_normal_cdf_log(q, p[1], exp(p[2]), p[3])),
     c(0.2, 0, 0), c(-2, -1, -1e-3), c(2, 1, 1e-3))
self("exgaussian (mu, log sigma, log beta)",
     function(p) both(ns$exgaussian_cdf_log(q, p[1], exp(p[2]), exp(p[3]))),
     c(0.2, 0, 0.3), c(-2, -2, -2), c(2, 1, 2))
self("asym_laplace (mu, log sigma, logit quantile)",
     function(p) {
       pq <- plogis(p[3])
       both(ns$asym_laplace_cdf_log(q, p[1], exp(p[2]), pq,
                                    list(l = log(pq), l1m = log1p(-pq))))
     }, c(0.1, 0, 0.4), c(-2, -1, -3), c(2, 1, 3))
qv <- c(-pi + 1e-3, -2.5, -1, 0, 0.7, 2, pi - 1e-3)
self("von_mises (2 atan(mu), log kappa)",
     function(p) both(ns$von_mises_cdf_log(qv, 2 * atan(p[1]), exp(p[2]))),
     c(0.3, log(4)), c(-3, -3), c(3, 8))
self("shifted_lognormal (mu, log sigma, log ndt)",
     function(p) {
       both(ns$shifted_lognormal_cdf_log(c(0.6, 1, 2, 5), p[1], exp(p[2]),
                                         exp(p[3])))
     }, c(0, -0.5, log(0.3)), c(-1, -1, -3), c(1, 1, log(0.5)))
self("beta_binomial (logit mu, log phi), 15 trials",
     function(p) {
       m <- plogis(p[1])
       ph <- exp(p[2])
       both(ns$beta_binomial_cdf_log(c(-1, 0, 3, 8, 14, 15), 15, m * ph,
                                     (1 - m) * ph))
     }, c(0.2, log(4)), c(-3, -2), c(3, 6))
self("compois (log mu, log nu)",
     function(p) {
       both(ns$compois_cdf_log(c(-1, 0, 2, 5, 12), exp(p[1]), exp(p[2])))
     }, c(log(3), 0.2), c(-1, -1.5), c(2.5, 1.5))
