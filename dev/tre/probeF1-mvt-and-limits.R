# Probe F. The pieces the registry entry needs, checked before any
# package code is written.
#
# F1. Probe E showed brms puts ONE mixing variable per LEVEL, shared
#     across the coefficients of that level, whether or not the block is
#     correlated. So a d > 1 block is a multivariate t, and a diagonal
#     block is a multivariate t with a diagonal scale matrix, NOT a
#     product of d independent univariate t's. F1 confirms that the
#     closed-form MVT density is the density of brms's construction, by
#     simulation.
# F2. RTMB numerics: does the closed-form MVT log-density tape, does it
#     agree with `dt` at d = 1, and does the gaussian limit hold to the
#     1e-6 the validation gate asks for?
# F3. `simulate()` has to DRAW from the latent. The draw and the density
#     must be the same distribution; checked by a Kolmogorov-Smirnov
#     test on the Mahalanobis form, which is exactly F/d-distributed
#     under the MVT.
#
# Run: Rscript dev/tre/probeF1-mvt-and-limits.R
suppressMessages(library(RTMB))
sink("dev/tre/probeF1.txt")

cat("Probe F: the multivariate-t construction and the numerics\n")
cat("=========================================================\n\n")

# closed-form multivariate-t log density, scale matrix Sigma
ldmvt <- function(x, Sigma, nu) {
  d <- ncol(Sigma)
  L <- chol(Sigma)
  q <- rowSums((x %*% solve(Sigma)) * x)
  lgamma((nu + d) / 2) - lgamma(nu / 2) - d / 2 * log(nu * pi) -
    sum(log(diag(L))) - (nu + d) / 2 * log1p(q / nu)
}

cat("F1. brms's construction against the closed-form MVT.\n")
cat("    b = sqrt(nu * u) * W z, u ~ inv-chi2(nu), z ~ N(0, I),\n")
cat("    W = diag(sd) chol(C)'. Compared by the distribution of the\n")
cat("    Mahalanobis form q = b' Sigma^-1 b, which is d * F(d, nu)\n")
cat("    for a multivariate t.\n\n")
set.seed(1)
cat(sprintf("%-4s %-6s %14s %14s %12s\n", "d", "nu", "KS D", "KS p",
            "mean q sim/thy"))
for (d in c(1, 2, 3)) {
  for (nu in c(3, 5, 10)) {
    sd_v <- c(0.7, 1.3, 0.4)[seq_len(d)]
    C <- diag(d)
    if (d > 1) { C[1, 2] <- C[2, 1] <- 0.5 }
    if (d > 2) { C[1, 3] <- C[3, 1] <- -0.3; C[2, 3] <- C[3, 2] <- 0.2 }
    Sigma <- C * (sd_v %o% sd_v)
    N <- 200000
    u <- 1 / rchisq(N, df = nu)            # inv-chi2(nu), as brms
    z <- matrix(rnorm(N * d), N, d)
    W <- t(chol(Sigma))
    b <- sqrt(nu * u) * (z %*% t(W))
    q <- rowSums((b %*% solve(Sigma)) * b)
    ks <- suppressWarnings(ks.test(q / d, "pf", d, nu))
    cat(sprintf("%-4d %-6s %14.5f %14.4f %6.3f/%5.3f\n", d, nu,
                ks$statistic, ks$p.value, mean(q),
                if (nu > 2) d * nu / (nu - 2) else NA))
  }
}

cat("\n    And the density itself: closed form vs a kernel-free check,\n")
cat("    the log-density ratio between two points, from the mixture\n")
cat("    integral computed numerically.\n\n")
for (d in c(2, 3)) {
  nu <- 4
  sd_v <- c(0.7, 1.3, 0.4)[seq_len(d)]
  C <- diag(d); C[1, 2] <- C[2, 1] <- 0.5
  if (d > 2) { C[1, 3] <- C[3, 1] <- -0.3; C[2, 3] <- C[3, 2] <- 0.2 }
  Sigma <- C * (sd_v %o% sd_v)
  x <- matrix(c(0.3, -0.8, 1.1)[seq_len(d)], 1, d)
  # p(x) = int N(x; 0, nu*u*Sigma) * inv-chi2(u; nu) du
  f <- function(u) {
    sapply(u, function(uu) {
      exp(-0.5 * (rowSums((x %*% solve(Sigma)) * x) / (nu * uu)) -
            d / 2 * log(2 * pi * nu * uu) -
            0.5 * determinant(Sigma)$modulus) *
        dchisq(1 / uu, df = nu) / uu^2
    })
  }
  num <- integrate(f, 0, Inf, rel.tol = 1e-12)$value
  cat(sprintf("d = %d   closed form %14.9f   mixture integral %14.9f",
              d, ldmvt(x, Sigma, nu), log(num)))
  cat(sprintf("   diff %.2e\n", abs(ldmvt(x, Sigma, nu) - log(num))))
}

cat("\n\nF2. RTMB numerics.\n\n")
# the taped form, as the registry entry would write it
ad_ldmvt <- function(B, sdv, Cchol, nu) {
  d <- nrow(B)
  W <- Cchol * 0
  q <- numeric(0)
  NULL
}
tp <- MakeTape(function(p) {
  nu <- p[1]
  sdv <- exp(p[2:3])
  rho <- p[4] / sqrt(1 + p[4]^2)
  C <- RTMB::matrix(c(1, rho, rho, 1), 2, 2)
  Sigma <- C * (RTMB::matrix(sdv, ncol = 1) %*%
                  RTMB::matrix(sdv, nrow = 1))
  B <- RTMB::matrix(c(0.4, -0.2, 1.1, 0.7), 2, 2)
  Si <- RTMB::solve(Sigma)
  q <- colSums(B * (Si %*% B))
  d <- 2
  ld <- lgamma((nu + d) / 2) - lgamma(nu / 2) - d / 2 * log(nu * pi) -
    0.5 * log(Sigma[1, 1] * Sigma[2, 2] - Sigma[1, 2] * Sigma[2, 1]) -
    (nu + d) / 2 * log(1 + q / nu)
  sum(ld)
}, c(4, 0, 0, 0.3))
p0 <- c(4, 0.1, -0.2, 0.3)
cat("taped MVT value          ", sprintf("%.9f", tp(p0)), "\n")
gfd <- sapply(seq_along(p0), function(i) {
  a <- b <- p0; a[i] <- a[i] + 1e-6; b[i] <- b[i] - 1e-6
  (tp(a) - tp(b)) / 2e-6
})
cat("AD gradient              ",
    paste(sprintf("%12.6f", tp$jacobian(p0)), collapse = ""), "\n")
cat("FD gradient              ",
    paste(sprintf("%12.6f", gfd), collapse = ""), "\n")
cat("max diff                 ",
    sprintf("%.2e", max(abs(tp$jacobian(p0) - gfd))), "\n\n")

cat("d = 1: closed-form MVT against RTMB::dt (scaled).\n")
for (nu in c(2.5, 3, 10, 100)) {
  for (s in c(0.3, 1, 4)) {
    x <- 0.77
    a <- ldmvt(matrix(x, 1, 1), matrix(s^2, 1, 1), nu)
    b <- stats::dt(x / s, nu, log = TRUE) - log(s)
    cat(sprintf("  nu %-6s s %-4s  %14.10f  %14.10f  diff %.2e\n",
                nu, s, a, b, abs(a - b)))
  }
}

cat("\nThe gaussian limit, |log dmvt(nu) - log dnorm| at d = 1 and 2.\n\n")
cat(sprintf("%-12s %16s %16s\n", "nu", "d = 1", "d = 2"))
Sig2 <- matrix(c(1, 0.4, 0.4, 2.25), 2, 2)
x2 <- matrix(c(0.6, -1.2), 1, 2)
lnorm2 <- as.vector(mvtnorm::dmvnorm(x2, sigma = Sig2, log = TRUE))
for (nu in c(1e3, 1e4, 1e5, 1e6, 1e7, 1e8)) {
  a1 <- ldmvt(matrix(0.77, 1, 1), matrix(1, 1, 1), nu)
  b1 <- dnorm(0.77, log = TRUE)
  a2 <- ldmvt(x2, Sig2, nu)
  cat(sprintf("%-12.0e %16.3e %16.3e\n", nu, abs(a1 - b1),
              abs(a2 - lnorm2)))
}
cat("\nSame through RTMB's own lgamma on the tape (the value that\n")
cat("would actually be computed):\n\n")
tl <- MakeTape(function(p) {
  nu <- p[1]
  lgamma((nu + 1) / 2) - lgamma(nu / 2) - 0.5 * log(nu * pi) -
    (nu + 1) / 2 * log(1 + 0.77^2 / nu)
}, 1e4)
for (nu in c(1e3, 1e4, 1e5, 1e6, 1e7, 1e8)) {
  cat(sprintf("  nu %-10.0e  |taped - dnorm| = %.3e\n", nu,
              abs(tl(nu) - dnorm(0.77, log = TRUE))))
}

cat("\n\nF3. Drawing from the latent: the simulate() path.\n")
cat("    b = sqrt(nu/chisq(nu)) * chol(Sigma)' z, one chi-square per\n")
cat("    LEVEL. KS of the Mahalanobis form against d*F(d, nu).\n\n")
set.seed(2)
for (d in c(1, 2)) {
  for (nu in c(3, 8)) {
    sd_v <- c(0.9, 1.4)[seq_len(d)]
    C <- diag(d); if (d > 1) C[1, 2] <- C[2, 1] <- -0.35
    Sigma <- C * (sd_v %o% sd_v)
    N <- 100000
    L <- chol(Sigma)                       # upper, as draw_b() uses
    z <- matrix(rnorm(N * d), N, d) %*% L
    b <- z * sqrt(nu / rchisq(N, nu))
    q <- rowSums((b %*% solve(Sigma)) * b)
    ks <- suppressWarnings(ks.test(q / d, "pf", d, nu))
    cat(sprintf("  d %-3d nu %-4s  KS D %8.5f  p %6.4f  ", d, nu,
                ks$statistic, ks$p.value))
    cat(sprintf("var ratio %6.3f (theory %6.3f)\n",
                var(b[, 1]), sd_v[1]^2 * nu / (nu - 2)))
  }
}

sink()
cat("done\n")
