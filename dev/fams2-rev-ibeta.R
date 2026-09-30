# Reviewer, claim 2: log_ibeta_half() accuracy, truncation and
# smoothness, aimed where the worker's random sweep is thin: x near the
# switch point m = (a + 1) / (a + b + 2) at large shapes, x near 0, and
# shapes down to 1e-3. Reference stats::pbeta(log.p = TRUE); where the
# two disagree, a 400-digit Rmpfr evaluation of the same continued
# fraction run to convergence decides.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({ library(RTMB); library(Rmpfr) })
lib <- frmtmb:::log_ibeta_half

# the incomplete beta by Lentz in mpfr, iterated to a relative change of
# 1e-60 (NR betacf), on the converging side, with the complement above m
mp_log_ibeta <- function(x, a, b, prec = 400, maxit = 200000) {
  x <- mpfr(x, prec); a <- mpfr(a, prec); b <- mpfr(b, prec)
  cf <- function(x, a, b) {
    qab <- a + b; qap <- a + 1; qam <- a - 1
    c <- mpfr(1, prec); d <- 1 / (1 - qab * x / qap); h <- d
    for (m in seq_len(maxit)) {
      m2 <- 2 * m
      aa <- m * (b - m) * x / ((qam + m2) * (a + m2))
      d <- 1 / (1 + aa * d); c <- 1 + aa / c; h <- h * d * c
      aa <- -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
      d <- 1 / (1 + aa * d); c <- 1 + aa / c; del <- d * c; h <- h * del
      if (abs(del - 1) < mpfr(1e-60, prec)) break
    }
    lf <- a * log(x) + b * log1p(-x) + lgamma(a + b) - lgamma(a) - lgamma(b) -
      log(a)
    list(v = lf + log(h), it = m)
  }
  if (x < (a + 1) / (a + b + 2)) {
    r <- cf(x, a, b)
    c(as.numeric(r$v), r$it)
  } else {
    r <- cf(1 - x, b, a)
    c(as.numeric(log1p(-exp(r$v))), r$it)
  }
}

err <- function(v, ref) abs(v - ref) / pmax(1, abs(ref))

cat("== A. grid: shapes 1e-3..1e6, x from 1e-300 to 0.49 ==\n")
sh <- 10^seq(-3, 6, by = 1)
xs <- c(1e-300, 1e-100, 1e-20, 1e-8, 1e-4, 1e-2, 0.1, 0.25, 0.4, 0.49)
G <- expand.grid(x = xs, a = sh, b = sh)
G$ref <- suppressWarnings(pbeta(G$x, G$a, G$b, log.p = TRUE))
G$ours <- lib(G$x, G$a, G$b)
G$err <- err(G$ours, G$ref)
ok <- is.finite(G$ref) & G$ref > -700
cat("points", nrow(G), "with -700 < log I:", sum(ok),
    " non-finite ours:", sum(!is.finite(G$ours[ok])), "\n")
cat(sprintf("max err %.2e; count > 1e-10: %d; > 1e-6: %d\n",
            max(G$err[ok]), sum(G$err[ok] > 1e-10), sum(G$err[ok] > 1e-6)))
print(head(G[ok & G$err > 1e-10, ][order(-G$err[ok & G$err > 1e-10]), ], 15))

cat("\n== B. x at the switch point m, shapes equal-ish and large ==\n")
# x = m (1 + delta), both branches' convergence is slowest here
B <- expand.grid(s = 10^(0:6), r = c(0.05, 0.2, 0.45),
                 delta = c(-1e-2, -1e-4, -1e-8, 0, 1e-8, 1e-4, 1e-2))
B$a <- B$r * B$s
B$b <- (1 - B$r) * B$s
B$m <- (B$a + 1) / (B$a + B$b + 2)
B$x <- B$m * (1 + B$delta)
B <- B[B$x < 0.5 & B$x > 0, ]
B$ref <- pbeta(B$x, B$a, B$b, log.p = TRUE)
B$ours <- lib(B$x, B$a, B$b)
B$ours200 <- lib(B$x, B$a, B$b, N = 200L)
B$ours2000 <- lib(B$x, B$a, B$b, N = 2000L)
B$err <- err(B$ours, B$ref)
B$err200 <- err(B$ours200, B$ref)
B$err2000 <- err(B$ours2000, B$ref)
op <- options(width = 200)
print(B[, c("s", "r", "delta", "x", "ref", "err", "err200", "err2000")],
      digits = 3, row.names = FALSE)
cat("\nworst by shape size (err at N = 50):\n")
print(tapply(B$err, B$s, max), digits = 3)

cat("\n== B2. mpfr check at the worst B points ==\n")
wb <- head(B[order(-B$err), ], 6)
for (i in seq_len(nrow(wb))) {
  mp <- mp_log_ibeta(wb$x[i], wb$a[i], wb$b[i])
  cat(sprintf("s %.0e r %.2f delta %+.0e: mpfr %.15g (it %d) pbeta %.15g ours %.15g\n",
              wb$s[i], wb$r[i], wb$delta[i], mp[1], as.integer(mp[2]),
              wb$ref[i], wb$ours[i]))
}

cat("\n== C. xbeta-shaped: q = kappa/(1+2kappa) near mu, phi large ==\n")
# P(Y = 0) = I_q(mu phi, (1 - mu) phi) is at its switch when q ~ mu
C <- expand.grid(kappa = c(0.01, 0.1, 0.5), phi = 10^(1:6),
                 off = c(-0.02, -0.002, 0, 0.002, 0.02))
C$q <- C$kappa / (1 + 2 * C$kappa)
C$mu <- pmin(pmax(C$q + C$off, 1e-3), 0.999)
C$a <- C$mu * C$phi
C$b <- (1 - C$mu) * C$phi
C$ref <- pbeta(C$q, C$a, C$b, log.p = TRUE)
C$ours <- lib(C$q, C$a, C$b)
C$err <- err(C$ours, C$ref)
C$Pzero <- exp(C$ref)
print(C[C$err > 1e-12, c("kappa", "phi", "off", "q", "mu", "Pzero", "err")],
      digits = 3, row.names = FALSE)
cat("count err > 1e-10:", sum(C$err > 1e-10), "of", nrow(C), "\n")
cat("max err by phi:\n")
print(tapply(C$err, C$phi, max), digits = 3)

cat("\n== D. smoothness across the switch (taped, in x) ==\n")
# derivatives in x at m - eps and m + eps against each other and
# against central differences of the reference
tapex <- function(a, b) {
  F <- MakeTape(function(p) lib(p[1], a, b), 0.1)
  J <- F$jacfun(); H <- J$jacfun(); T3 <- H$jacfun()
  list(F = F, J = J, H = H, T3 = T3)
}
refd <- function(x, a, b, k) {
  # derivatives of pbeta(log) in x from the density: d/dx log I = f / I
  f <- function(x) pbeta(x, a, b, log.p = TRUE)
  h <- 1e-4 * x
  switch(k,
         dbeta(x, a, b) / exp(f(x)),
         (f(x + h) - 2 * f(x) + f(x - h)) / h^2,
         (f(x + 2 * h) - 2 * f(x + h) + 2 * f(x - h) - f(x - 2 * h)) /
           (2 * h^3))
}
for (ab in list(c(2, 5), c(0.5, 3), c(30, 70), c(300, 700), c(3000, 7000),
                c(3e4, 7e4))) {
  a <- ab[1]; b <- ab[2]; m <- (a + 1) / (a + b + 2)
  if (m >= 0.5) next
  M <- tapex(a, b)
  for (e in c(1e-3, 1e-7, 1e-12)) {
    lo <- m * (1 - e); hi <- m * (1 + e)
    v <- c(M$F(lo), M$F(hi)); g <- c(M$J(lo), M$J(hi))
    h2 <- c(M$H(lo), M$H(hi)); t3 <- c(M$T3(lo), M$T3(hi))
    cat(sprintf(paste0("a %g b %g eps %.0e: value jump %.2e (ref %.2e) | ",
                       "grad %.10g / %.10g (ref %.10g) | ",
                       "d2 %.6g / %.6g | d3 %.6g / %.6g (fd %.4g)\n"),
                a, b, e, diff(v),
                diff(pbeta(c(lo, hi), a, b, log.p = TRUE)),
                g[1], g[2], refd(m, a, b, 1), h2[1], h2[2], t3[1], t3[2],
                refd(m, a, b, 3)))
  }
}

cat("\n== E. derivative spike from the switch weight near t = 0 ==\n")
# w = (t + |t|) / (2|t| + 1e-300): dw/dt = 2e-300 / (2|t| + 1e-300)^2
a <- 2; b <- 5; m <- (a + 1) / (a + b + 2)
M <- tapex(a, b)
for (t in c(1e-12, 1e-100, 1e-150, 1e-155, 1e-160, 0)) {
  x <- m - t
  cat(sprintf("t = %.0e: x == m? %s  value %.15g  grad %.6g  d2 %.6g  d3 %.6g\n",
              t, x == m, M$F(x), M$J(x), M$H(x), M$T3(x)))
}
cat("(t below ~1e-16 * m rounds x to m itself; the spike needs |t| ~ 1e-150,",
    "which a double x near m = 0.33 cannot express)\n")

cat("\n== F. third derivatives, in (x, log a, log b), finite on a grid ==\n")
F3 <- MakeTape(function(p) lib(p[1], exp(p[2]), exp(p[3])), c(0.1, 0, 0))
T3 <- F3$jacfun()$jacfun()$jacfun()
set.seed(20260930)
nf <- 0; nn <- 0; nfr <- 0
FR <- MakeTape(function(p) log(pbeta(p[1], exp(p[2]), exp(p[3]))), c(0.1, 0, 0))
T3R <- FR$jacfun()$jacfun()$jacfun()
for (i in 1:3000) {
  p <- c(runif(1, 1e-6, 0.4999), runif(1, log(1e-3), log(1e6)),
         runif(1, log(1e-3), log(1e6)))
  lv <- pbeta(p[1], exp(p[2]), exp(p[3]), log.p = TRUE)
  if (!is.finite(lv) || lv < -600) next
  nn <- nn + 1
  if (!all(is.finite(T3(p)))) nf <- nf + 1
  if (!all(is.finite(T3R(p)))) nfr <- nfr + 1
}
cat("points", nn, "; non-finite third derivatives: ours", nf, " RTMB::pbeta", nfr, "\n")
options(op)
