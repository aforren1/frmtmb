# Reviewer, punch round 1, item 1: is lbeta_ad() needed, and is it no
# worse than RTMB::lbeta() in value and in first to third derivatives?
# Reference: log B(a, b) and its partials psi^(n)(a) - psi^(n)(a + b)
# (and psi^(n)(b) - ...) in 256-bit Rmpfr, polygammas by recurrence up
# to 60 plus the asymptotic series.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({ library(RTMB); library(Rmpfr) })
prec <- 256
B2n <- c(1, -1, 1, -1, 5, -691, 7, -3617, 43867, -174611, 854513, -236364091, 8553103, -23749461029)
B2d <- c(6, 30, 42, 30, 66, 2730, 6, 510, 798, 330, 138, 2730, 6, 870)
# psi^(n)(x) in mpfr, n = 0 (digamma), 1, 2
mp_poly <- function(n, x) {
  x <- mpfr(x, prec)
  acc <- mpfr(0, prec)
  while (x < 60) {
    acc <- acc + switch(n + 1, -1 / x, 1 / x^2, -2 / x^3)
    x <- x + 1
  }
  s <- switch(n + 1,
    log(x) - 1 / (2 * x),
    1 / x + 1 / (2 * x^2),
    -1 / x^2 - 1 / x^3)
  for (k in seq_along(B2n)) {
    b <- mpfr(B2n[k], prec) / B2d[k]
    tk <- 2 * k
    s <- s + switch(n + 1,
      -b / (tk * x^tk),
      b / x^(tk + 1),
      -b * (tk + 1) / x^(tk + 2))
  }
  s + acc
}
mp_lbeta <- function(a, b) {
  a <- mpfr(a, prec); b <- mpfr(b, prec)
  lgamma(a) + lgamma(b) - lgamma(a + b)
}
# the reference derivatives of lbeta(a, b) in (a, b)
ref_der <- function(a, b) {
  pa <- function(n) as.numeric(mp_poly(n, a) - mp_poly(n, mpfr(a, prec) + b))
  pb <- function(n) as.numeric(mp_poly(n, b) - mp_poly(n, mpfr(a, prec) + b))
  pab <- function(n) as.numeric(-mp_poly(n, mpfr(a, prec) + b))
  list(v = as.numeric(mp_lbeta(a, b)),
       g = c(pa(0), pb(0)),
       h = c(pa(1), pab(1), pb(1)),
       t = c(pa(2), pab(2), pab(2), pb(2)))   # aaa, aab, abb, bbb
}
# sanity: the mpfr digamma of Rmpfr itself
stopifnot(abs(as.numeric(mp_poly(0, 0.37) - digamma(mpfr(0.37, prec)))) < 1e-35)
lbeta_ad <- frmtmb:::lbeta_ad
mk <- function(f) {
  F <- MakeTape(function(p) f(p[1], p[2]), c(1, 2))
  J <- F$jacfun(); H <- J$jacfun(); T3 <- H$jacfun()
  list(F = F, J = J, H = H, T3 = T3)
}
O <- mk(lbeta_ad)
R <- mk(function(a, b) RTMB::lbeta(a, b))
pick <- function(M, a, b) {
  p <- c(a, b)
  h <- M$H(p); t3 <- M$T3(p)
  list(v = M$F(p), g = as.vector(M$J(p)), h = c(h[1, 1], h[1, 2], h[2, 2]),
       t = c(t3[1, 1], t3[1, 2], t3[2, 2], t3[4, 2]))
}
rel <- function(x, r) max(abs(x - r) / pmax(abs(r), 1e-300))

cat("== 1. the point the worker names: a = 0.01, b = 5e6 ==\n")
for (ab in list(c(0.01, 5e6), c(5e6, 0.01), c(1e-3, 1e7), c(0.3, 8e5))) {
  r <- ref_der(ab[1], ab[2])
  o <- pick(O, ab[1], ab[2]); q <- pick(R, ab[1], ab[2])
  cat(sprintf(paste0("a %g b %g\n  value  rel err  lbeta_ad %.2e  RTMB %.2e\n",
                     "  grad   rel err  lbeta_ad %.2e  RTMB %.2e   (ref %s)\n",
                     "  hess   rel err  lbeta_ad %.2e  RTMB %.2e\n",
                     "  third  rel err  lbeta_ad %.2e  RTMB %.2e\n"),
              ab[1], ab[2], rel(o$v, r$v), rel(q$v, r$v),
              rel(o$g, r$g), rel(q$g, r$g),
              paste(format(r$g, digits = 17), collapse = ", "),
              rel(o$h, r$h), rel(q$h, r$h), rel(o$t, r$t), rel(q$t, r$t)))
}

cat("\n== 2. grid: a, b in 10^(-3..7) by 0.5, plus a == b and a = b (1 +- 1e-8) ==\n")
sh <- 10^seq(-3, 7, by = 0.5)
G <- rbind(expand.grid(a = sh, b = sh), data.frame(a = sh, b = sh),
           data.frame(a = sh, b = sh * (1 + 1e-8)),
           data.frame(a = sh, b = sh * (1 - 0.1)))
res <- t(vapply(seq_len(nrow(G)), function(i) {
  a <- G$a[i]; b <- G$b[i]
  r <- ref_der(a, b); o <- pick(O, a, b); q <- pick(R, a, b)
  c(ov = rel(o$v, r$v), rv = rel(q$v, r$v),
    og = rel(o$g, r$g), rg = rel(q$g, r$g),
    oh = rel(o$h, r$h), rh = rel(q$h, r$h),
    ot = rel(o$t, r$t), rt = rel(q$t, r$t),
    ofin = all(is.finite(c(o$v, o$g, o$h, o$t))),
    rfin = all(is.finite(c(q$v, q$g, q$h, q$t))))
}, numeric(10)))
G <- cbind(G, res)
cat("points", nrow(G), "; all finite: lbeta_ad", sum(G$ofin == 1),
    " RTMB", sum(G$rfin == 1), "\n")
s <- function(v) sprintf("max %.2e  q99 %.2e  median %.2e", max(v),
                         quantile(v, 0.99), median(v))
for (k in c("v", "g", "h", "t")) {
  cat(sprintf("%-5s lbeta_ad %s\n      RTMB     %s\n", k,
              s(G[[paste0("o", k)]]), s(G[[paste0("r", k)]])))
}
cat("\npoints where lbeta_ad is worse than RTMB by more than 10x (and above 1e-12):\n")
for (k in c("v", "g", "h", "t")) {
  o <- G[[paste0("o", k)]]; r <- G[[paste0("r", k)]]
  bad <- o > 10 * r & o > 1e-12
  cat(sprintf("  %s: %d", k, sum(bad)))
  if (any(bad)) {
    w <- which(bad)[order(-o[bad])][1:min(4, sum(bad))]
    cat("  worst:", paste(sprintf("(a %g, b %g: %.1e vs %.1e)", G$a[w], G$b[w],
                                  o[w], r[w]), collapse = " "))
  }
  cat("\n")
}
cat("points where RTMB is worse than lbeta_ad by more than 10x (and above 1e-12):\n")
for (k in c("v", "g", "h", "t")) {
  o <- G[[paste0("o", k)]]; r <- G[[paste0("r", k)]]
  bad <- r > 10 * o & r > 1e-12
  cat(sprintf("  %s: %d", k, sum(bad)))
  if (any(bad)) {
    w <- which(bad)[order(-r[bad])][1:min(4, sum(bad))]
    cat("  worst:", paste(sprintf("(a %g, b %g: %.1e vs %.1e)", G$a[w], G$b[w],
                                  r[w], o[w]), collapse = " "))
  }
  cat("\n")
}
cat("\n== 3. the blend of lbeta_ad: derivative in a at a == b, and across the band ==\n")
for (s0 in c(0.01, 1, 100, 1e5)) {
  for (b in s0 * c(1, 1.1 / 0.9, 1.2, 1.3)) {
    r <- ref_der(s0, b); o <- pick(O, s0, b)
    cat(sprintf("  a %g b %g: grad rel %.1e  hess rel %.1e  third rel %.1e\n",
                s0, b, rel(o$g, r$g), rel(o$h, r$h), rel(o$t, r$t)))
  }
}
