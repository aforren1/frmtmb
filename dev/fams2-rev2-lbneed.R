# Reviewer, item 1: does lbeta_ad() matter to log_ibeta_half()'s
# gradient? The package function against a copy that uses
# RTMB::lbeta(), both against a 300-bit Rmpfr central difference of the
# converged continued fraction, at a small shape beside a large one.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({ library(RTMB); library(Rmpfr) })
ns <- asNamespace("frmtmb")
lib <- ns$log_ibeta_half
# the same function with RTMB::lbeta() swapped in for lbeta_ad()
cf_r <- ns$log_ibeta_cf
body(cf_r) <- do.call(substitute, list(body(cf_r),
                                       list(lbeta_ad = quote(RTMB::lbeta))))
lib_r <- ns$log_ibeta_half
environment(lib_r) <- list2env(list(log_ibeta_cf = cf_r), parent = ns)
prec <- 300
mp_log_ibeta <- function(x, a, b, maxit = 400000) {
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
    a * log(x) + b * log1p(-x) + lgamma(a + b) - lgamma(a) - lgamma(b) -
      log(a) + log(h)
  }
  if (x < (a + 1) / (a + b + 2)) cf(x, a, b) else log1p(-exp(cf(1 - x, b, a)))
}
ref_grad <- function(x, a, b) {
  x <- mpfr(x, prec); a <- mpfr(a, prec); b <- mpfr(b, prec)
  h <- mpfr(1e-30, prec)
  c(as.numeric((mp_log_ibeta(x * (1 + h), a, b) - mp_log_ibeta(x * (1 - h), a, b)) / (2 * h * x)),
    as.numeric((mp_log_ibeta(x, a * (1 + h), b) - mp_log_ibeta(x, a * (1 - h), b)) / (2 * h * a)),
    as.numeric((mp_log_ibeta(x, a, b * (1 + h)) - mp_log_ibeta(x, a, b * (1 - h))) / (2 * h * b)))
}
Tn <- MakeTape(function(p) lib(p[1], p[2], p[3]), c(0.1, 1, 2))$jacfun()
Tr <- MakeTape(function(p) lib_r(p[1], p[2], p[3]), c(0.1, 1, 2))$jacfun()
cat("value check, package vs swapped copy at (0.1, 0.01, 5e6):",
    lib(0.1, 0.01, 5e6), lib_r(0.1, 0.01, 5e6), "\n")
for (ab in list(c(0.01, 5e6), c(5e6, 0.01), c(0.001, 1e7), c(0.3, 8e5), c(2, 3))) {
  a <- ab[1]; b <- ab[2]; m <- (a + 1) / (a + b + 2)
  for (x in unique(c(m * c(0.5, 2, 20), 1e-3, 0.1, 0.4))) {
    if (x >= 0.5) next
    r <- ref_grad(x, a, b)
    gn <- Tn(c(x, a, b)); gr <- Tr(c(x, a, b))
    re <- function(g) abs(g - r) / pmax(abs(r), 1e-300)
    cat(sprintf("a %g b %g x %.3g: d/dx,d/da,d/db rel err  lbeta_ad %s | RTMB::lbeta %s\n",
                a, b, x, paste(sprintf("%.1e", re(gn)), collapse = " "),
                paste(sprintf("%.1e", re(gr)), collapse = " ")))
  }
}
