# Tail precision of the RTMB primitives and of frmtmb's log_ibeta_half()
# over the whole unit interval. Seed 3. Output: dev/cdfs-out/probe3.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
lib <- asNamespace("frmtmb")
lih <- lib$log_ibeta_half
rel <- function(a, b) abs(a - b) / pmax(1, abs(b))
cat("== RTMB::ppois left tail, log, against stats::ppois(log.p = TRUE)\n")
for (lam in c(5, 50, 200, 700)) for (k in c(0, 2, 10)) {
  v <- log(RTMB::ppois(k, lam))
  r <- stats::ppois(k, lam, log.p = TRUE)
  cat(sprintf("lambda %4g k %3g  rtmb %.15g  ref %.15g  rel %.2e\n", lam, k,
              v, r, rel(v, r)))
}
cat("== RTMB::pgamma lower tail (poisson upper via P(Y > k) = P(G(k+1) < lam))\n")
for (lam in c(0.5, 5)) for (k in c(10, 50, 200, 700)) {
  v <- log(RTMB::pgamma(lam, shape = k + 1))
  r <- stats::ppois(k, lam, lower.tail = FALSE, log.p = TRUE)
  cat(sprintf("lambda %4g k %3g  rtmb %.15g  ref %.15g  rel %.2e\n", lam, k,
              v, r, rel(v, r)))
}
cat("== RTMB::pgamma lower tail at small x\n")
for (a in c(0.3, 3, 30, 300)) for (x in c(1e-3, 0.1) * a) {
  v <- log(RTMB::pgamma(x, shape = a))
  r <- stats::pgamma(x, a, log.p = TRUE)
  cat(sprintf("a %5g x %8.3g rtmb %.15g ref %.15g rel %.2e\n", a, x, v, r,
              rel(v, r)))
}
cat("== log_ibeta_half over (0, 1): worst relative log error by shape pair\n")
set.seed(3)
for (a in c(0.05, 0.5, 1, 2, 7.3, 40, 300, 5000)) {
  for (b in c(0.05, 0.5, 1, 3.1, 25, 900)) {
    x <- c(1e-12, 1e-6, 1e-3, seq(0.01, 0.99, length.out = 99), 1 - 1e-3,
           1 - 1e-6)
    v <- lih(x, a, b)
    r <- stats::pbeta(x, a, b, log.p = TRUE)
    ok <- is.finite(r) & r > -700
    e <- rel(v[ok], r[ok])
    cat(sprintf("a %7g b %6g  worst %.2e at x = %.6g  nonfinite %d\n", a, b,
                max(e), x[ok][which.max(e)], sum(!is.finite(v[ok]))))
  }
}
