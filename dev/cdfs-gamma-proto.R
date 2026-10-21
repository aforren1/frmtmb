# Prototype: log Q(a, x) by the Legendre continued fraction, evaluated
# bottom up to N terms, and log P(a, x) by the series to N terms. Where
# does each hold, by shape and standardized distance? Output:
# dev/cdfs-out/gamma-proto.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
lq_cf <- function(x, a, N) {
  # Q = e^-x x^a / Gamma(a) / (b1 + a1 / (b2 + a2 / (...))),
  # b_n = x + 2n - 1 - a, a_n = -n (n - a)
  t <- x + 2 * N - 1 - a
  for (n in rev(seq_len(N - 1L))) {
    t <- (x + 2 * n - 1 - a) - n * (n - a) / t
  }
  -x + a * log(x) - lgamma(a) - log(t)
}
lp_ser <- function(x, a, N) {
  # P = e^-x x^a / Gamma(a + 1) * sum_n x^n / ((a + 1) ... (a + n))
  s <- 1
  for (n in rev(seq_len(N))) s <- 1 + s * x / (a + n)
  -x + a * log(x) - lgamma(a + 1) + log(s)
}
rel <- function(v, r) abs(v - r) / pmax(1, abs(r))
cat("== CF for log Q at x = a + d * sqrt(a) + c\n")
for (N in c(30, 60, 100)) {
  for (a in c(0.01, 0.3, 1, 3, 30, 300, 3000, 3e4)) {
    out <- character()
    for (d in c(0, 1, 2, 3, 4, 6, 10, 30)) {
      x <- a + d * sqrt(a) + 1
      e <- rel(lq_cf(x, a, N), pgamma(x, a, lower.tail = FALSE, log.p = TRUE))
      out <- c(out, sprintf("%g:%.0e", d, e))
    }
    cat(sprintf("N %3d a %7g  %s\n", N, a, paste(out, collapse = " ")))
  }
}
cat("== series for log P at x = rho * (a + 1)\n")
for (N in c(30, 60, 100)) {
  for (a in c(0.01, 0.3, 1, 3, 30, 300, 3000, 3e4)) {
    out <- character()
    for (rho in c(1e-6, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9)) {
      x <- rho * (a + 1)
      e <- rel(lp_ser(x, a, N), pgamma(x, a, log.p = TRUE))
      out <- c(out, sprintf("%g:%.0e", rho, e))
    }
    cat(sprintf("N %3d a %7g  %s\n", N, a, paste(out, collapse = " ")))
  }
}
