# Reviewer 2: the Rmpfr series of r2-ref.R, for reuse (copied verbatim).
PREC <- 1200L
mp <- function(x) mpfr(x, PREC)
# log Phi(x) in mpfr, for any x
lPhi <- function(x) log(pnorm(x))

lower_images <- function(t, v, a, w, K = 60L) {
  t <- mp(t); v <- mp(v); a <- mp(a); w <- mp(w)
  z <- a * w
  st <- sqrt(t)
  tot <- mp(0)
  for (k in -K:K) {
    ck <- z + 2 * k * a
    s <- if (k >= 0) 1 else -1   # sign of c_k, since 0 < z < a
    cc <- abs(ck)
    mu <- -v * s
    G <- pnorm((mu * t - cc) / st) + exp(2 * mu * cc) * pnorm((-mu * t - cc) / st)
    tot <- tot + s * exp(v * (ck - z)) * G
  }
  tot
}
p_lower <- function(v, a, w) {
  v <- mp(v); a <- mp(a); w <- mp(w)
  if (v == 0) return(1 - w)
  (exp(-2 * v * a * w) - exp(-2 * v * a)) / (1 - exp(-2 * v * a))
}
tail_eigen <- function(t, v, a, w) {
  u <- mp(t) / mp(a)^2
  va <- mp(v) * mp(a); w <- mp(w)
  # terms fall like exp(-k^2 pi^2 u / 2); stop at 2^-(PREC + 64) of k = 1
  K <- ceiling(sqrt(2 * (PREC + 64) * log(2) / (pi^2 * as.numeric(u)))) + 5
  k <- mp(seq_len(K))
  den <- va^2 + k^2 * Const("pi", PREC)^2
  s <- sum(k * sin(k * Const("pi", PREC) * w) * exp(-den * u / 2) / den)
  2 * Const("pi", PREC) * exp(-va * w) * s
}
