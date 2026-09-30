# Punch 1, B2: the candidate log incomplete beta, defined here so it can
# be swept before it goes into R/families.R. Sourced by the sweep scripts.
# `x - pos(x - m)` is min(x, m) and `x + pos(m - x)` is max(x, m), each
# returning x bit for bit on its own side.
pos <- function(u) 0.5 * (u + abs(u))
cmin <- function(x, m) x - pos(x - m)
cmax <- function(x, m) x + pos(m - x)
# 0 below 0, 1 above 1, and the quintic smoothstep between: value, first
# and second derivatives continuous
smooth01 <- function(t) {
  t <- cmax(cmin(t, 1), 0)
  t * t * t * (10 - 15 * t + 6 * t * t)
}
cand <- function(x, a, b, N = 50L, S0 = 150, S1 = 450, K0 = 3, K1 = 4) {
  cf <- frmtmb:::log_ibeta_cf
  s <- a + b
  m <- (a + 1) / (s + 2)
  sd <- sqrt(a * b / (s * s * (s + 1)))
  # the fraction and its complement, blended over m +- sd
  wd <- smooth01((m + sd - x) / (2 * sd))
  ld <- cf(cmin(x, m + sd), a, b, N)
  lc <- cf(1 - cmax(x, m - sd), b, a, N)
  cap <- log1p(-2^-53)
  # the cap side must come out exactly: cap - pos(cap - lc) is cap there
  lc <- cap - pos(cap - lc)
  lcf <- wd * ld + (1 - wd) * log1p(-exp(lc))
  # near the mean at large shapes, RTMB::pbeta(), at clamped arguments
  # that keep it where its third derivatives are finite when unused
  k <- abs(x - m) / sd
  mn <- cmin(a, b)
  ws <- smooth01((log(mn) - log(S0)) / (log(S1) - log(S0)))
  wk <- smooth01((K1 - k) / (K1 - K0))
  wp <- ws * wk
  a2 <- cmax(a, S0)
  b2 <- cmax(b, S0)
  s2 <- a2 + b2
  m2 <- a2 / s2
  sd2 <- sqrt(a2 * b2 / (s2 * s2 * (s2 + 1)))
  x2 <- cmax(cmin(x, m2 + 4.5 * sd2), m2 - 4.5 * sd2)
  lp <- log(RTMB::pbeta(x2, a2, b2))
  (1 - wp) * lcf + wp * lp
}
