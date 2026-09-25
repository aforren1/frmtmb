# Reviewer 2, item 2: my own high-precision reference for the Wiener
# response-time distribution, and the lane's ddm_rt_lcdf2() /
# ddm_rt_lcdf_b() against it.
#
# Two independent series, both in Rmpfr:
#  (I) images: F_lower(t) = sum_k exp(v (c_k - z)) sgn(c_k) G(|c_k|, mu_k, t),
#      c_k = z + 2 k a, G the first-passage CDF of a drifted Brownian
#      motion to level c (Phi((mu t - c)/sqrt t) + e^{2 mu c}
#      Phi((-mu t - c)/sqrt t)), mu_k = -v sgn(c_k). Derived here from
#      the small-time density by completing the square.
#  (E) eigenfunctions: tail_lower(t) = 2 pi e^{-v a w} sum_k k sin(k pi w)
#      exp(-(v^2 a^2 + k^2 pi^2) u / 2) / (v^2 a^2 + k^2 pi^2),
#      F_lower = P_lower - tail_lower, S = tail_lower + tail_upper.
# The upper boundary is the lower one with v -> -v, w -> 1 - w.
# Where both converge they are compared; the reference for each
# quantity is taken from the series that does not cancel.
source("dev/phase3b-review2/r2-prelude.R")
r2_lib("lane")
suppressPackageStartupMessages({library(Rmpfr); library(frmtmb.eam)})
ns <- asNamespace("frmtmb.eam")
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
ref_point <- function(t, v, a, w) {
  u <- t / a^2
  tl <- tail_eigen(t, v, a, w); tu <- tail_eigen(t, -v, a, 1 - w)
  S_e <- tl + tu
  Fl_e <- p_lower(v, a, w) - tl
  Fu_e <- p_lower(-v, a, 1 - w) - tu
  if (u <= 3) {
    Fl_i <- lower_images(t, v, a, w); Fu_i <- lower_images(t, -v, a, 1 - w)
    S_i <- 1 - Fl_i - Fu_i
    agreeS <- as.numeric(abs(log(S_i) - log(S_e)))
    agreeFl <- as.numeric(abs(log(Fl_i) - log(Fl_e)))
  } else {
    Fl_i <- Fl_e; Fu_i <- Fu_e; agreeS <- NA; agreeFl <- NA
  }
  # take each quantity from the series that does not cancel: images for
  # the defective functions at small u, eigen for S at large u
  Fl <- if (u <= 1) Fl_i else Fl_e
  Fu <- if (u <= 1) Fu_i else Fu_e
  S <- S_e
  c(lS = as.numeric(log(S)), lF = as.numeric(log(Fl + Fu)),
    lFl = as.numeric(log(Fl)), lFu = as.numeric(log(Fu)),
    agreeS = agreeS, agreeFl = agreeFl)
}

set.seed(20260924)
pts <- rbind(
  # targeted: tiny u, bias near the edges, a small and large, v near 0
  expand.grid(u = c(1e-5, 1e-4, 1e-3), w = c(1e-4, 0.001, 0.5, 0.999, 0.9999),
              a = c(0.1, 1.4), v = c(-3, 0, 2)),
  expand.grid(u = c(0.01, 0.3, 2, 8, 40, 200), w = c(0.001, 0.3, 0.9999),
              a = c(0.1, 1.4, 10), v = c(-12, -1e-3, 0, 1e-7, 0.02, 0.04, 0.06,
                                         0.1, 5, 25)),
  # |v| a past the worker's grid: 150 to 400
  expand.grid(u = c(0.003, 0.03, 0.3, 3), w = c(0.2, 0.5, 0.8),
              a = c(8, 20), v = c(-20, 20)),
  data.frame(u = exp(runif(120, log(1e-4), log(50))), w = runif(120, 1e-3, 0.999),
             a = exp(runif(120, log(0.1), log(10))), v = runif(120, -30, 30)))
pts$t <- pts$u * pts$a^2
cat("points:", nrow(pts), "\n")
t0 <- proc.time()[["elapsed"]]
refs <- t(vapply(seq_len(nrow(pts)), function(i) {
  ref_point(pts$t[i], pts$v[i], pts$a[i], pts$w[i])
}, numeric(6)))
cat("reference time", round(proc.time()[["elapsed"]] - t0), "s\n")
res <- cbind(pts, refs)
lane2 <- ns$ddm_rt_lcdf2(res$t, res$v, res$a, res$w)
res$g_lS <- lane2$lS; res$g_lF <- lane2$lF
res$g_lFl <- ns$ddm_rt_lcdf_b(res$t, res$v, res$a, res$w, 0)
res$g_lFu <- ns$ddm_rt_lcdf_b(res$t, res$v, res$a, res$w, 1)
saveRDS(res, "dev/phase3b-review2/ref-points.rds")
