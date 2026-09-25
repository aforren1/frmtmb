# Punch round 1: a wide high-precision reference for the Wiener
# response-time distribution, over the ranges the review probed:
# v in -20..20, a in 0.3..6, w in 0.001..0.999, u = t / a^2 in 1e-3..10.
#
# Four quantities, all as natural LOGS so that nothing underflows:
#   lS   log P(T > t), both boundaries
#   lF   log P(T <= t)
#   lFl  log P(T <= t, lower boundary), the defective distribution fn
#   lFu  log P(T <= t, upper boundary)
#
# Primary route: the METHOD OF IMAGES, at 700 bits, 60 images each side,
# vectorized over the images.
#   S  = sum_j [ e^{2vja} M(z + 2ja) - e^{2v(ja - z)} M(2ja - z) ]
#        with M(c) = P(0 < c + vt + sqrt(t) N < a), z = a w,
#        which is the killed driftless transition density (images)
#        integrated against the Girsanov weight e^{v(x - z) - v^2 t / 2}.
#   Fl = the inverse-Gaussian image series of Blurton et al. (2012).
#   Fu = Fl at (-v, 1 - w).
# Independent route, where it converges (u >= 0.05): the EIGENFUNCTION
# series at the same precision, with P_b the gambler's-ruin probability:
#   Fl = P_l - 2 pi e^{-vaw} sum_k k sin(k pi w) e^{-lambda_k t} / den_k
#   S  = tail_l + tail_u.
# The two share no algebra. Reported: their worst relative disagreement
# on each quantity, and the identity Fl + Fu + S = 1 on the images.
#
# Run: Rscript dev/phase3b-cdf-reference-v3.R <chunk> <nchunks>
# Output: dev/phase3b-log/cdf-ref3-chunk<k>.csv, merged by
# dev/phase3b-cdf-merge3.R.
.libPaths(c("C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(Rmpfr))
args <- commandArgs(trailingOnly = TRUE)
CH <- as.integer(args[[1]]); NCH <- as.integer(args[[2]])
PREC <- 700L
J <- 60L
mp <- function(x) mpfr(x, PREC)
PI <- Const("pi", PREC)
SQ2 <- sqrt(mp(2))

# P(lo < N < hi) for standard normal N, from the tail that keeps digits
pmass <- function(lo, hi) {
  out <- mp(rep(0, length(lo)))
  up <- asNumeric(lo) > 0
  if (any(up)) {
    out[up] <- (erfc(lo[up] / SQ2) - erfc(hi[up] / SQ2)) / 2
  }
  if (any(!up)) {
    out[!up] <- (erfc(-hi[!up] / SQ2) - erfc(-lo[!up] / SQ2)) / 2
  }
  out
}

img_S <- function(t, v, a, w) {
  z <- a * w
  st <- sqrt(t)
  j <- mp((-J):J)
  M <- function(cc) pmass((-cc - v * t) / st, (a - cc - v * t) / st)
  c1 <- z + 2 * j * a
  c2 <- 2 * j * a - z
  sum(exp(2 * v * j * a) * M(c1) - exp(2 * v * (j * a - z)) * M(c2))
}

pn <- function(x) erfc(-x / SQ2) / 2
img_Fl <- function(t, v, a, w) {
  st <- sqrt(t)
  jp <- mp(0:J)
  jn <- mp((-J):(-1))
  cp <- a * (w + 2 * jp)
  cn <- a * (w + 2 * jn)
  sum(exp(-2 * v * a * (w + jp)) * pn((v * t - cp) / st) +
        exp(2 * v * a * jp) * pn(-(v * t + cp) / st)) -
    sum(exp(-2 * v * a * (w + jn)) * pn(-(v * t - cn) / st) +
          exp(2 * v * a * jn) * pn((v * t + cn) / st))
}

p_lower <- function(v, a, w) {
  if (asNumeric(v) == 0) return(1 - w)
  (exp(-2 * v * a * w) - exp(-2 * v * a)) / (1 - exp(-2 * v * a))
}
eig_tail_l <- function(t, v, a, w) {
  s <- mp(0)
  k <- 1L
  repeat {
    den <- v * v * a * a + k * k * PI * PI
    ex <- -den * t / (2 * a * a)
    term <- (k * sin(k * PI * w) / den) * exp(-v * a * w + ex)
    s <- s + term
    bound <- (k / den) * exp(-v * a * w + ex)
    if (k > 5L && bound < abs(s) * mp(2)^(-PREC)) break
    k <- k + 1L
    if (k > 50000L) stop("eigen series did not converge")
  }
  2 * PI * s
}

gr <- expand.grid(v = c(-20, -12, -5, -1, 0, 1, 5, 12, 20),
                  a = c(0.3, 1, 2, 4, 6),
                  w = c(0.001, 0.05, 0.3, 0.5, 0.7, 0.95, 0.999),
                  u = c(1e-3, 5e-3, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 3, 10))
gr$t <- gr$u * gr$a^2
rows <- which(seq_len(nrow(gr)) %% NCH == CH - 1L)
fmt <- function(x) formatMpfr(x, digits = 25)
res <- vector("list", length(rows))
t0 <- proc.time()[["elapsed"]]
for (ii in seq_along(rows)) {
  i <- rows[ii]
  t <- mp(gr$t[i]); v <- mp(gr$v[i]); a <- mp(gr$a[i]); w <- mp(gr$w[i])
  S <- img_S(t, v, a, w)
  Fl <- img_Fl(t, v, a, w)
  Fu <- img_Fl(t, -v, a, 1 - w)
  ident <- asNumeric(abs(Fl + Fu + S - 1))
  eS <- eFl <- eFu <- NA_real_
  if (gr$u[i] >= 0.05) {
    tl <- eig_tail_l(t, v, a, w)
    tu <- eig_tail_l(t, -v, a, 1 - w)
    Pl <- p_lower(v, a, w)
    Pu <- 1 - Pl
    Se <- tl + tu
    eS <- asNumeric(abs(Se - S) / S)
    eFl <- asNumeric(abs((Pl - tl) - Fl) / Fl)
    eFu <- asNumeric(abs((Pu - tu) - Fu) / Fu)
  }
  res[[ii]] <- data.frame(
    v = gr$v[i], a = gr$a[i], w = gr$w[i], u = gr$u[i], t = gr$t[i],
    lS = fmt(log(S)), lF = fmt(log(Fl + Fu)), lFl = fmt(log(Fl)),
    lFu = fmt(log(Fu)), ident = ident, agreeS = eS, agreeFl = eFl,
    agreeFu = eFu)
}
out <- do.call(rbind, res)
dir.create("dev/phase3b-log", showWarnings = FALSE)
utils::write.csv(out, sprintf("dev/phase3b-log/cdf-ref3-chunk%d.csv", CH),
                 row.names = FALSE)
cat("chunk", CH, "rows", nrow(out), "seconds",
    round(proc.time()[["elapsed"]] - t0), "\n")
