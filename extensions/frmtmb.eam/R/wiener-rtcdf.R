# The distribution of the Wiener response time, for cens() and trunc():
# the survival P(T > t), the distribution function P(T <= t) over both
# boundaries, and the defective distribution function of ONE boundary,
# P(T <= t, boundary b). Everything is a natural LOG, from the first
# term to the last, so that a survival of 1e-500 is a number rather
# than an underflow.
#
# WHY THIS FILE REPLACED THE FIRST VERSION. The first version formed
# the survival's small-time route as 1 - F, floored at 1e-300. Where
# |v| a is large the process has almost surely finished at a small
# normalized time, F rounds to one, and 1 - F is the floor: log -690.8
# for a true log survival of -74 (review of 2026-09-24, 500-bit
# reference, |v| a up to 120). A blend weight of 1e-3 on -690 is enough
# to destroy the answer. Every route here now computes its own quantity
# directly, never as a complement, and each sum is anchored at a term
# known to be of the order of the answer.
#
# THE TWO SURVIVAL ROUTES.
#
# Small time: the killed driftless transition density by the method of
# images, integrated against the Girsanov weight. Each image is a
# normal probability mass on an interval, so
#
#   S = sum_j [ e^{2vja} M(z + 2ja) - e^{2v(ja - z)} M(2ja - z) ],
#   M(c) = P(0 < c + vt + sqrt(t) N < a),  z = a w,
#
# with M taken as a log from whichever normal tail keeps its digits.
# The j = 0 positive term is the anchor: it is the free particle's mass
# inside the interval and is of the order of S wherever this route is
# used.
#
# Large time: the eigenfunction series, with the k = 1 term, which is
# positive, factored out.
#
# THE DEFECTIVE DISTRIBUTION FUNCTION of the lower boundary: the image
# series of Blurton, Kesselmeier and Vollrath (2012) in log form for
# small times, and P_lower - tail_lower for large ones, with the tail
# again an anchored eigenfunction series. The upper boundary is the
# lower one reflected.

# Image counts for the two small-time routes. The blend gives them
# weight up to about nine times their center, and the images decay
# like exp(-(2j - 1)^2 / (2u)), so at u = 3 ten images leave exp(-60).
ddm_rt_ks <- 10L
ddm_rt_kl <- 30L

#' Log of the standard normal distribution function.
#'
#' `RTMB::pnorm(x, log.p = TRUE)` itself. Its value is right far into
#' both tails, but its tape derivative is not: 4.85e8 at x = -2e8 where
#' the truth is 2e8, and NaN at x = 8.3e9. Those arguments arise only at
#' a normalized decision time near zero, so the callers hold that time
#' at `ddm_rt_u_floor` instead of clamping here: at that floor no
#' argument exceeds about 2e6 in size, where the derivative is finite and
#' right to 2e-5, on terms that are exp(-1e12) and contribute nothing
#' (dev/phase3b-nan-hunt4.R and -hunt5.R). Clamping here was tried
#' and measured worse: any clamp spelled without a comparison rounds its
#' argument to the ulp of the clamp's bound, and the image sums below
#' subtract such values, which cost the survival a factor of 10 to 100
#' in accuracy (dev/phase3b-log/cdf-validate3.txt records the run).
#'
#' @noRd
ddm_lpnorm <- function(x) RTMB::pnorm(x, log.p = TRUE)

# The normalized decision time the distribution functions are held at.
# The quantities are 0 or 1 to every bit a double holds there (the
# lower tail's leading term is exp(-w^2 / (2u)) = exp(-5e9 w^2)).
ddm_rt_u_floor <- 1e-10

#' Log of P(lo < N < hi) for a standard normal N, lo < hi.
#'
#' Two spellings of the same number: from the lower tail,
#' log Phi(hi) + log(1 - Phi(lo) / Phi(hi)), which keeps its digits
#' while the interval is left of zero, and from the upper tail, which
#' keeps them while it is right of zero. A tape cannot choose between
#' them on a parameter, so they are blended on the interval's midpoint
#' with a logistic that saturates, to every bit a double holds, before
#' the wrong spelling has lost more than a few digits.
#'
#' @noRd
ddm_lpmass <- function(lo, hi) {
  pl <- ddm_lpnorm
  lhi <- pl(hi)
  llo <- pl(lo)
  L <- lhi + log(ddm_floor(-expm1(llo - lhi), 1e-300))
  ulo <- pl(-lo)
  uhi <- pl(-hi)
  U <- ulo + log(ddm_floor(-expm1(uhi - ulo), 1e-300))
  om <- 0.5 * (1 + ddm_tanh_s(0.5 * (lo + hi) / 0.5))
  (1 - om) * L + om * U
}

#' Log survival, small-time route: images, anchored at the free mass.
#'
#' @noRd
ddm_lsurv_small <- function(t, v, a, w, K = ddm_rt_ks) {
  z <- a * w
  st <- sqrt(t)
  vt <- v * t
  lp <- function(cc) ddm_lpmass((-cc - vt) / st, (a - cc - vt) / st)
  A <- lp(z)
  s <- 1 - exp(ddm_smin(-2 * v * z + lp(-z) - A, ddm_cdf_cap))
  for (j in seq_len(K)) {
    for (jj in c(j, -j)) {
      s <- s + exp(ddm_smin(2 * v * jj * a + lp(z + 2 * jj * a) - A,
                            ddm_cdf_cap)) -
        exp(ddm_smin(2 * v * (jj * a - z) + lp(2 * jj * a - z) - A,
                     ddm_cdf_cap))
    }
  }
  A + log(ddm_floor(s, 1e-300))
}

#' Log survival, large-time route: eigenfunctions, k = 1 factored out.
#'
#' The two boundaries' remaining masses share every eigenvalue, and
#' sin(k pi (1 - w)) is (-1)^(k+1) sin(k pi w), so they are one sum.
#' The larger of the two boundary exponents is factored out as well, so
#' no exponential here can overflow however large |v| a is.
#'
#' @noRd
ddm_lsurv_large <- function(t, v, a, w, K = ddm_rt_kl) {
  va <- v * a
  u <- t / (a * a)
  el <- -va * w
  eu <- va * (1 - w)
  mx <- ddm_floor(el, eu)
  e_l <- exp(el - mx)
  e_u <- exp(eu - mx)
  s1 <- sin(pi * w)
  den1 <- va * va + pi * pi
  s <- 0 * t
  for (k in 1:K) {
    denk <- va * va + k * k * pi * pi
    s <- s + k * (sin(k * pi * w) / s1) * (den1 / denk) *
      exp(-(k * k - 1) * pi * pi * u / 2) * (e_l + (-1)^(k + 1) * e_u)
  }
  log(2 * pi) + log(s1 / den1) - den1 * u / 2 + mx +
    log(ddm_floor(s, 1e-300))
}

#' Log defective distribution function of the LOWER boundary, small time.
#'
#' `ddm_lower_cdf_small()` with every term kept as a log and the sum
#' anchored at the larger of the two j = 0 terms, both of which are
#' positive.
#'
#' @noRd
ddm_llower_small <- function(t, v, a, w, K = ddm_rt_ks) {
  st <- sqrt(t)
  pl <- ddm_lpnorm
  c0 <- a * w
  a0 <- -2 * v * a * w + pl((v * t - c0) / st)
  b0 <- pl(-(v * t + c0) / st)
  A <- ddm_floor(a0, b0)
  s <- exp(a0 - A) + exp(b0 - A)
  for (j in seq_len(K)) {
    cj <- a * (w + 2 * j)
    s <- s + exp(ddm_smin(-2 * v * a * (w + j) + pl((v * t - cj) / st) - A,
                          ddm_cdf_cap)) +
      exp(ddm_smin(2 * v * a * j + pl(-(v * t + cj) / st) - A, ddm_cdf_cap))
    cn <- a * (w - 2 * j)
    s <- s - exp(ddm_smin(-2 * v * a * (w - j) + pl(-(v * t - cn) / st) - A,
                          ddm_cdf_cap)) -
      exp(ddm_smin(-2 * v * a * j + pl((v * t + cn) / st) - A, ddm_cdf_cap))
  }
  A + log(ddm_floor(s, 1e-300))
}

#' Log defective distribution function of the LOWER boundary, large time.
#'
#' `P_lower - tail_lower`, with the tail's k = 1 term factored out and
#' the difference taken as a log.
#'
#' @noRd
ddm_llower_large <- function(t, v, a, w, K = ddm_rt_kl) {
  ltail <- ddm_ltail_lower_large(t, v, a, w, K)
  lP <- ddm_llower_prob_s(v, a, w)
  # ltail cannot exceed lP; outside this route's regime it can, and a
  # -Inf reaching ddm_floor() would be NaN even at a zero blend weight
  lP + log(ddm_floor(-expm1(ddm_smin(ltail - lP, 0)), 1e-300))
}

#' Log of P(T > t, lower boundary): the lower boundary's mass still to
#' come, large-time route, with the k = 1 term factored out.
#'
#' @noRd
ddm_ltail_lower_large <- function(t, v, a, w, K = ddm_rt_kl) {
  va <- v * a
  u <- t / (a * a)
  s1 <- sin(pi * w)
  den1 <- va * va + pi * pi
  s <- 0 * t
  for (k in 1:K) {
    denk <- va * va + k * k * pi * pi
    s <- s + k * (sin(k * pi * w) / s1) * (den1 / denk) *
      exp(-(k * k - 1) * pi * pi * u / 2)
  }
  log(2 * pi) + log(s1 / den1) - den1 * u / 2 - va * w +
    log(ddm_floor(s, 1e-300))
}

#' Log of P(t1 < T <= t2, boundary `up`): an interval-censored row.
#'
#' Three spellings of the same mass, each accurate somewhere, and the
#' one used is chosen by which is accurate rather than by time:
#'
#' - F(t2) (1 - F(t1) / F(t2)), from the two log distribution functions.
#'   Its error is the logs' error over `dE = 1 - F(t1) / F(t2)`, so it
#'   fails when the interval holds a tiny share of what came before it.
#' - T(t1) (1 - T(t2) / T(t1)), from the two log masses still to come,
#'   `T = P - F`, with the same failure over `dL = 1 - T(t2) / T(t1)`
#'   when the interval holds a tiny share of what comes after.
#' - Gauss-Legendre quadrature of the density on `[t1, t2]`, 16 nodes,
#'   which is accurate exactly when both of the above fail: an interval
#'   holding a tiny share of the mass on both sides is short against the
#'   scale on which the density changes.
#'
#' The first version blended the first two on the lower edge's
#' normalized time, and a weight of 1e-6 on a route that had cancelled
#' to the 1e-300 floor put the log mass 25 units out on the review's
#' Rmpfr intervals (review of 2026-09-24, interval.log). A choice on
#' `dL / dE` cannot do that: each route is weighted only where it is the
#' more accurate, and the logistic saturates a factor of about 40 from
#' the crossing.
#'
#' `lF1`, the log defective function at `t1`, may be passed in when the
#' caller already has it.
#'
#' @noRd
ddm_rt_linterval_b <- function(t1, t2, v, a, w, up, lF1 = NULL) {
  vv <- v * (1 - 2 * up)
  ww <- w + up * (1 - 2 * w)
  f1 <- ddm_floor(t1, ddm_rt_u_floor * a * a)
  f2 <- ddm_floor(t2, ddm_rt_u_floor * a * a)
  if (is.null(lF1)) lF1 <- ddm_rt_lcdf_b(t1, v, a, w, up)
  lF2 <- ddm_rt_lcdf_b(t2, v, a, w, up)
  dE <- ddm_floor(-expm1(ddm_smin(lF1 - lF2, 0)), 1e-300)
  early <- lF2 + log(dE)
  # T = P - F, through the large-time tail series wherever its 30 terms
  # have converged (u past 0.008), and through P - F only before that.
  # P - F loses digits as F nears P, and a share of it survived the
  # first version's time blend: on the review's intervals at u = 0.17 and
  # 0.3 it put the log mass 6e-3 and 7e-5 out. Before u = 0.008 the tail
  # is still taken where F is within 1 percent of P, since there P - F
  # has cancelled worse than the unconverged tail errs.
  lP <- ddm_llower_prob_s(vv, a, ww)
  gate <- function(lu, u0) 0.5 * (1 + ddm_tanh_s((lu - log(u0)) / 0.05))
  lTt <- function(ff, lF) {
    lu <- log(ff / (a * a))
    d <- ddm_floor(-expm1(ddm_smin(lF - lP, 0)), 1e-300)
    gd <- 0.5 * (1 + ddm_tanh_s((log(1e-2) - log(d)) / 0.5))
    g <- gate(lu, 0.008)
    wl <- g + (1 - g) * gd * gate(lu, 0.002)
    (1 - wl) * (lP + log(d)) + wl * ddm_ltail_lower_large(ff, vv, a, ww)
  }
  lT1 <- lTt(f1, lF1)
  lT2 <- lTt(f2, lF2)
  dL <- ddm_floor(-expm1(ddm_smin(lT2 - lT1, 0)), 1e-300)
  late <- lT1 + log(dL)
  om <- 0.5 * (1 + ddm_tanh_s((log(dL) - log(dE)) / 0.5))
  diffs <- (1 - om) * early + om * late
  # quadrature, anchored at the largest node so no exponential overflows
  gl <- ddm_gl16
  h <- f2 - f1
  lq <- vector("list", length(gl[["x"]]))
  for (k in seq_along(gl[["x"]])) {
    lq[[k]] <- log(gl[["w"]][[k]]) +
      ddm_lpdf_lower(f1 + h * gl[["x"]][[k]], vv, a, ww)
  }
  m <- lq[[1L]]
  for (k in seq_along(lq)[-1L]) m <- ddm_floor(m, lq[[k]])
  acc <- 0
  for (k in seq_along(lq)) acc <- acc + exp(lq[[k]] - m)
  quad <- log(ddm_floor(h, 1e-300)) + m + log(acc)
  # the quadrature where both differences hold under a thousandth of the
  # mass on their side
  qw <- 0.5 * (1 + ddm_tanh_s((log(1e-3) - log(ddm_floor(dE, dL))) / 0.5))
  (1 - qw) * diffs + qw * quad
}

# Gauss-Legendre nodes on [0, 1] for the interval quadrature, fixed when
# the package is built, so no node count depends on a parameter.
ddm_gl16 <- ddm_gauss_legendre(16L)

#' The blend weight on the large-time route, logistic in log(u).
#'
#' @noRd
ddm_rt_lam <- function(lu, u0) {
  0.5 * (1 + ddm_tanh_s((lu - log(u0)) / ddm_rtcdf_us))
}

#' log(exp(x) + exp(y)), anchored at the larger.
#'
#' @noRd
ddm_lse2 <- function(x, y) {
  m <- ddm_floor(x, y)
  m + log(exp(x - m) + exp(y - m))
}

#' Log distribution function and log survival, both boundaries.
#'
#' `t` is decision time. Returns `list(lF = log P(T <= t),
#' lS = log P(T > t))`. The distribution function's small route is the
#' sum of the two defective ones and its large route is `1 - S` from
#' the large-time survival, which keeps its digits because S is small
#' wherever that route has any weight.
#'
#' `t` is floored at `ddm_rt_u_floor * a^2`: a censoring or truncation bound
#' at or below the non-decision time has F = 0 and S = 1.
#'
#' @noRd
ddm_rt_lcdf2 <- function(t, v, a, w) {
  tt <- ddm_floor(t, ddm_rt_u_floor * a * a)
  lu <- log(ddm_floor(tt / (a * a), ddm_u_floor))
  lamF <- ddm_rt_lam(lu, ddm_rtcdf_u0F)
  lamS <- ddm_rt_lam(lu, ddm_rtcdf_u0S)
  lSl <- ddm_lsurv_large(tt, v, a, w)
  lFs <- ddm_lse2(ddm_llower_small(tt, v, a, w),
                  ddm_llower_small(tt, -v, a, 1 - w))
  lF <- (1 - lamF) * lFs +
    lamF * log(ddm_floor(-expm1(ddm_smin(lSl, 0)), 1e-300))
  lS <- (1 - lamS) * ddm_lsurv_small(tt, v, a, w) + lamS * lSl
  list(lF = lF, lS = lS)
}

#' Log defective distribution function of the boundary `up` names.
#'
#' `up` is 1 for the upper boundary and 0 for the lower, as data. This
#' is what a LEFT- or INTERVAL-censored row reads: such a trial reached
#' a boundary, and which one is known.
#'
#' @noRd
ddm_rt_lcdf_b <- function(t, v, a, w, up) {
  vv <- v * (1 - 2 * up)
  ww <- w + up * (1 - 2 * w)
  tt <- ddm_floor(t, ddm_rt_u_floor * a * a)
  lu <- log(ddm_floor(tt / (a * a), ddm_u_floor))
  lam <- ddm_rt_lam(lu, ddm_rtcdf_u0B)
  (1 - lam) * ddm_llower_small(tt, vv, a, ww) +
    lam * ddm_llower_large(tt, vv, a, ww)
}

# Centers and scale of the blend weights in log(u), one per output,
# because the routes cancel in different places for each. Measured
# against the 700-bit reference of dev/phase3b-cdf-reference-v3.R;
# dev/phase3b-log/cdf-validate3.txt has the sweep.
ddm_rtcdf_u0F <- 0.1
ddm_rtcdf_u0S <- 0.05
ddm_rtcdf_u0B <- 0.2
ddm_rtcdf_us <- 0.12

#' Log of `sinh(x) / x`, with a right second derivative at zero.
#'
#' `ddm_lsinhc()` regularizes `|x|` as `sqrt(x^2 + 1e-20)`, whose second
#' derivative at zero is 1e10. The value is right, but the Laplace
#' approximation differentiates twice, and at a drift of exactly zero,
#' which is every fit's starting value, the inner Hessian came out at
#' -1e15 and the inner Newton step failed (dev/phase3b-hess-check.R).
#' So below |x| = 0.035 this takes the Taylor series,
#' x^2/6 - x^4/180 + x^6/2835 - x^8/37800, exact there to 1e-17, and
#' above 0.07 the exact form, through a logistic in x^2 that saturates
#' to every bit on either side. `ddm_lsinhc()` itself is left as it is,
#' because wiener_gng() uses it and its fits must not move.
#'
#' @noRd
ddm_lsinhc_s <- function(x) {
  x2 <- x * x
  ser <- x2 * (1 / 6 + x2 * (-1 / 180 + x2 * (1 / 2835 - x2 / 37800)))
  om <- 0.5 * (1 + ddm_tanh_s((x2 - 0.0025) / 1e-4))
  (1 - om) * ser + om * ddm_lsinhc(x)
}

#' `ddm_llower_prob()` with `ddm_lsinhc_s()`.
#'
#' @noRd
ddm_llower_prob_s <- function(v, a, w) {
  -v * a * w + log1p(-w) +
    ddm_lsinhc_s(v * a * (1 - w)) - ddm_lsinhc_s(v * a)
}

#' `tanh()` with its argument held to `[-40, 40]`.
#'
#' RTMB's second derivative of `tanh(x)` is NaN once `|x|` passes about
#' 700, where `cosh(x)` overflows, and the blends here reach arguments
#' of 1e8. The Laplace approximation differentiates twice, so the inner
#' Hessian was NaN wherever a blend saturated that far
#' (dev/phase3b-hess-hunt.R). `tanh(40)` is 1 to every bit a double
#' holds, so the clamp changes no value; the clamp is spelled with
#' `abs()`, whose derivatives are 1 or 0.
#'
#' @noRd
ddm_tanh_s <- function(x) {
  y <- -40 + 0.5 * ((x + 40) + abs(x + 40))
  tanh(40 - 0.5 * ((40 - y) + abs(40 - y)))
}

#' Evaluate a distribution-function slot on some rows only.
#'
#' `f(q, dpars, aterms)` runs on the rows `idx` names, with every
#' per-row entry of `dpars` and `aterms` cut to those rows, and the
#' result is scattered into a full-length vector that holds `fill`
#' elsewhere. `idx` comes from the censoring code, which is data, so
#' the selection leaves no branch on the tape. An entry of length one
#' is a constant and is passed as it is.
#'
#' @noRd
ddm_on_rows <- function(f, q, dpars, aterms, idx, fill) {
  "[<-" <- RTMB::ADoverload("[<-")
  n <- length(q)
  out <- rep(fill, n)
  if (!length(idx)) return(out)
  cut <- function(x) if (length(x) == n) x[idx] else x
  out[idx] <- f(q[idx], lapply(dpars, cut), lapply(aterms, cut))
  out
}
