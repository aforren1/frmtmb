# Distribution functions for the core families that had none, so that
# cens() and trunc() reach them. Each family's `lcdf` returns the CDF
# on the probability scale and its `lccdf` the log survivor function
# (frmtmb_family()); both are built here from one log-scale pair per
# family, so the two tails cannot disagree with each other. Everything
# runs on the tape: no branch reads a parameter, and a quantity a
# branch would skip is evaluated at a clamped point and weighted zero.
# The probability `q` a censored row is scored at is DATA in every core
# caller (the recorded response, a censoring or truncation bound), so
# its support edges are resolved off the tape.

#' `log(1 - exp(l))` for `l <= 0`, held finite at `l = 0`.
#'
#' A log-probability that rounds to exactly 0 has a complement below
#' the doubles; the floor returns `log(1e-300)` there instead of `-Inf`,
#' so that a branch weighted zero cannot turn a sum into `NaN`. The floor
#' is added to a positive part rather than taken through `ad_min()`,
#' which forms `l - (l + 1e-300)` above it and returns 0 there, the
#' `log(0)` the floor exists to avoid.
#'
#' Below `l = -log 2` it is `log1p(-exp(l))` instead (Maechler's
#' `log1mexp()` switch): there the answer is near zero, and
#' `log(-expm1(l))` takes the log of a number near one, which keeps
#' its absolute accuracy but not its relative accuracy. That relative
#' accuracy is what an interval's mass in the upper tail is formed from
#' (`log_window()`): `log(-expm1(-36.7))` is `-1.1102e-16` with every
#' digit after the first wrong, and its derivative 10 percent off. Each
#' arm is read at its argument clamped into its own region.
#'
#' @noRd
log1m_exp_ad <- function(l) {
  near0 <- log(-expm1(-(1e-300 + ad_pos(-l - 1e-300))))
  w <- 0.5 * (1 + sign_pos(l + log(2)))
  # the weights, not ad_min(), keep the far arm finite above -log 2:
  # ad_min() of -Inf, a CDF below its support, is Inf - Inf
  far <- log1p(-exp((1 - w) * l - w * log(2)))
  w * near0 + (1 - w) * far
}

#' `sign(x)` with zero sent to `+1`: a two-way switch that has no third
#' state at an exact tie, where `sign()` would weight both arms by one
#' half.
#'
#' @noRd
sign_pos <- function(x) {
  s <- sign(x)
  s + (1 - abs(s))
}

#' Standard normal log CDF and log survivor at `z`, each from
#' `RTMB::pnorm()`'s own log tail, so neither is the complement of the
#' other.
#'
#' @noRd
normal_cdf_log <- function(z) {
  list(lF = RTMB::pnorm(z, log.p = TRUE),
       lS = RTMB::pnorm(z, lower.tail = FALSE, log.p = TRUE))
}

#' Gauss-Legendre nodes and weights on `[0, 1]`, by Golub and Welsch.
#'
#' @noRd
gl_unit <- function(n) {
  k <- seq_len(n - 1L)
  b <- k / sqrt(4 * k * k - 1)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1L)] <- b
  J[cbind(k + 1L, k)] <- b
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(x = (e$values[o] + 1) / 2, w = e$vectors[1L, o]^2)
}

#' Both tails of the regularized incomplete beta: `lF = log I_x(a, b)`
#' and `lS = log(1 - I_x(a, b)) = log I_y(b, a)`, where `y = 1 - x` is
#' passed in because a caller often holds it more exactly than `1 - x`
#' would (`z^2 / (nu + z^2)` beside `nu / (nu + z^2)`).
#'
#' `log_ibeta_half()` is exact to 1.1e-14 relative on `x <= 0.6` for
#' shapes from 0.01 to 1e6 and `x` down to 1e-300 (dev/cdfs-probe4.R).
#' Above 0.6 it stays exact in the tail below the mean (5.0e-16 at
#' `x = 0.8` with shapes (50000, 25), where `log I = -10991`, against
#' 300-bit `Rmpfr::pbetaI()`; dev/cdfs-dbg4.R), but not where `x` is
#' within 1e-5 of 1 and both shapes are small: 0.21 relative at shapes
#' (0.01, 0.01), where `log I = -0.44` (dev/cdfs-dbg5.R). So each tail
#' is taken directly below `x = 0.5` and as the complement of the other
#' tail above 0.6, blended in between, where both direct forms hold. The
#' complement of a log near zero keeps every digit of a small
#' probability down to about 1e-300, where the log itself rounds to
#' zero; below that the direct form is kept on either side, as a second
#' weight in the value of its log.
#'
#' @noRd
log_ibeta_two <- function(x, y, a, b) {
  # the incomplete beta is NaN at exactly one; the largest double below
  # it reads the same limit
  top <- 1 - 2^-53
  d <- log_ibeta_half(ad_min(x, top), a, b)
  e <- log_ibeta_half(ad_min(y, top), b, a)
  wF <- ad_smooth01((x - 0.5) / 0.1) * ad_smooth01((d + 680) / 40)
  wS <- ad_smooth01((y - 0.5) / 0.1) * ad_smooth01((e + 680) / 40)
  list(lF = (1 - wF) * d + wF * log1m_exp_ad(e),
       lS = (1 - wS) * e + wS * log1m_exp_ad(d))
}

#' `log(x^a e^(-x) / Gamma(a))`, the prefactor of both tails of the
#' incomplete gamma.
#'
#' Formed directly its terms are each about `a` in size and cancel to
#' the size of the answer, which costs about `a * 1e-16` of it: 3e-11
#' relative at `a = 3e4` (dev/cdfs-gamma-proto.R). Once `a` passes 12,
#' where `lgamma_binet()` is exact, it is formed around `x = a` instead,
#' as `a (log(x / a) - (x - a) / a) + log(a / (2 pi)) / 2 - binet(a)`,
#' with the log ratio from `log1p()` of the relative distance near the
#' mean, as `log_ibeta_pre()` does for the beta.
#'
#' @noRd
log_igamma_pre <- function(x, a) {
  naive <- a * log(x) - x - lgamma(a)
  a2 <- ad_max(a, 12)
  d <- (x - a2) / a2
  w <- ad_smooth01(4 * abs(d) - 1)
  lr <- (1 - w) * log1p(ad_max(d, -0.75)) + w * (log(x) - log(a2))
  near <- a2 * (lr - d) + 0.5 * log(a2 / (2 * pi)) - lgamma_binet(a2)
  wa <- ad_smooth01((a - 12) / 12)
  (1 - wa) * naive + wa * near
}

#' Both tails of the regularized incomplete gamma: `lP = log P(a, x)`
#' and `lQ = log Q(a, x)`, on the tape, with finite derivatives to the
#' third order.
#'
#' `RTMB::pgamma()` is not used. Its log is `-Inf` once `P < 1e-308`,
#' its upper tail and `log.p` do not tape (they reach
#' `stats::pgamma()` and stop at "Non-numeric argument to mathematical
#' function"), and its third derivatives are not finite at 960 of 2000
#' random points (dev/cdfs-probe1.R, dev/cdfs-probe3.R,
#' dev/cdfs-dbg9.R), which a Laplace approximation needs. Each region
#' has its own form instead:
#'
#' - the power series of `P` to 170 terms, which holds below
#'   `x = 0.8 (a + 1)` for every shape, and up to three standard
#'   deviations above the mode while `a` is below about 185;
#' - the Legendre continued fraction of `Q` to 60 terms, evaluated
#'   bottom up, good to 1e-11 from two standard deviations above
#'   `a + 1` (dev/cdfs-gamma-proto.R);
#' - for `a` above 60, between the two, the tail on the far side of
#'   `x` from the mode by 40-point Gauss-Legendre quadrature of the
#'   integrand divided by its value at `x` (after Numerical Recipes'
#'   `gammpapprox()`, with a window that stops where the integrand has
#'   fallen by `exp(-40)`), good to 1e-14 from `a = 60` to 1e4
#'   (dev/cdfs-gamma-gl.R).
#'
#' Neighboring forms are blended where both hold, each read at its
#' argument clamped into its own region, and the two schemes in the
#' shape from 60 to 100.
#'
#' @noRd
log_igamma_two <- function(x, a, N = 60L, Ns = 170L) {
  ap1 <- a + 1
  sa <- sqrt(a)
  top <- ap1 + 3 * sa
  # the series, to three standard deviations above the mode
  xs <- ad_min(x, top)
  s <- 1
  for (n in rev(seq_len(Ns))) s <- 1 + s * xs / (a + n)
  lser <- log_igamma_pre(xs, a) - log(a) + log(s)
  # the continued fraction, from two standard deviations above
  lo_cf <- ap1 + 2 * sa
  wc <- ad_smooth01((x - lo_cf) / sa)
  xc <- ad_max(x, lo_cf)
  t <- xc + (2 * N - 1) - a
  for (n in rev(seq_len(N - 1L))) t <- (xc + (2 * n - 1) - a) - n * (n - a) / t
  lcf <- log_igamma_pre(xc, a) - log(t)
  # a below 100: the series meets the fraction
  lP_s <- (1 - wc) * lser + wc * log1m_exp_ad(lcf)
  lQ_s <- wc * lcf + (1 - wc) * log1m_exp_ad(lser)
  # a above 60: the series to 0.8 (a + 1), quadrature, the fraction
  gl <- igamma_gl_tails(ad_min(ad_max(x, 0.7 * ap1), top), ad_max(a, 60))
  ws <- ad_smooth01((0.8 - x / ap1) / 0.1)
  lP_l <- ws * lser + (1 - ws) * ((1 - wc) * gl[["lP"]] +
                                    wc * log1m_exp_ad(lcf))
  lQ_l <- wc * lcf + (1 - wc) * (ws * log1m_exp_ad(lser) +
                                   (1 - ws) * gl[["lQ"]])
  wa <- ad_smooth01((a - 60) / 40)
  list(lP = (1 - wa) * lP_s + wa * lP_l, lQ = (1 - wa) * lQ_s + wa * lQ_l)
}

# The 40-point Gauss-Legendre rule igamma_gl_tails() reads
gl40 <- gl_unit(40L)

#' Both tails of the incomplete gamma near the mode, for `a` of 60 and
#' more, by quadrature of the tail on the far side of `x` from the mode
#' `a1 = a - 1`: the lower tail for `x` below it, the upper above, and
#' each other tail as the complement, blended within half a unit of the
#' mode. The integrand `t^a1 e^-t` is taken relative to its value at
#' `x`, so it is at most one over the window and does not underflow far
#' from the mode; the prefactor `a1^a1 e^-a1 / Gamma(a)` is formed from
#' Binet's remainder, without the cancellation of `lgamma(a)`.
#'
#' @noRd
igamma_gl_tails <- function(x, a) {
  a1 <- a - 1
  s <- sqrt(a1)
  lpre <- a1 * log1p(-1 / a) + 1 - 0.5 * log(2 * pi * a) - lgamma_binet(a)
  lf <- function(t) {
    d <- (t - a1) / a1
    a1 * (log1p(d) - d)
  }
  lfx <- lf(x)
  # Each window ends where the integrand, near gaussian at these shapes,
  # has fallen by exp(-50) from its value at x: w sd past x with
  # k w + w^2 / 2 = 50 for x k sd short of the mode on the window's
  # side, so w = sqrt(100 + k^2) - k, which is 10 at the mode, crosses
  # it when k < 0, and tends to 50 / k far off, where the slope sets it.
  # Numerical Recipes' fixed windows lost 6e-7 relative at x one sd
  # below the mode of a = 1e4 (dev/cdfs-acc.R)
  k <- (a1 - x) / s
  lo <- ad_max(x - s * 100 / (sqrt(100 + k * k) + k), 0)
  wup <- s * 100 / (sqrt(100 + k * k) - k)
  # A window that crosses the mode, which happens once x lies more
  # than about 37 sd from it, has an integrand that rises above its
  # value at x; past a = 12,600 the rise at the caller's clamp of
  # x to 0.7 (a + 1) passes exp(709), and the arm the blend weights
  # zero became 0 * Inf = NaN on the whole region below the mode. The
  # cap is active only where the arm is weighted zero: 40 terms of
  # exp(600) times the window stay finite, where exp(700) does not
  # (dev/reviews/2026-10-20-cdfs.md, finding 1).
  sl <- 0
  su <- 0
  for (k in seq_along(gl40$x)) {
    sl <- sl + gl40$w[k] *
      exp(ad_min(lf(lo + (x - lo) * gl40$x[k]) - lfx, 600))
    su <- su + gl40$w[k] *
      exp(ad_min(lf(x + wup * gl40$x[k]) - lfx, 600))
  }
  glP <- lpre + lfx + log((x - lo) * sl)
  glQ <- lpre + lfx + log(wup * su)
  w <- ad_smooth01(a1 + 0.5 - x)
  list(lP = w * glP + (1 - w) * log1m_exp_ad(glQ),
       lQ = (1 - w) * glQ + w * log1m_exp_ad(glP))
}

#' Rows of a discrete CDF outside the support, resolved off the tape.
#'
#' `lp` is the log pair evaluated at a clamped, in-support `q`;
#' `below` and `above` mark rows whose `q` lies under the support
#' (`F = 0`, `S = 1`) or at or over its top (`F = 1`, `S = 0`).
#'
#' @noRd
cdf_edges <- function(lp, below, above = NULL) {
  "[<-" <- RTMB::ADoverload("[<-")
  lF <- lp[["lF"]]
  lS <- lp[["lS"]]
  ib <- which(below)
  if (length(ib)) {
    lF[ib] <- -Inf
    lS[ib] <- 0
  }
  ia <- which(above %||% FALSE)
  if (length(ia)) {
    lF[ia] <- 0
    lS[ia] <- -Inf
  }
  list(lF = lF, lS = lS)
}

#' Broadcast a scalar parameter to one value per row before rows are
#' picked, the way `[<-` on the tape needs it.
#'
#' @noRd
cdf_rows <- function(v, n) if (length(v) == n) v else v + numeric(n)

#' Student-t log CDF and log survivor at the standardized `z`.
#'
#' The tail beyond `|z|` is `I_x(nu / 2, 1 / 2) / 2` with
#' `x = nu / (nu + z^2)`, and the body is one minus that, which is
#' `(1 + I_y(1 / 2, nu / 2)) / 2` with `y = z^2 / (nu + z^2)`; both
#' come from one `log_ibeta_two()`, so the tail is exact far out.
#' `RTMB::pbeta()` is not used: its derivative in a shape is exactly 0
#' at an integer shape (0 against -0.2524 at `pbeta(0.3, 2, 5)`,
#' dev/cdfs-probe2.R), and `nu / 2` is an integer at frmtmb's own start
#' of `nu = 4`.
#'
#' Both forms read `z` through `z^2`, whose derivative is zero at
#' `z = 0`, so near zero the CDF is taken as `1 / 2 + z int_0^1 f(z u)
#' du` by `gl20` instead, with `f` the density: odd in `z`, with the
#' density as its slope at zero. The two are blended for `|z|` from 0.3
#' to 0.5, where `F` is between 0.3 and 0.7 and both hold.
#'
#' @noRd
student_cdf_log <- function(z, nu) {
  z2 <- z * z
  den <- nu + z2
  ib <- log_ibeta_two(nu / den, ad_max(z2 / den, 1e-300), nu / 2, 0.5)
  ltail <- ib[["lF"]] - log(2)
  lbody <- log1p(exp(ib[["lS"]])) - log(2)
  s <- sign(z)
  zn <- ad_max(ad_min(z, 0.5), -0.5)
  m <- 0
  for (k in seq_along(gl20$x)) {
    m <- m + gl20$w[k] * exp(dt_stable(zn * gl20$x[k], nu))
  }
  w <- ad_smooth01((0.5 - abs(z)) / 0.2)
  list(lF = w * log(0.5 + zn * m) +
         (1 - w) * (0.5 * (1 - s) * ltail + 0.5 * (1 + s) * lbody),
       lS = w * log(0.5 - zn * m) +
         (1 - w) * (0.5 * (1 + s) * ltail + 0.5 * (1 - s) * lbody))
}

# The 20-point Gauss-Legendre rule that owen_t_ad() reads
gl20 <- gl_unit(20L)

# The tanh-sinh rule on [0, 1] at step 0.1 to |k h| = 3, whose nodes
# crowd both ends doubly exponentially. Its last nodes round to the end
# points, where the integrand is finite, with weights below 1e-14.
ts_rule <- local({
  k <- seq(-30L, 30L) * 0.1
  x <- tanh(pi / 2 * sinh(k))
  list(u = (x + 1) / 2,
       w = 0.05 * pi / 2 * cosh(k) / cosh(pi / 2 * sinh(k))^2)
})

#' Owen's `T(h, a) = (1 / (2 pi)) int_0^a exp(-h^2 (1 + x^2) / 2) /
#' (1 + x^2) dx` by `gl20`, for `|a| <= 1.1`, where the integrand is
#' smooth on the interval for the small `h` it is used at.
#'
#' @noRd
owen_t_gl <- function(h, a) {
  tot <- 0
  for (k in seq_along(gl20$x)) {
    x <- a * gl20$x[k]
    tot <- tot + gl20$w[k] * exp(-0.5 * h * h * (1 + x * x)) / (1 + x * x)
  }
  a * tot / (2 * pi)
}

#' Owen's `T(h, a)` for every `a`, in absolute accuracy: `owen_t_gl()`
#' for `|a| <= 0.9`, and above 1.1 the identity
#' `T(h, a) = Phi(h) / 2 + Phi(a h) / 2 - Phi(h) Phi(a h) - T(a h, 1 / a)`
#' for `a > 0`, odd in `a`, blended in between where both hold.
#'
#' @noRd
owen_t_ad <- function(h, a) {
  aa <- abs(a)
  wi <- ad_smooth01((aa - 0.9) / 0.2)
  b <- ad_max(aa, 0.9)
  ph <- RTMB::pnorm(h)
  pbh <- RTMB::pnorm(b * h)
  ident <- sign(a) * (0.5 * ph + 0.5 * pbh - ph * pbh - owen_t_gl(b * h, 1 / b))
  (1 - wi) * owen_t_gl(h, ad_max(ad_min(a, 1.1), -1.1)) + wi * ident
}

#' `log G(z, a)` by its angle integral, the peak divided out.
#'
#' The integrand of `(1 / pi) int_atan(a)^(pi / 2) exp(-z^2 / (2
#' cos^2 t)) dt` peaks at `t = atan(max(a, 0))` and falls by `exp(-40)`
#' at `atan(sqrt(a+^2 + 80 / z^2))` above the peak and, when `a < 0`, at
#' `-atan(sqrt(80) / |z|)` below it. Each of those two pieces takes
#' `ts_rule`, after the peak value is taken out on the log scale, so the
#' far tail neither underflows nor loses digits. The ends are written as
#' `pi / 2 - atan(.)` so that `z = 0` does not divide by zero.
#'
#' @noRd
owen_g_angle <- function(z, a) {
  ap <- ad_pos(a)
  az <- abs(z)
  z2 <- z * z
  tp <- atan(ap)
  rgt <- pi / 2 - atan(az / sqrt(ap * ap * z2 + 80))
  lft <- ad_max(atan(a), -(pi / 2 - atan(az / sqrt(80))))
  tot <- 0
  for (k in seq_along(ts_rule$u)) {
    t1 <- tp + (rgt - tp) * ts_rule$u[k]
    t2 <- tp - (tp - lft) * ts_rule$u[k]
    g1 <- exp(-0.5 * z2 * (tan(t1)^2 - ap * ap))
    g2 <- exp(-0.5 * z2 * (tan(t2)^2 - ap * ap))
    tot <- tot + ts_rule$w[k] * ((rgt - tp) * g1 + (tp - lft) * g2)
  }
  -log(pi) - 0.5 * z2 * (1 + ap * ap) + log(tot)
}

#' `log G(z, a)` with `G(z, a) = Phi(-|z|) - 2 T(z, a)`, Owen's `T`: the
#' skew-normal mass of the tail on the side of `z` away from zero.
#'
#' Where the tail is shallow, `t = |z| sqrt(1 + max(a, 0)^2) <= 1.2`,
#' `G` is at least 0.1, so the subtraction through `owen_t_ad()` costs
#' nothing relative. Deeper, the angle integral (`owen_g_angle()`),
#' whose integrand is positive, so the far tail is exact where the
#' subtraction would have cancelled. The two are blended over `t`
#' from 1.2 to 1.6, the shallow form read at `z` held inside its region.
#' Against `stats::integrate()` on the angle at `rel.tol = 1e-13`, and
#' the closed forms at `a = 0` and `a = +-1`, the worst relative error
#' on 281 points (`z` to +-38, `a` to +-1000) is 3.0e-13
#' (dev/cdfs-owen-proto3.R). The angle integral alone with this rule
#' lost 5.4e-11 at `|z| = 0.8`, and a geometric Gauss-Legendre rule
#' 1.4e-3 near `z = 0`, where the drop sits at the end of the interval
#' (dev/cdfs-owen-proto2.R, dev/cdfs-owen-proto.R).
#'
#' @noRd
owen_g_log <- function(z, a) {
  ap <- ad_pos(a)
  sq <- sqrt(1 + ap * ap)
  t <- abs(z) * sq
  ws <- ad_smooth01((1.6 - t) / 0.4)
  zs <- ad_max(ad_min(z, 1.6 / sq), -1.6 / sq)
  shallow <- log(RTMB::pnorm(-abs(zs)) - 2 * owen_t_ad(zs, a))
  ws * shallow + (1 - ws) * owen_g_angle(z, a)
}

#' Skew-normal log CDF and log survivor in brms's parameterization
#' (`mu` the mean, `sigma` the standard deviation, `alpha` the shape).
#'
#' For `z < 0` the CDF is `G(z, alpha)`, and for `z >= 0` the survivor
#' is `G(z, -alpha)` (`owen_g_log()`), so the tail on the side of `z` is
#' always the integral and the other is its complement. Stan's and
#' brms's `Phi(z) - 2 T(z, alpha)` cancels in that tail instead.
#'
#' @noRd
skew_normal_cdf_log <- function(q, mu, sigma, alpha) {
  delta <- alpha / sqrt(1 + alpha * alpha)
  omega <- sigma / sqrt(1 - 2 * delta * delta / pi)
  xi <- mu - omega * delta * sqrt(2 / pi)
  z <- (q - xi) / omega
  s <- sign_pos(z)
  lg <- owen_g_log(z, -s * alpha)
  lc <- log1m_exp_ad(lg)
  list(lF = 0.5 * (1 - s) * lg + 0.5 * (1 + s) * lc,
       lS = 0.5 * (1 + s) * lg + 0.5 * (1 - s) * lc)
}

#' `log M(t)`, Mills's ratio `M(t) = Phi(-t) / phi(t)`, on the tape.
#'
#' Below `t = 6` it is `log Phi(-t) + t^2 / 2 + log(2 pi) / 2`, whose
#' terms cancel by about `t^2 / 2` and cost `1e-16 t^2` of it; above 7
#' it is Laplace's continued fraction `1 / (t + 1 / (t + 2 / (t + ...)))`
#' to 60 terms, evaluated bottom up, whose terms are all positive. The
#' two are blended from 6 to 7, each read at its argument clamped into
#' its own region.
#'
#' @noRd
log_mills <- function(t, N = 60L) {
  td <- ad_min(t, 7)
  direct <- RTMB::pnorm(-td, log.p = TRUE) + 0.5 * td * td +
    0.5 * log(2 * pi)
  tc <- ad_max(t, 6)
  cf <- tc
  for (k in rev(seq_len(N))) cf <- tc + k / cf
  w <- ad_smooth01(t - 6)
  (1 - w) * direct - w * log(cf)
}

#' The ex-gaussian's exponential term `log(exp(u^2 / 2 - u z) Phi(z - u))`,
#' which is the density times `beta`, at the gaussian standardization
#' `z` and `u = sigma / beta`.
#'
#' It equals `log phi(z) + log M(u - z)` with `M` Mills's ratio
#' (`log_mills()`), where the `u^2` terms have cancelled exactly, on
#' paper. Formed as written, `u^2 / 2 - u z` and `log Phi(z - u)` are
#' each about `u^2 / 2` and cancel to `O(log u)`: the value keeps its
#' digits to about `u = 1e3`, but its derivatives were 5.8e-2 wrong
#' there and 2e4 at `u = 2.2e4` (dev/reviews/2026-10-20-cdfs.md,
#' finding 5). `RTMBdist::dexgauss()` switches to a gaussian
#' approximation from `u` of about 20, which is 0.16 wrong in the log
#' density at `u = 20` (RTMBdist-4 in dev/upstream-bugs.md), so the
#' family's density reads this function too.
#'
#' The Mills form cancels the other way, `z^2 / 2` against `t^2 / 2`
#' with `t = u - z`, once `t` is far below zero: the exponential part's
#' long right tail, where `z` is large and `u` small, which put the log
#' survivor 1.5e-2 off at `sigma / beta = 1e-6`. The written form holds
#' there, so each form is used on its side of `t = 0`, blended over `t`
#' from -1 to 1, where both hold.
#'
#' @noRd
exgauss_lx <- function(z, u) {
  t <- u - z
  mills <- -0.5 * z * z - 0.5 * log(2 * pi) + log_mills(t)
  direct <- 0.5 * u * u - u * z + RTMB::pnorm(z - u, log.p = TRUE)
  w <- ad_smooth01((t + 1) / 2)
  w * mills + (1 - w) * direct
}

#' Ex-gaussian log CDF and log survivor in brms's parameterization
#' (`mu` the mean, `beta` the exponential scale, the gaussian at
#' `mu - beta`).
#'
#' With `z` the gaussian standardization and `u = sigma / beta`,
#' `F = Phi(z) - exp(u^2 / 2 - u z) Phi(z - u)` and
#' `S = Phi(-z) + exp(u^2 / 2 - u z) Phi(z - u)`. The survivor is a sum
#' of positive terms, taken in logs; the CDF is the gaussian's times one
#' minus a ratio below one. The exponential factor overflows on its own
#' long before the product does, so it never leaves the log scale.
#'
#' @noRd
exgaussian_cdf_log <- function(q, mu, sigma, beta) {
  u <- sigma / beta
  z <- (q - mu + beta) / sigma
  lx <- exgauss_lx(z, u)
  lphi <- RTMB::pnorm(z, log.p = TRUE)
  # lx - log Phi(z), which F is one minus the exponential of. Where the
  # Mills form holds it is log M(u - z) - log M(-z), log phi(z) having
  # cancelled on paper: as a difference of two logs near -z^2 / 2 it
  # lost 3e-6 of log F deep in the lower tail
  t <- u - z
  w <- ad_smooth01((t + 1) / 2)
  d <- w * (log_mills(t) - log_mills(-z)) + (1 - w) * (lx - lphi)
  list(lF = lphi + log1m_exp_ad(d),
       lS = RTMB::logspace_add(RTMB::pnorm(z, lower.tail = FALSE,
                                           log.p = TRUE), lx))
}

#' Asymmetric Laplace log CDF and log survivor (brms's
#' `asym_laplace_lcdf()` and `_lccdf()`, both arms). `gq` is
#' `dpar_log_complement()` of the quantile, so a quantile near 0 or 1
#' keeps its logs.
#'
#' @noRd
asym_laplace_cdf_log <- function(q, mu, sigma, p, gq) {
  u <- (q - mu) / sigma
  s <- sign_pos(u)
  lo <- gq[["l"]] + (1 - p) * ad_min(u, 0)
  hi <- gq[["l1m"]] - p * ad_max(u, 0)
  list(lF = 0.5 * (1 - s) * lo + 0.5 * (1 + s) * log1m_exp_ad(hi),
       lS = 0.5 * (1 + s) * hi + 0.5 * (1 - s) * log1m_exp_ad(lo))
}

# The tanh-sinh rule at step 0.075 to |k h| = 3, which vm_tail_log()
# reads: at step 0.1 the tail lost 1.0e-9 at u = 0.001, k = 3e4, and at
# 0.075 at most 5.7e-13 (dev/cdfs-vm-proto.R)
ts_rule_vm <- local({
  k <- seq(-40L, 40L) * 0.075
  x <- tanh(pi / 2 * sinh(k))
  list(u = (x + 1) / 2,
       w = 0.0375 * pi / 2 * cosh(k) / cosh(pi / 2 * sinh(k))^2)
})

#' `log P(T > u)` for a centered von Mises `T` and `u` in `[0, pi]`.
#'
#' As the angle integral `exp(-k (1 - cos u)) / (2 pi I0e(k))
#' int_u^hi exp(-2 k sin((s + u) / 2) sin((s - u) / 2)) ds`, whose
#' integrand is at most one, so the far side of a concentrated
#' distribution neither underflows nor cancels; `I0e()` is the
#' exponentially scaled Bessel function the density already reads. The
#' exponent is at least `2 k sin((s - u) / 2)^2`, so it passes 40 within
#' `pi sqrt(20 / k)` of `u`, and the range stops there or at `pi`. A
#' Fourier series (Stan's below `k = 50`) or a normal approximation
#' (Hill 1977, Stan's above) gives the mass on the far side only as a
#' difference of numbers near one; Hill's had 0.59 relative error at
#' `k = 10.5`, a millirad from `-pi` (dev/cdfs-acc.R, first version).
#' Against adaptive quadrature of the density, 168 points with `k` from
#' 1e-8 to 3e4, the worst relative error is 5.7e-13
#' (dev/cdfs-vm-proto.R).
#'
#' @noRd
vm_tail_log <- function(u, k) {
  hi <- pi - ad_pos(pi - u - pi * sqrt(20 / (k + 1e-300)))
  L <- hi - u
  J <- 0
  for (j in seq_along(ts_rule_vm$u)) {
    s <- u + L * ts_rule_vm$u[j]
    J <- J + ts_rule_vm$w[j] *
      exp(-2 * k * sin((s + u) / 2) * sin((s - u) / 2))
  }
  li0 <- log(RTMB::besselI(k, 0, expon.scaled = TRUE))
  -2 * k * sin(u / 2)^2 - log(2 * pi) - li0 + log(1e-300 + L * J)
}

#' The centered von Mises log CDF and log survivor at `t` in `(-pi, pi]`:
#' the tail beyond `|t|` is `vm_tail_log()`, and it is the CDF for
#' `t < 0` and the survivor for `t >= 0`, the other side its complement.
#'
#' @noRd
vm_tails <- function(t, k) {
  tl <- vm_tail_log(abs(t), k)
  tc <- log1m_exp_ad(tl)
  s <- sign_pos(t)
  list(lF = 0.5 * (1 - s) * tl + 0.5 * (1 + s) * tc,
       lS = 0.5 * (1 + s) * tl + 0.5 * (1 - s) * tc)
}

#' `t` moved into `(-pi, pi]` and the number of turns that took, for
#' `|t| < 5 pi`, with `sign()` so that nothing compares on the tape.
#'
#' @noRd
vm_wrap <- function(t) {
  m <- 0
  for (j in 1:2) {
    edge <- (2 * j - 1) * pi
    m <- m + 0.5 * (sign_pos(t - edge) + 1) - 0.5 * (sign_pos(-edge - t) + 1)
  }
  list(t = t - 2 * pi * m, m = m)
}

#' `log P(a < T <= b)` for a centered von Mises `T`, `0 <= b - a <= 2 pi`,
#' from the tails at the two ends after each is wrapped into
#' `(-pi, pi]`. An arc that crosses the antipode is the sum of two
#' tails; one that stays on one side of the mean is the difference of
#' two tails on that side, taken as a ratio; one that holds the mean is
#' one minus two tails. No branch subtracts two numbers near one half,
#' so the far side keeps its digits; but the one-sided branch does
#' subtract two numbers near EACH OTHER when the arc is short.
#'
#' A SHORT arc on one side is still a difference of two tails that
#' nearly agree: the arc `(-pi, -pi + 1e-9]` was good to 4e-8 relative
#' that way, where every other arc held 1.2e-12. So where the log of
#' the density changes by less than about 4 across the arc, which is
#' `V = k L (|sin m| + L)` with `L` the length and `m` the midpoint,
#' the mass is the density integrated over the arc by `gl20` directly,
#' its value at `m` divided out, and the two are blended for `V` from 2
#' to 4. `L` is data in every caller.
#'
#' @noRd
vm_lmass <- function(a, b, k) {
  L <- b - a
  m <- 0.5 * (a + b)
  V <- k * L * (abs(sin(m)) + L)
  wq <- ad_smooth01((4 - V) / 2)
  tot <- 0
  for (j in seq_along(gl20$x)) {
    t <- a + L * gl20$x[j]
    tot <- tot + gl20$w[j] * exp(ad_min(k * (cos(t) - cos(m)), 600))
  }
  li0 <- log(RTMB::besselI(k, 0, expon.scaled = TRUE))
  quad <- log(1e-300 + L * tot) + k * (cos(m) - 1) - log(2 * pi) - li0
  (1 - wq) * vm_lmass_tails(a, b, k) + wq * quad
}

#' `vm_lmass()` from the tails at the two ends.
#'
#' @noRd
vm_lmass_tails <- function(a, b, k) {
  wa <- vm_wrap(a)
  wb <- vm_wrap(b)
  ta <- vm_tails(wa[["t"]], k)
  tb <- vm_tails(wb[["t"]], k)
  d <- wb[["m"]] - wa[["m"]]
  across <- RTMB::logspace_add(ta[["lS"]], tb[["lF"]])
  wpos <- 0.5 * (1 + sign_pos(wa[["t"]]))
  wneg <- (1 - wpos) * 0.5 * (1 + sign_pos(-wb[["t"]]))
  wmid <- 1 - wpos - wneg
  pos <- ta[["lS"]] + log1m_exp_ad(tb[["lS"]] - ta[["lS"]])
  neg <- tb[["lF"]] + log1m_exp_ad(ta[["lF"]] - tb[["lF"]])
  mid <- log1m_exp_ad(RTMB::logspace_add(ta[["lF"]], tb[["lS"]]))
  d * across + (1 - d) * (wpos * pos + wneg * neg + wmid * mid)
}

#' von Mises log CDF and log survivor on `[-pi, pi]`, Stan's and brms's
#' convention: `F(y) = P(-pi <= Y <= y)` for `Y` with mean direction
#' `mu`, which is the centered mass between `-pi - mu` and `y - mu`.
#'
#' @noRd
von_mises_cdf_log <- function(q, mu, k) {
  # the centered angle and kappa ride along for an arc between two
  # recorded angles, which log_window() takes from vm_lmass() itself
  list(lF = vm_lmass(-pi - mu, q - mu, k),
       lS = vm_lmass(q - mu, pi - mu, k),
       ctr = q - mu, kappa = k + 0 * q,
       lmass = function(pl, ph, r) {
         vm_lmass(pl[["ctr"]][r], ph[["ctr"]][r], pl[["kappa"]][r])
       })
}

#' Inverse-gaussian log CDF and log survivor.
#'
#' With `rt = sqrt(shape / q)`,
#' `F = Phi(rt (q / mu - 1)) + e^(2 shape / mu) Phi(-rt (q / mu + 1))`
#' and `S` the same with the second term subtracted. The exponential
#' factor overflows at `shape / mu = 355`
#' while the product is tiny, so both ride in logs; the survivor is the
#' first term times one minus a ratio below one.
#' `RTMBdist::pinvgauss()` forms `F` with the factor outside the log,
#' and its `lower.tail = FALSE` is `log(1 - F)`, which is `-Inf` from
#' `log S = -34`.
#'
#' @noRd
inverse_gaussian_cdf_log <- function(q, mu, shape) {
  r <- sqrt(shape / q)
  a1 <- r * (q / mu - 1)
  l2 <- 2 * shape / mu + RTMB::pnorm(-r * (q / mu + 1), log.p = TRUE)
  l1s <- RTMB::pnorm(a1, lower.tail = FALSE, log.p = TRUE)
  list(lF = RTMB::logspace_add(RTMB::pnorm(a1, log.p = TRUE), l2),
       lS = l1s + log1m_exp_ad(l2 - l1s))
}

#' Negative-binomial log CDF and log survivor at a count `q` (data),
#' from the log mean and log shape. `P(Y <= q) = I_p(shape, q + 1)` with
#' `p = shape / (shape + mu)`, and the survivor is `I_(1-p)(q + 1,
#' shape)`; both probabilities are formed from the two logs, so a mean
#' whose inverse link has underflowed keeps its digits. A count below 0
#' has `F = 0`.
#'
#' negbinomial(), its zero-inflated and hurdle forms, geometric(), and
#' every family built on negbinomial()'s slots read this one function.
#'
#' @noRd
nbinom_cdf_log <- function(q, lmu, lsh) {
  "[<-" <- RTMB::ADoverload("[<-")
  n <- max(length(q), length(lmu), length(lsh))
  q <- rep_len(q, n)
  below <- q < 0
  lmu <- cdf_rows(lmu, n)
  lsh <- cdf_rows(lsh, n)
  # log p and log(1 - p) as log-sums about 0, not as differences from
  # log(shape + mu): at shape 1e8 and mu 1e-8 that sum rounds to
  # log(shape), and p^shape at a count of 0 came out 36 relative off
  lp <- -RTMB::logspace_add(0 * lmu, lmu - lsh)
  l1p <- -RTMB::logspace_add(0 * lmu, lsh - lmu)
  lF <- 0 * lp
  lS <- 0 * lp
  # at a count of 0, which is data, F is p^shape in closed form: every
  # row of a zero-truncated model's normalizer, where the incomplete
  # beta costs about 300 times RTMB::ppois() per row
  # (dev/cdfs-fix1-cost.R)
  i0 <- which(q == 0)
  if (length(i0)) {
    lF[i0] <- exp(lsh[i0]) * lp[i0]
    lS[i0] <- log1m_exp_ad(lF[i0])
  }
  i1 <- which(q >= 1)
  if (length(i1)) {
    lp2 <- log_ibeta_two(exp(lp[i1]), exp(l1p[i1]), exp(lsh[i1]), q[i1] + 1)
    lF[i1] <- lp2[["lF"]]
    lS[i1] <- lp2[["lS"]]
  }
  cdf_edges(list(lF = lF, lS = lS), below)
}

#' Binomial log CDF and log survivor at a count `q` (data) out of
#' `size` trials: `P(Y <= q) = I_(1-p)(size - q, q + 1)` and the
#' survivor `I_p(q + 1, size - q)`. `pq` is `dpar_complement()` of the
#' success probability. Both shapes are data, so the integer-shape
#' derivative of `RTMB::pbeta()` (`student_cdf_log()`) does not arise;
#' `log_ibeta_two()` is used for its tails.
#'
#' @noRd
binomial_cdf_log <- function(q, size, pq) {
  n <- max(length(q), length(size), length(pq[["p"]]))
  q <- rep_len(q, n)
  size <- rep_len(size, n)
  below <- q < 0
  above <- q >= size
  inside <- !below & !above
  qc <- ifelse(inside, q, 0)
  sz <- ifelse(inside, size, 1)
  lp2 <- log_ibeta_two(cdf_rows(pq[["q"]], n), cdf_rows(pq[["p"]], n),
                       sz - qc, qc + 1)
  cdf_edges(lp2, below, above)
}

#' Poisson log CDF and log survivor at a count `q` (data). `RTMB::ppois()`
#' is not used: its upper tail does not tape, and its lower tail is 0
#' once the CDF passes below 1e-308 (`ppois(0, 3000)`), where its log is
#' `-Inf` with a `NaN` derivative.
#'
#' The count is data, so each row takes the cheapest exact form its
#' count allows, chosen off the tape:
#'
#' - `q = 0`: `log F = -mu` and `log S = log(1 - exp(-mu))`, closed form.
#'   This is every row of a zero-truncated model's normalizer.
#' - `1 <= q <= 30`: the CDF as the sum of the masses `0..q` in logs, and
#'   the survivor as the sum of the next 40 masses plus a geometric
#'   remainder at the last term's ratio `mu / (j + 1)`, which bounds the
#'   rest from above because that ratio falls with `j`; the survivor is
#'   the complement of the CDF where that is above about 1e-3, as
#'   `compois_cdf_log()` takes it.
#' - `q > 30`: both tails from one `log_igamma_two()`,
#'   `P(Y <= q) = Q(q + 1, mu)` and `P(Y > q) = P(q + 1, mu)`.
#'
#' The incomplete gamma costs about 170 times `RTMB::ppois()` per row
#' (dev/cdfs-fix1-cost.R): on every row of an importance-sampled
#' zero-truncated fit it took 428 seconds and 27 GB where 0.72.0 took
#' 5 seconds, and with these forms it takes 11 (dev/cdfs-fix1-imp.R).
#'
#' @noRd
poisson_cdf_log <- function(q, mu, qmax = 30L, K = 40L) {
  "[<-" <- RTMB::ADoverload("[<-")
  n <- max(length(q), length(mu))
  q <- rep_len(q, n)
  below <- q < 0
  mu <- cdf_rows(mu, n)
  # every row at 0, as in a zero-truncated model's normalizer: the closed
  # form alone, with no sub-assignment. log(-expm1(-mu)) is the survivor to
  # its absolute accuracy, which is what a normalizer adds; it is the
  # tape of 0.72.0's log(1 - ppois(0, mu)), where log1m_exp_ad()'s two
  # arms doubled the time of an importance-sampled fit
  if (all(q == 0)) return(list(lF = -mu, lS = log(-expm1(-mu))))
  lF <- 0 * mu
  lS <- 0 * mu
  i0 <- which(q == 0)
  if (length(i0)) {
    lF[i0] <- -mu[i0]
    lS[i0] <- log1m_exp_ad(-mu[i0])
  }
  is <- which(q >= 1 & q <= qmax)
  if (length(is)) {
    m <- mu[is]
    lm <- log(m)
    qs <- q[is]
    lpmf <- function(j) j * lm - m - lgamma(j + 1)
    s <- -m
    for (j in seq_len(max(qs))) {
      add <- as.numeric(j <= qs)
      s <- add * RTMB::logspace_add(s, lpmf(pmin(j, qs))) + (1 - add) * s
    }
    up <- lpmf(qs + 1)
    for (k in seq.int(2L, K)) up <- RTMB::logspace_add(up, lpmf(qs + k))
    lr <- ad_min(lm - log(qs + K + 1), -0.01)
    ldir <- RTMB::logspace_add(up, lpmf(qs + K) + lr - log1m_exp_ad(lr))
    lc <- log1m_exp_ad(s)
    w <- ad_smooth01((log(1e-3) - lc) / 2)
    # where the survivor is the small tail, F is its complement too: as a
    # running log sum it keeps only absolute accuracy, which a window in
    # the upper tail reads relatively ((30, 31] at mu 3 was -690.8 for
    # -47.0)
    lF[is] <- w * log1m_exp_ad(ldir) + (1 - w) * s
    lS[is] <- w * ldir + (1 - w) * lc
  }
  il <- which(q > qmax)
  if (length(il)) {
    ig <- log_igamma_two(mu[il], q[il] + 1)
    lF[il] <- ig[["lQ"]]
    lS[il] <- ig[["lP"]]
  }
  cdf_edges(list(lF = lF, lS = lS), below)
}

#' A discrete CDF by summing the log mass function over the support.
#'
#' `lpmf(j)` is the vector of log masses at support point `j` (one
#' value per row; `j` may be a per-row vector). The lower sum runs over
#' `0..q` and the upper over `q + 1..top`, each in log space with one
#' `logspace_add()` per support point, so neither tail is a complement
#' and neither underflows. `top` is the support's last point per row
#' (the trials of a binomial-type count); with `top = NULL` there is no
#' upper sum and the survivor is the complement of the CDF.
#'
#' @noRd
discrete_sum_cdf_log <- function(q, lpmf, top = NULL) {
  "[<-" <- RTMB::ADoverload("[<-")
  n <- length(q)
  below <- q < 0
  above <- if (!is.null(top)) q >= top else rep(FALSE, n)
  inside <- !below & !above
  qc <- ifelse(inside, q, 0)
  jmax <- max(c(qc, if (!is.null(top)) top[inside]), 0)
  lF <- cdf_rows(lpmf(0), n)
  if (jmax >= 1) {
    for (j in seq_len(max(qc))) {
      add <- as.numeric(j <= qc)
      if (any(add > 0)) {
        lF <- add * RTMB::logspace_add(lF, lpmf(pmin(j, qc))) + (1 - add) * lF
      }
    }
  }
  if (is.null(top)) {
    lS <- log1m_exp_ad(lF)
  } else {
    tc <- ifelse(inside, top, 1)
    lS <- cdf_rows(lpmf(pmin(qc + 1, tc)), n)
    for (j in seq_len(jmax)) {
      add <- as.numeric(j > qc + 1 & j <= tc)
      if (any(add > 0)) {
        jj <- pmin(pmax(j, qc + 1), tc)
        lS <- add * RTMB::logspace_add(lS, lpmf(jj)) + (1 - add) * lS
      }
    }
  }
  # a sum of masses whose logs each carry lgamma() rounding can pass
  # zero by that rounding, 3e-10 at shapes of 1e5 over 1000 trials; a
  # log probability is held at or below it
  cdf_edges(list(lF = lF - ad_pos(lF), lS = lS - ad_pos(lS)), below, above)
}

#' Beta-binomial log CDF and log survivor at a count `q` (data) out of
#' `size` trials, with shapes `a` and `b`, by summing the mass function
#' (`discrete_sum_cdf_log()`): the CDF is a generalized hypergeometric
#' function with no closed form, and the support is finite.
#'
#' @noRd
beta_binomial_cdf_log <- function(q, size, a, b) {
  n <- max(length(q), length(size), length(a), length(b))
  q <- rep_len(q, n)
  size <- rep_len(size, n)
  a <- cdf_rows(a, n)
  b <- cdf_rows(b, n)
  lnorm <- lgamma(a + b) - lgamma(a) - lgamma(b) - lgamma(size + a + b) +
    lgamma(size + 1)
  lpmf <- function(j) {
    lnorm + lgamma(j + a) + lgamma(size - j + b) - lgamma(j + 1) -
      lgamma(size - j + 1)
  }
  discrete_sum_cdf_log(q, lpmf, top = size)
}

#' Zero-inflated log CDF and log survivor from the base family's pair
#' `lb` at `q` (data): a point mass `zi` at zero on top of the base,
#' `F(q) = zi 1{q >= 0} + (1 - zi) F_b(q)` and
#' `S(q) = zi 1{q < 0} + (1 - zi) S_b(q)`. `g` is
#' `dpar_log_complement()` of `zi`. `base_zero` marks the rows where the
#' base has no mass at or below `q`, which the caller knows from the
#' data and the base's support; there `F` is the point mass alone,
#' written directly rather than as `logspace_add()` of `-Inf`.
#'
#' The rows are as many as the longest of `q`, the base pair and `g`:
#' a constant `trunc()` bound reaches here as one `q` beside one base
#' value per row, and sizing the rows by `q` gave the point mass to
#' row 1 alone (dev/reviews/2026-10-20-cdfs.md, finding 2).
#'
#' Beside the pair it returns the base pair, the base part's log
#' weight `lw = log(1 - zi)` and the rows with `q >= 0`, from which
#' `log_window()` forms a window that holds none of the point mass.
#'
#' @noRd
zi_cdf_log <- function(q, lb, g, base_zero = q < 0) {
  "[<-" <- RTMB::ADoverload("[<-")
  n <- max(length(q), length(lb[["lF"]]), length(lb[["lS"]]),
           length(g[["l"]]), length(g[["l1m"]]))
  base_zero <- rep_len(base_zero, n)
  q <- rep_len(q, n)
  l <- cdf_rows(g[["l"]], n)
  l1m <- cdf_rows(g[["l1m"]], n)
  lF <- l1m + cdf_rows(lb[["lF"]], n)
  lS <- l1m + cdf_rows(lb[["lS"]], n)
  i0 <- which(q >= 0 & !base_zero)
  if (length(i0)) lF[i0] <- RTMB::logspace_add(l[i0], lF[i0])
  iz <- which(q >= 0 & base_zero)
  if (length(iz)) lF[iz] <- l[iz]
  ineg <- which(q < 0)
  if (length(ineg)) lS[ineg] <- RTMB::logspace_add(l[ineg], lS[ineg])
  list(lF = lF, lS = lS, lw = l1m,
       base = list(lF = cdf_rows(lb[["lF"]], n), lS = cdf_rows(lb[["lS"]], n)),
       in0 = q >= 0)
}

#' Hurdle log CDF and log survivor from the base family's pair `lb` at
#' `q` (data): a mass `hu` at zero and the base truncated to above zero,
#' `S(q) = (1 - hu) S_b(q) / (1 - f_b(0))` for `q >= 0` (`lnz` is
#' `log(1 - f_b(0))`; 0 for a continuous base, which has no mass at
#' zero), and `F(q) = hu + (1 - hu) (1 - S_b(q) / (1 - f_b(0)))`. `g` is
#' `dpar_log_complement()` of `hu`. Every row with `q < 0` has `F = 0`
#' and `S = 1`; a hurdle response is non-negative.
#'
#' `F` is formed directly where `S >= 1/2` and as the complement of `S`
#' above that. As `log1m_exp(lS)` alone it was `hu` to an absolute
#' error of the doubles, `7.8e-6` relative at `hu = 1e-12` and `q = 0`;
#' at `q = 0`, which is data, `F` is `hu` and `S` is `1 - hu`, set
#' directly.
#'
#' It returns the base pair and the base part's log weight
#' `log(1 - hu) - lnz` beside the pair, as `zi_cdf_log()` does.
#'
#' @noRd
hurdle_cdf_log <- function(q, lb, g, lnz = 0) {
  n <- max(length(q), length(lb[["lF"]]), length(lb[["lS"]]),
           length(g[["l"]]), length(g[["l1m"]]), length(lnz))
  q <- rep_len(q, n)
  l1m <- cdf_rows(g[["l1m"]], n)
  lSb <- cdf_rows(lb[["lS"]], n) - lnz
  lS <- l1m + lSb
  direct <- RTMB::logspace_add(cdf_rows(g[["l"]], n),
                               l1m + log1m_exp_ad(lSb))
  w <- 0.5 * (1 + sign_pos(lS - log(0.5)))
  lF <- w * direct + (1 - w) * log1m_exp_ad(lS)
  # at q = 0 the base part holds nothing, so F is hu and S is 1 - hu
  # exactly; through the base's tails it is their difference with lnz,
  # which is zero only to the doubles
  "[<-" <- RTMB::ADoverload("[<-")
  i0 <- which(q == 0)
  if (length(i0)) {
    lF[i0] <- cdf_rows(g[["l"]], n)[i0]
    lS[i0] <- l1m[i0]
  }
  out <- cdf_edges(list(lF = lF, lS = lS), q < 0)
  out$lw <- cdf_rows(l1m - lnz, n)
  out$base <- list(lF = cdf_rows(lb[["lF"]], n), lS = cdf_rows(lb[["lS"]], n))
  out$in0 <- q >= 0
  out
}

#' Beta log CDF and log survivor at `q` (data) with shapes `a` and `b`.
#' A `q` at or below 0 has `F = 0`, and at or above 1 has `F = 1`.
#'
#' @noRd
beta_cdf_log <- function(q, a, b) {
  below <- q <= 0
  above <- q >= 1
  qc <- ifelse(below | above, 0.5, q)
  n <- max(length(q), length(a), length(b))
  lp <- log_ibeta_two(rep_len(qc, n), rep_len(1 - qc, n), cdf_rows(a, n),
                      cdf_rows(b, n))
  cdf_edges(lp, rep_len(below, n), rep_len(above, n))
}

#' Gamma log CDF and log survivor at `q` (data), mean `mu` and shape
#' `shape`: the incomplete gamma at `q shape / mu`. A `q` at or below 0
#' has `F = 0`.
#'
#' @noRd
gamma_cdf_log <- function(q, mu, shape) {
  n <- max(length(q), length(mu), length(shape))
  q <- rep_len(q, n)
  below <- q <= 0
  qc <- ifelse(below, 1, q)
  ig <- log_igamma_two(qc * cdf_rows(shape, n) / cdf_rows(mu, n),
                       cdf_rows(shape, n))
  cdf_edges(list(lF = ig[["lP"]], lS = ig[["lQ"]]), below)
}

#' Lognormal log CDF and log survivor at `q` (data). A `q` at or below
#' 0 has `F = 0`; it is resolved here rather than through `log(0)`,
#' whose derivative in `sigma` is `Inf * 0`.
#'
#' @noRd
lognormal_cdf_log <- function(q, mu, sigma) {
  n <- max(length(q), length(mu), length(sigma))
  q <- rep_len(q, n)
  below <- q <= 0
  z <- (log(ifelse(below, 1, q)) - mu) / sigma
  cdf_edges(list(lF = cdf_rows(RTMB::pnorm(z, log.p = TRUE), n),
                 lS = cdf_rows(RTMB::pnorm(z, lower.tail = FALSE,
                                           log.p = TRUE), n)), below)
}

#' Shifted lognormal log CDF and log survivor: the lognormal at
#' `q - ndt`, where `ndt` is a parameter, so the edge cannot be resolved
#' off the tape. A `q - ndt` at or below 0 is read at 1e-300, where the
#' CDF is `Phi((-690.8 - mu) / sigma)`, zero to the doubles at any
#' `sigma` below about 18, and its derivative in `ndt` is zero.
#'
#' @noRd
shifted_lognormal_cdf_log <- function(q, mu, sigma, ndt) {
  z <- (log(ad_max(q - ndt, 1e-300)) - mu) / sigma
  list(lF = RTMB::pnorm(z, log.p = TRUE),
       lS = RTMB::pnorm(z, lower.tail = FALSE, log.p = TRUE))
}

#' `inverse_gaussian_cdf_log()` with a `q` (data) at or below 0 sent to
#' `F = 0`.
#'
#' @noRd
inverse_gaussian_cdf_edges <- function(q, mu, shape) {
  n <- max(length(q), length(mu), length(shape))
  q <- rep_len(q, n)
  below <- q <= 0
  lp <- inverse_gaussian_cdf_log(ifelse(below, 1, q), cdf_rows(mu, n),
                                 cdf_rows(shape, n))
  cdf_edges(lp, below)
}

#' COM-Poisson log CDF and log survivor at a count `q` (data), in
#' frmtmb's mean parameterization, by summing the mass function
#' (`discrete_sum_cdf_log()`). The mass at 0 and the log rate come from
#' `RTMB::dcompois2()` at 0 and 1, the rest from the recursion
#' `log f(j) = log f(0) + j log(lambda) - nu log(j!)`, so the
#' normalizing constant is formed twice per row rather than once per
#' support point.
#'
#' The support has no top. brms's `com_poisson_lccdf()` is the
#' complement of the CDF, which is `-Inf` once the CDF rounds to one;
#' with the floor of `log1m_exp_ad()` it was `-690.8` with a zero
#' gradient (`-31.29` at `q = 10`, `mu = 0.1`, `nu = 0.3`). So where
#' the complement is below about `1e-3` the survivor is the upper sum
#' itself: `K` terms past `q`, a window whose length is data, then the
#' rest as a geometric series at the ratio `lambda / (j + 1)^nu` of the
#' last term, which bounds it from above because that ratio falls with
#' `j`.
#'
#' @noRd
compois_cdf_log <- function(q, mu, nu, K = 100L) {
  n <- max(length(q), length(mu), length(nu))
  q <- rep_len(q, n)
  mu <- cdf_rows(mu, n)
  nu <- cdf_rows(nu, n)
  lp0 <- RTMB::dcompois2(numeric(n), mu, nu, log = TRUE)
  llam <- RTMB::dcompois2(rep(1, n), mu, nu, log = TRUE) - lp0
  lpmf <- function(j) lp0 + j * llam - nu * lgamma(j + 1)
  lp <- discrete_sum_cdf_log(q, lpmf, top = NULL)
  qc <- pmax(q, -1)
  up <- lpmf(qc + 1)
  for (k in seq.int(2L, K)) up <- RTMB::logspace_add(up, lpmf(qc + k))
  lr <- ad_min(llam - nu * log(qc + K + 1), -0.01)
  ldir <- RTMB::logspace_add(up, lpmf(qc + K) + lr - log1m_exp_ad(lr))
  w <- ad_smooth01((log(1e-3) - lp[["lS"]]) / 2)
  # F mirrors the survivor: the summed CDF is accurate only in absolute
  # terms, and an upper-tail window needs it relatively
  list(lF = w * log1m_exp_ad(ldir) + (1 - w) * lp[["lF"]],
       lS = w * ldir + (1 - w) * lp[["lS"]])
}

# --- The per-family pairs the family objects' slots read ---------------
#
# Each takes the slot's own arguments, plus the links a robust density
# input needs, and returns list(lF, lS), which is the family's
# `cdf_log` slot; frmtmb_family() derives `lcdf` and `lccdf` from it.

#' student()'s pair, at the scale a known `se()` sets.
#'
#' @noRd
student_cdf_pair <- function(q, dpars, aterms) {
  s <- resid_sd(dpars[["sigma"]], aterms)
  student_cdf_log((q - dpars[["mu"]]) / s, dpars[["nu"]])
}

#' The log mean and log shape a negative binomial CDF reads, with
#' `rate()` folded into both as brms folds it (`mu * d`, `shape * d`).
#' `lk_shape = NULL` is a shape held at 1.
#'
#' @noRd
nbinom_log_pars <- function(dpars, aterms, lk, lk_shape) {
  lmu <- robust_logmu(dpars, lk) %||% log(dpars[["mu"]])
  lsh <- if (is.null(lk_shape)) 0 else dpar_log(dpars, "shape", lk_shape)
  d <- aterms[["rate"]]
  if (!is.null(d)) {
    lmu <- lmu + log(d)
    lsh <- lsh + log(d)
  }
  list(lmu = lmu, lsh = lsh)
}

#' negbinomial()'s pair; `lk_shape = NULL` is geometric().
#'
#' @noRd
negbinomial_cdf_pair <- function(q, dpars, aterms, lk, lk_shape) {
  lp <- nbinom_log_pars(dpars, aterms, lk, lk_shape)
  nbinom_cdf_log(q, lp[["lmu"]], lp[["lsh"]])
}

#' zero_inflated_negbinomial()'s pair.
#'
#' @noRd
zi_negbinomial_cdf_pair <- function(q, dpars, aterms, lk, lk_shape, lk_zi) {
  base <- negbinomial_cdf_pair(q, dpars, aterms, lk, lk_shape)
  zi_cdf_log(q, base, dpar_log_complement(dpars, "zi", lk_zi))
}

#' hurdle_negbinomial()'s pair. The normalizer is `log(1 - P(0))` with
#' `log P(0) = -shape log(1 + mu / shape)`, the ratio in logs, as the
#' density forms it.
#'
#' @noRd
hurdle_negbinomial_cdf_pair <- function(q, dpars, aterms, lk, lk_shape,
                                        lk_hu) {
  lp <- nbinom_log_pars(dpars, aterms, lk, lk_shape)
  base <- nbinom_cdf_log(q, lp[["lmu"]], lp[["lsh"]])
  lp0 <- -exp(lp[["lsh"]]) *
    RTMB::logspace_add(0 * lp[["lmu"]], lp[["lmu"]] - lp[["lsh"]])
  hurdle_cdf_log(q, base, dpar_log_complement(dpars, "hu", lk_hu),
                 RTMB::logspace_sub(0 * lp0, lp0))
}

#' poisson()'s pair, at the mean `rate()` sets.
#'
#' @noRd
poisson_cdf_pair <- function(q, dpars, aterms, lk) {
  poisson_cdf_log(q, rate_mu(dpars, lk, aterms))
}

#' zero_inflated_poisson()'s pair.
#'
#' @noRd
zi_poisson_cdf_pair <- function(q, dpars, aterms, lk_zi) {
  base <- poisson_cdf_log(q, dpars[["mu"]])
  zi_cdf_log(q, base, dpar_log_complement(dpars, "zi", lk_zi))
}

#' hurdle_poisson()'s pair. The positive part's normalizer is
#' `log(1 - exp(-mu))`, through `expm1()` as the density takes it.
#'
#' @noRd
hurdle_poisson_cdf_pair <- function(q, dpars, aterms, lk_hu) {
  mu <- dpars[["mu"]]
  base <- poisson_cdf_log(q, mu)
  hurdle_cdf_log(q, base, dpar_log_complement(dpars, "hu", lk_hu),
                 log(-expm1(-mu)))
}

#' binomial()'s pair; `name` is the family the trials refusal names.
#'
#' @noRd
binomial_cdf_pair <- function(q, dpars, aterms, lk, name) {
  binomial_cdf_log(q, aterm_trials(aterms, name),
                   dpar_complement(dpars, "mu", lk))
}

#' zero_inflated_binomial()'s pair.
#'
#' @noRd
zi_binomial_cdf_pair <- function(q, dpars, aterms, lk, lk_zi) {
  base <- binomial_cdf_pair(q, dpars, aterms, lk, "zero_inflated_binomial")
  zi_cdf_log(q, base, dpar_log_complement(dpars, "zi", lk_zi))
}

#' beta_binomial()'s pair, with the shapes off the log-odds as the
#' density takes them.
#'
#' @noRd
beta_binomial_cdf_pair <- function(q, dpars, aterms, lk, name) {
  mp <- dpar_complement(dpars, "mu", lk)
  beta_binomial_cdf_log(q, aterm_trials(aterms, name),
                        mp[["p"]] * dpars[["phi"]], mp[["q"]] * dpars[["phi"]])
}

#' zero_inflated_beta_binomial()'s pair.
#'
#' @noRd
zi_beta_binomial_cdf_pair <- function(q, dpars, aterms, lk, lk_zi) {
  base <- beta_binomial_cdf_pair(q, dpars, aterms, lk,
                                 "zero_inflated_beta_binomial")
  zi_cdf_log(q, base, dpar_log_complement(dpars, "zi", lk_zi))
}

#' Beta()'s pair, with the shapes off the log-odds.
#'
#' @noRd
beta_cdf_pair <- function(q, dpars, aterms, lk) {
  mp <- dpar_complement(dpars, "mu", lk)
  beta_cdf_log(q, mp[["p"]] * dpars[["phi"]], mp[["q"]] * dpars[["phi"]])
}

#' zero_inflated_beta()'s pair. The point mass is at 0, where the beta
#' has none, so a row at `q = 0` has the mass alone.
#'
#' @noRd
zi_beta_cdf_pair <- function(q, dpars, aterms, lk, lk_zi) {
  zi_cdf_log(q, beta_cdf_pair(q, dpars, aterms, lk),
             dpar_log_complement(dpars, "zi", lk_zi), base_zero = q <= 0)
}

#' hurdle_gamma()'s pair.
#'
#' @noRd
hurdle_gamma_cdf_pair <- function(q, dpars, aterms, lk_hu) {
  base <- gamma_cdf_log(q, dpars[["mu"]], dpars[["shape"]])
  hurdle_cdf_log(q, base, dpar_log_complement(dpars, "hu", lk_hu))
}

#' hurdle_lognormal()'s pair.
#'
#' @noRd
hurdle_lognormal_cdf_pair <- function(q, dpars, aterms, lk_hu) {
  base <- lognormal_cdf_log(q, dpars[["mu"]], dpars[["sigma"]])
  hurdle_cdf_log(q, base, dpar_log_complement(dpars, "hu", lk_hu))
}

#' asym_laplace()'s pair.
#'
#' @noRd
asym_laplace_cdf_pair <- function(q, dpars, aterms, lk_quantile) {
  asym_laplace_cdf_log(q, dpars[["mu"]], dpars[["sigma"]],
                       dpars[["quantile"]],
                       dpar_log_complement(dpars, "quantile", lk_quantile))
}

#' zero_inflated_asym_laplace()'s pair. brms's
#' `zero_inflated_asym_laplace_lccdf()` adds the point mass at `y <= 0`,
#' so at `y = 0` it is `P(Y >= 0)` where its `lcdf` is `P(Y <= 0)`, and
#' the two sum to `1 + zi`. Here the survivor is `P(Y > y)` at every
#' `y`, the convention every other family's `lccdf` follows.
#'
#' @noRd
zi_asym_laplace_cdf_pair <- function(q, dpars, aterms, lk_quantile, lk_zi) {
  base <- asym_laplace_cdf_pair(q, dpars, aterms, lk_quantile)
  zi_cdf_log(q, base, dpar_log_complement(dpars, "zi", lk_zi),
             base_zero = rep(FALSE, length(q)))
}

# --- Truncated means ---------------------------------------------------
#
# fitted() of a trunc() response is E[Y | lb <= Y <= ub]. These are the
# families brms 2.23.0 has a posterior_epred_trunc_*() for; brms refuses
# the rest, and so does response_mean(). They run off the tape.

#' Kim (2008), as brms's `posterior_epred_trunc_student()`, with the
#' gamma ratio and the powers in logs so that a large `nu` does not
#' overflow `nu^(nu / 2)`.
#'
#' @noRd
student_trunc_mean <- function(mu, s, nu, lb, ub) {
  zl <- (lb - mu) / s
  zu <- (ub - mu) / s
  # the window's mass from the tail it lies in, as pnorm_diff() does
  P <- ifelse(zl > 0, stats::pt(-zl, nu) - stats::pt(-zu, nu),
              stats::pt(zu, nu) - stats::pt(zl, nu))
  lg <- lgamma((nu - 1) / 2) - lgamma(nu / 2) - lgamma(0.5) - log(2 * P)
  pw <- function(z) exp(lg + nu / 2 * log(nu) - (nu - 1) / 2 * log(nu + z^2))
  mu + s * (pw(zl) - pw(zu))
}

#' brms's `posterior_epred_trunc_gamma()`: the partial first moment is
#' `mu` times the incomplete gamma with the shape raised by one.
#'
#' @noRd
gamma_trunc_mean <- function(mu, shape, lb, ub) {
  lo <- pmax(lb, 0) * shape / mu
  hi <- ub * shape / mu
  mu * (stats::pgamma(hi, shape + 1) - stats::pgamma(lo, shape + 1)) /
    (stats::pgamma(hi, shape) - stats::pgamma(lo, shape))
}

#' `E[Y | lb <= Y <= ub]` of a binomial count, the lower bound INCLUSIVE
#' as the likelihood reads it: `y f(y; n, p) = n p f(y - 1; n - 1, p)`,
#' so the partial sum is `n p (F_(n-1)(ub - 1) - F_(n-1)(lb - 2))`.
#'
#' @noRd
binomial_trunc_mean <- function(p, size, lb, ub) {
  num <- size * p * (stats::pbinom(ub - 1, size - 1, p) -
                       stats::pbinom(lb - 2, size - 1, p))
  num / (stats::pbinom(ub, size, p) - stats::pbinom(lb - 1, size, p))
}

#' `E[Y | lb <= Y <= ub]` of a negative binomial count, the lower bound
#' inclusive: `y f(y; r, mu) = mu f(y - 1; r + 1, mu (r + 1) / r)`, the
#' same success probability with one more success.
#'
#' @noRd
nbinom_trunc_mean <- function(mu, size, lb, ub) {
  mu1 <- mu * (size + 1) / size
  num <- mu * (stats::pnbinom(ub - 1, size = size + 1, mu = mu1) -
                 stats::pnbinom(lb - 2, size = size + 1, mu = mu1))
  num / (stats::pnbinom(ub, size = size, mu = mu) -
           stats::pnbinom(lb - 1, size = size, mu = mu))
}

# --- Truncated draws by inversion ---------------------------------------

#' Both log tails of a family's CDF at `q`, numerically, for draws: from
#' `cdf_log` where the family declares it, else from `lcdf` and `lccdf`.
#' A tail that is `NaN` (a `q` outside a support the slot does not
#' guard) is read as the edge of the support on that side.
#'
#' @noRd
sim_cdf_pair <- function(fam, q, dp, av, extra) {
  n <- length(q)
  if (!is.null(fam[["cdf_log"]])) {
    lp <- fam_log_pair(fam, q, dp, av, extra, n)
    lF <- as.numeric(lp[["lF"]])
    lS <- as.numeric(lp[["lS"]])
  } else {
    Fp <- as.numeric(fam_lcdf(fam, q, dp, av, extra)) + numeric(n)
    lS <- if (has_lccdf(fam)) {
      as.numeric(fam_lccdf(fam, q, dp, av, extra)) + numeric(n)
    } else {
      log1p(-Fp)
    }
    lF <- ifelse(Fp > 0.5, log(-expm1(lS)), log(Fp))
  }
  list(lF = ifelse(is.nan(lF), -Inf, lF), lS = ifelse(is.nan(lS), 0, lS))
}

#' Draws of the rows `idx` from their family truncated to `[lb, ub]`, by
#' inverting the CDF: a uniform point of the window's mass, then
#' bisection on `q`. The point is placed in the tail the window's lower
#' edge leaves more of, on the log scale, as `log_window()` forms the
#' window, so a window deep in an upper tail is still drawn to its own
#' digits. A count is drawn as the smallest value whose CDF reaches the
#' point. `sim_response()` falls back to this for the rows rejection
#' sampling has not filled: a valid truncated fit may put a row's mean
#' far outside its window, and there rejection accepts a few draws in a
#' hundred.
#'
#' @noRd
sim_trunc_inverse <- function(fam, dpars, aterms, n, idx, lb, ub,
                              extra = NULL) {
  m <- length(idx)
  dp <- subset_obs(dpars, idx, n)
  av <- subset_obs(aterms, idx, n)
  disc <- identical(fam[["type"]], "discrete")
  lo <- if (disc) lb[idx] - 1 else lb[idx]
  if (disc) lo[lo < 0] <- -Inf
  hi <- ub[idx]
  pair <- function(q) sim_cdf_pair(fam, q, dp, av, extra)
  pl <- pair(ifelse(is.finite(lo), lo, 0))
  ph <- pair(ifelse(is.finite(hi), hi, 0))
  lFl <- ifelse(is.finite(lo), pl$lF, -Inf)
  lSl <- ifelse(is.finite(lo), pl$lS, 0)
  lFh <- ifelse(is.finite(hi), ph$lF, 0)
  lSh <- ifelse(is.finite(hi), ph$lS, -Inf)
  v <- stats::runif(m)
  useS <- lSl < log(0.5)
  # the point, as a log tail: S(lo) - v (S(lo) - S(hi)) or
  # F(lo) + v (F(hi) - F(lo))
  tS <- lSl + log1p(-v * (-expm1(lSh - lSl)))
  wF <- lFh + log(-expm1(pmin(lFl - lFh, 0)))
  tF <- ifelse(is.finite(lFl), pmax(lFl, log(v) + wF) +
                 log1p(exp(-abs(lFl - log(v) - wF))), log(v) + wF)
  # past the point: S below it, or F at or above it
  past <- function(q) {
    p <- pair(q)
    ifelse(useS, p$lS <= tS, p$lF >= tF)
  }
  # a bracket [a, b] with a short of the point and b past it, widened
  # outward from a finite edge where the window is open
  a <- ifelse(is.finite(lo), lo, ifelse(is.finite(hi), hi - 1, 0))
  b <- ifelse(is.finite(hi), hi, a + 1)
  step <- rep(1, m)
  for (k in 1:200) {
    out <- !is.finite(hi) & !past(b)
    if (!any(out)) break
    a[out] <- b[out]
    step[out] <- 2 * step[out]
    b[out] <- b[out] + step[out]
  }
  step <- rep(1, m)
  for (k in 1:200) {
    out <- !is.finite(lo) & past(a)
    if (!any(out)) break
    b[out] <- a[out]
    step[out] <- 2 * step[out]
    a[out] <- a[out] - step[out]
  }
  if (disc) {
    a <- floor(a)
    b <- ceiling(b)
    for (k in 1:200) {
      open <- b - a > 1
      if (!any(open)) break
      mid <- floor((a + b) / 2)
      p <- past(mid)
      b <- ifelse(open & p, mid, b)
      a <- ifelse(open & !p, mid, a)
    }
    return(b)
  }
  for (k in 1:200) {
    mid <- 0.5 * (a + b)
    p <- past(mid)
    b <- ifelse(p, mid, b)
    a <- ifelse(p, a, mid)
    if (all(b - a <= 1e-13 * pmax(1, abs(a)))) break
  }
  pmin(pmax(0.5 * (a + b), lb[idx]), ub[idx])
}
