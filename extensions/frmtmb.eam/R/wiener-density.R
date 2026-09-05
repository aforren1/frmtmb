# The Wiener first-passage time density, Navarro and Fuss (2009).
#
# Two series converge to the same function. The small-time series is
# cheap when the normalized time u = t / a^2 is small and useless when
# it is large; the large-time series is the other way round. Navarro
# and Fuss pick one per evaluation with a criterion on u.
#
# An AD tape cannot make that choice. u depends on the boundary
# separation, which is a PARAMETER, so "if (u < kappa)" is a branch on
# a parameter: RTMB advectors carry no comparison operators, so it does
# not merely give a wrong gradient, it fails to tape at all. A branch on
# the response alone would be tape-safe but wrong, because the same t
# lands on either side of the criterion depending on a.
#
# So both series are evaluated at a fixed generous truncation and
# combined with a smooth weight in log(u). The weight is Navarro and
# Fuss's criterion with the step replaced by a logistic: it saturates to
# exactly 0 and exactly 1 in double precision outside a narrow band, and
# inside the band both series are accurate, so any convex combination of
# them is accurate too. Blending the LOG densities rather than the
# densities avoids underflow and makes the combination exact wherever
# the two logs agree.
#
# The cost is that every evaluation pays for both series. That is the
# price of a density that is correct over the whole range of normalized
# times rather than over the part of it a fixed truncation of one series
# happens to cover.

# ------------------------------------------------------- tuning constants
#
# ddm_ks: the small-time series runs k = -ddm_ks .. ddm_ks.
# ddm_kl: the large-time series runs k = 1 .. ddm_kl.
# ddm_u0, ddm_us: center and scale of the logistic weight in log(u).
#
# The values are measured, not guessed. See the accuracy tests in
# tests/testthat/test-density.R and the report in dev-findings.md: both
# series are accurate to better than 1e-12 relative on the whole band
# where the weight is not saturated.
ddm_ks <- 12L
ddm_kl <- 30L
ddm_u0 <- 0.35
ddm_us <- 0.12

#' A smooth positive floor.
#'
#' The series sums are positive wherever their own series is the
#' accurate one and may be anything at all, including negative, where it
#' is not. `log()` of a negative number is `NaN`, and `0 * NaN` is `NaN`,
#' so the saturated weight would not be enough on its own to keep the
#' garbage branch out of the answer. This maps the sum to its magnitude
#' without a comparison and without a kink at zero. In the regime where
#' a series is the accurate one its sum is of order one, so the floor is
#' inactive there.
#'
#' @noRd
ddm_pos <- function(x) sqrt(x * x + 1e-40)

#' Log of the normalized small-time density, `2 * ddm_ks + 1` terms.
#'
#' The series is written with the `k = 0` exponent factored out. The raw
#' spelling multiplies `(2 * pi * u^3)^(-1/2)`, which overflows for
#' small u, by a sum that underflows there; this form has the overflow
#' and the underflow cancel algebraically before either is evaluated, so
#' the sum tends to `w` rather than to `0 / Inf` as u goes to zero.
#'
#' Plain arithmetic only: no `c()`, no `[<-`, so it needs no
#' `frmtmb_ad_overload()` bindings of its own and it runs on advectors
#' and on doubles alike.
#'
#' @noRd
ddm_log_gs <- function(u, w, K = ddm_ks) {
  s <- w                            # the k = 0 term, exactly
  for (k in seq_len(K)) {
    s <- s + (w + 2 * k) * exp(-2 * k * (w + k) / u) +
      (w - 2 * k) * exp(2 * k * (w - k) / u)
  }
  -0.5 * log(2 * pi) - 1.5 * log(u) - w * w / (2 * u) + log(ddm_pos(s))
}

#' Log of the normalized large-time density, `ddm_kl` terms.
#'
#' The `k = 1` exponent is factored out for the same reason: the whole
#' sum decays like `exp(-pi^2 u / 2)` and would underflow for large u,
#' while the factored sum tends to `sin(pi * w)`, which is bounded away
#' from zero for a bias strictly inside the boundaries.
#'
#' @noRd
ddm_log_gl <- function(u, w, K = ddm_kl) {
  h <- pi * pi * u / 2
  s <- sin(pi * w)                  # the k = 1 term, exactly
  for (k in 2:K) {
    s <- s + k * exp(-(k * k - 1) * h) * sin(k * pi * w)
  }
  log(pi) - h + log(ddm_pos(s))
}

#' Hold a quantity above a floor, without a comparison.
#'
#' `max(x, lo)` written as `lo + ((x - lo) + |x - lo|) / 2`, which keeps
#' `lo` OUTSIDE the cancelling sum. The other maximum spelling,
#' `(x + lo + |x - lo|) / 2`, is wrong at the magnitudes the floors here
#' use: when `lo` is far below the rounding of `x`, the `+ lo` and the
#' `- lo` inside the absolute value both vanish, the expression collapses
#' to `(x + |x|) / 2`, and a negative `x` returns exactly zero, the value
#' the floor exists to prevent. Measured with `lo = 1e-300`, that form
#' returns 0 at `x = -1e-17`; this one returns `lo`.
#'
#' This one helper serves every family in the package. Where `lo` is a
#' denormal or 1e-300, it is bit-for-bit inert on every value in range,
#' so a density is unchanged wherever it was already right. Where `lo`
#' is 1e-10 or 1e-12, as in the accumulator race's normalizers, it is an
#' exact maximum: `x` passes through untouched above the floor, where a
#' positive part plus the floor would have added `lo` to every value.
#'
#' `abs()` rather than a comparison because the argument descends from
#' a parameter, and an RTMB tape records no comparisons. The derivative
#' is 0 below the floor and 1 above it, which is what makes a row the
#' model cannot reach flat rather than undefined: its log density is a
#' large finite negative number, its contribution to a mixture's
#' log-sum-exp exponentiates to exactly 0, and its gradient is exactly
#' 0. `-Inf` would not do, because `-Inf` differentiates to `NaN`.
#'
#' @noRd
ddm_floor <- function(x, lo) lo + 0.5 * ((x - lo) + abs(x - lo))

# The floor the normalized time is held at. It is the smallest positive
# double there is, which is not a fudge factor but the whole mechanism:
# at u = ddm_u_floor the term -w^2 / (2u) overflows to -Inf while
# -1.5 log(u) stays finite near +1100 and every series term underflows
# to zero, so the log density is exactly -Inf rather than NaN.
ddm_u_floor <- 1e-320

# The floor the surviving share of the non-decision-time range is held
# at. Not the same magnitude as ddm_u_floor and not for the same reason:
# this one is passed to log(), and log() of a denormal has a derivative
# of 1e320, which overflows to Inf and then meets a zero from the
# floor's own flat branch. Inf times zero is NaN, which is the failure
# the floor exists to prevent. At 1e-300 the derivative is finite, and
# the value it produces is still far below anything a fit can reach.
ddm_share_floor <- 1e-300

#' Log first-passage density at the LOWER boundary.
#'
#' `t` is decision time (response time less the non-decision time), `v`
#' the drift rate, `a` the boundary separation and `w` the relative
#' start point in (0, 1). The upper boundary is not a second series: it
#' is this one at `v -> -v`, `w -> 1 - w`.
#'
#' `t <= 0` gives `-Inf`. The density there is zero, so a log density of
#' `-Inf` is the answer, and it is the one a mixture's log-sum-exp can
#' use: `NaN`, which is what the unfloored arithmetic produces, is not
#' an answer at all and propagates through every component. The family's
#' bounded non-decision-time link is what keeps a bare fit away from
#' there; see `wiener()`.
#'
#' @noRd
ddm_lpdf_lower <- function(t, v, a, w) {
  ur <- t / (a * a)
  u <- ddm_floor(ur, ddm_u_floor)
  # Navarro and Fuss's series criterion, with the step replaced by a
  # logistic so that the choice differentiates. tanh rather than the
  # exp spelling of the same curve: tanh saturates to exactly +/-1 in
  # double precision without ever forming Inf, and an Inf on the tape
  # poisons the reverse pass even where its weight is zero.
  lam <- 0.5 * (1 + tanh((log(u) - log(ddm_u0)) / ddm_us))
  lg <- (1 - lam) * ddm_log_gs(u, w) + lam * ddm_log_gl(u, w)
  -v * a * w - v * v * t / 2 - 2 * log(a) + lg
}

#' The same density with the boundary supplied as data.
#'
#' `up` is 1 for a response at the upper boundary and 0 for the lower
#' one. It reaches the density as an addition term, so it is data, and
#' the reflection is arithmetic on it rather than a branch.
#'
#' @noRd
ddm_lpdf_both <- function(t, v, a, w, up) {
  ddm_lpdf_lower(t, v * (1 - 2 * up), a, w + up * (1 - 2 * w))
}

# --------------------------------------------- across-trial variability
#
# Ratcliff's full diffusion model draws the drift rate, the start point
# and the non-decision time afresh on every trial, so the likelihood of
# one response time is the Wiener density averaged over those three
# distributions. The average is not a random effect: nothing is shared
# between trials, there is no level to estimate, and the model has one
# of these integrals per ROW whether or not it has any grouping factor
# at all. It therefore belongs inside the density, which is where it is.
#
# Three integrals, three different answers:
#
#   drift        Gaussian, and the drift enters the log density
#                LINEARLY in one term and QUADRATICALLY in one other.
#                That is an exponential-quadratic against a normal,
#                which integrates in closed form. Exact, no nodes.
#   start point  uniform over an interval strictly inside (0, 1), with
#                an analytic integrand on it: Gauss-Legendre.
#   non-decision uniform, and the only one of the three whose range is
#                cut by the data, because the density is zero for a
#                non-decision time at or past the response time.

#' A minimum written without a comparison.
#'
#' `min(x, cap)` where `cap` is data and `x` is a parameter. `abs()` is
#' the only order-sensitive operation an RTMB tape records, so it is the
#' only way to write one. The kink at `x = cap` costs nothing here: it
#' is always the upper end of a range over which the Wiener density
#' vanishes to every order, so the quantity it feeds is flat there.
#'
#' @noRd
ddm_smin <- function(x, cap) 0.5 * (x + cap - abs(cap - x))

#' Log first-passage density at the LOWER boundary with the drift rate
#' integrated over `N(v, sv^2)`.
#'
#' The drift appears in the log density only as `-v a w - v^2 t / 2`:
#' the series that carries the shape of the distribution is the
#' DRIFTLESS one, and the drift multiplies it by an
#' exponential-quadratic. So the trial-to-trial average over the drift
#' is a Gaussian integral of `exp(-A v - B v^2)`, which completes the
#' square:
#'
#'   mean of exp(-v a w - v^2 t / 2) =
#'     (1 + t sv^2)^(-1/2) exp( ((a w)^2 sv^2 - 2 v a w - v^2 t)
#'                              / (2 (1 + t sv^2)) )
#'
#' At `sv = 0` this is `exp(-v a w - v^2 t / 2)` exactly, term for term,
#' so the variability family reduces to the plain one with no tolerance
#' to argue about. The result is EXACT, not a quadrature: there is no
#' node count to choose for the drift and no accuracy to trade.
#'
#' `t` is floored the same way the plain density floors it, and the
#' floored value is what the drift factor uses as well, so that a
#' decision time at or below zero gives `-Inf` rather than the `NaN`
#' that `log(1 + t sv^2)` would produce at a sufficiently negative `t`.
#'
#' @noRd
ddm_lpdf_lower_sv <- function(t, v, a, w, sv) {
  aa <- a * a
  ur <- t / aa
  u <- ddm_floor(ur, ddm_u_floor)
  te <- u * aa
  lam <- 0.5 * (1 + tanh((log(u) - log(ddm_u0)) / ddm_us))
  lg <- (1 - lam) * ddm_log_gs(u, w) + lam * ddm_log_gl(u, w)
  s2 <- sv * sv
  d <- 1 + te * s2
  z <- a * w
  (z * z * s2 - 2 * v * z - v * v * te) / (2 * d) - 0.5 * log(d) -
    2 * log(a) + lg
}

#' The full-DDM log density for one row, boundary supplied as data.
#'
#' `y` is the response time, not the decision time: the non-decision
#' time is integrated over here, so it cannot be subtracted before the
#' integral the way the plain density subtracts it.
#'
#' The reflection to the upper boundary is applied once, at the top, and
#' commutes with all three variability distributions: a normal drift
#' reflects to a normal drift with the opposite mean, and a uniform
#' start point at `w0` reflects to a uniform of the same width at
#' `1 - w0`.
#'
#' `nd` carries the node sets, `st_on` says whether the non-decision
#' time is integrated at all, and `delta` is the margin by which the
#' non-decision time range is held below each row's own response time.
#' All three are decided when the family is built.
#'
#' @noRd
ddm_lpdf_var <- function(y, v, a, w, t0, sv, sz, st, up, nd, st_on,
                         delta) {
  vv <- v * (1 - 2 * up)
  w0 <- w + up * (1 - 2 * w)

  # The live part of the non-decision time's range. Past `y` the
  # decision time would be negative and the density is identically
  # zero, so the integral's upper limit is capped rather than its
  # integrand tested. The cap sits `delta` below the response time so
  # that every node is strictly inside the support: an unreachable row
  # then gets a log density of about -delta^-1, which exponentiates to
  # exactly zero in a mixture and differentiates to exactly zero, where
  # a true -Inf would differentiate to NaN.
  cap <- y - delta
  half <- 0.5 * st
  lo <- ddm_smin(t0 - half, cap)
  hi <- ddm_smin(t0 + half, cap)
  span <- hi - lo
  # The share of the uniform's own range that survives the cap. The
  # quadratures below average rather than integrate, so this is the one
  # place the width of the non-decision time distribution appears.
  #
  # A row the whole range has moved past keeps none of it, and log(0) is
  # the -Inf whose gradient is NaN. The share is therefore held at
  # `ddm_share_floor`, which turns that row's log density into a finite
  # number - zero when exponentiated, and flat, so its gradient is zero
  # rather than undefined. Held on the SHARE and not on the width,
  # because the share of an uncut range is exactly 1 whatever the width
  # is, so the floor cancels on every row that keeps its range rather
  # than perturbing a narrow one.
  #
  # `span` is a difference of two independently rounded evaluations of
  # the same cap, so on a row whose whole range has moved past it the
  # share is not reliably zero: it is zero or a tiny NEGATIVE, and a
  # negative is what breaks a smooth maximum. `ddm_floor()` is the form
  # that survives it.
  lfrac <- if (st_on) {
    log(ddm_floor(span / st, ddm_share_floor))
  } else 0

  # A log-sum-exp over the node grid, with the reference taken as the
  # MAXIMUM over the nodes rather than at any chosen one of them.
  #
  # The reference cancels algebraically, so it changes nothing but the
  # conditioning - and the conditioning is the whole difficulty. Both
  # failure modes are one-sided, and picking a corner only ever fixes
  # one of them. A reference below the largest node overflows: the log
  # density falls by v^2 span / 2 across the range, which passes 709 at
  # a drift near 120 for a range of a tenth of a second, and exp() of
  # the difference is then +Inf, so the log likelihood comes out +Inf
  # and an optimizer chases it. A reference above every node underflows
  # the lot to zero, and log(0) is -Inf with a NaN gradient.
  #
  # The maximum has neither problem by construction: every exponent is
  # at most zero, and one of them is exactly zero, so the sum lies in
  # [w_max, 1] whatever the parameters do. It is computed by pairwise
  # smooth maxima, which is exact and needs no comparison, and the node
  # values are kept rather than recomputed.
  lps <- vector("list", length(nd$st$x) * length(nd$sz$x))
  wts <- numeric(length(lps))
  m <- 0L
  for (j in seq_along(nd$st$x)) {
    tau <- lo + span * nd$st$x[[j]]
    for (i in seq_along(nd$sz$x)) {
      om <- w0 + sz * (nd$sz$x[[i]] - 0.5)
      m <- m + 1L
      lps[[m]] <- ddm_lpdf_lower_sv(y - tau, vv, a, om, sv)
      wts[m] <- nd$st$w[[j]] * nd$sz$w[[i]]
    }
  }
  lref <- lps[[1L]]
  for (m in seq_along(lps)[-1L]) {
    d <- lref - lps[[m]]
    lref <- 0.5 * (lref + lps[[m]] + abs(d))
  }
  acc <- 0
  for (m in seq_along(lps)) acc <- acc + wts[m] * exp(lps[[m]] - lref)
  lref + log(acc) + lfrac
}
