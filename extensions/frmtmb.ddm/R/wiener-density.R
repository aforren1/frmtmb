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

#' Log first-passage density at the LOWER boundary.
#'
#' `t` is decision time (response time less the non-decision time), `v`
#' the drift rate, `a` the boundary separation and `w` the relative
#' start point in (0, 1). The upper boundary is not a second series: it
#' is this one at `v -> -v`, `w -> 1 - w`.
#'
#' `t <= 0` gives `NaN`, which the optimizer reads as a rejected step.
#' The family's bounded non-decision-time link is what keeps the fit
#' away from there; see `wiener()`.
#'
#' @noRd
ddm_lpdf_lower <- function(t, v, a, w) {
  u <- t / (a * a)
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
