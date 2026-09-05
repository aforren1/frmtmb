# The natural cubic spline basis of Royston and Parmar (2002), written
# so that it runs on an automatic-differentiation tape.
#
# The basis is a function of the RESPONSE, which is data, so nothing here
# has to differentiate in the ordinary case. It is written branch-free
# anyway, because one-step-ahead residuals promote the response to a
# parameter and a basis built with `pmax()` or `ifelse()` would stop the
# tape there: RTMB refuses a comparison on an AD type outright ("Comparison
# is generally unsafe for AD types"), and it has no CondExp.
#
# The truncated power is therefore spelled
#
#     (e)_+ = 0.5 * (e + abs(e))
#
# which `abs()` supports on an advector and which is exactly the positive
# part in double precision, with the derivative at e = 0 taken as zero.
# Divided differences of these truncated cubes agree with
# splines::splineDesign to 1e-12 with AD in the input, which is the
# measurement recorded in dev/spline-seam-proposal.md.

#' The positive part, branch free.
#'
#' @noRd
sp_pospart <- function(e) 0.5 * (e + abs(e))

#' A smooth positive floor, for a quantity that must stay positive but
#' that no link constrains.
#'
#' The derivative of the log cumulative hazard has to be positive for the
#' hazard to exist, and nothing in this parameterization holds it there
#' (flexsurv does not hold it either). Where it goes negative the log
#' density must still be a FINITE number: `-Inf` differentiates to `NaN`
#' and stops the fit, and inside a mixture one `NaN` takes every other
#' component with it. This is the drift-diffusion package's floor idiom,
#' spelled smoothly so that it also has a derivative.
#'
#' The obvious spelling, `0.5 * (u + sqrt(u^2 + eps2))`, is right on
#' paper and WRONG in double precision on the side it exists for: at
#' `u = -1e10` the square root rounds to `1e10` and the sum cancels to
#' exactly zero, so the floor returns 0 and `log()` returns `-Inf`
#' after all. Measured: 632 of 686 rows of the test data reached `-Inf`
#' that way while the optimizer was still finding the scale.
#'
#' The identity `(u + s)(s - u) = eps2` gives the cancellation-free form
#' on each side: `u + s` is `a + s` when `u` is positive and
#' `eps2 / (a + s)` when it is negative. `sign()` then SELECTS between
#' the two rather than combining them, because adding the tiny branch to
#' the large one loses it to rounding just as surely as the original sum
#' did. The two agree at `u = 0`, where both are `sqrt(eps2) / 2`, so the
#' result is continuous and `sign()`'s own zero derivative is the right
#' one.
#'
#' `eps2 = 1e-14` is the measured compromise. The floor changes
#' `log(u)` by 2.4e-15 at `u = 1` and by 1.1e-13 at `u = 0.1`, and it
#' returns 7e-17 rather than 0 at `u = -35.75`, which is the value that
#' produced 647 rows of `-Inf` before it was written this way.
#'
#' @noRd
sp_floor_pos <- function(u, eps2 = 1e-14) {
  sg <- sign(u)
  s <- sqrt(u * u + eps2)
  w <- abs(u) + s
  0.25 * ((1 + sg) * w + (1 - sg) * (eps2 / w))
}

#' A smooth ceiling, branch free.
#'
#' `min(u, cap)` written so that a tape can record it. Exact in double
#' precision wherever `|cap - u|` exceeds 6.7e-08, which is everywhere
#' except at the cap itself.
#'
#' @noRd
sp_cap <- function(u, cap, eps2 = 1e-30) {
  d <- cap - u
  cap - 0.5 * (d + sqrt(d * d + eps2))
}

#' `log(1 + exp(u))` without the overflow.
#'
#' Written as `max(u, 0) + log1p(exp(-|u|))`, so the exponential's
#' argument is never positive and the result is exact at both tails.
#' Spelled out rather than left as `log(1 + exp(u))` because that form
#' returns `Inf` past `u = 709` and an `Inf` in a log likelihood is a
#' `NaN` gradient one line later.
#'
#' @noRd
sp_log1pexp <- function(u) {
  a <- abs(u)
  0.5 * (u + a) + log1p(exp(-a))
}

#' Hold a probability strictly inside (0, 1), branch free.
#'
#' frmtmb forms the right-censored contribution as `log(Fub - F(y))`,
#' which is `log(1 - F)` without truncation, so a family that hands back
#' an `F` of exactly 1 hands back `-Inf`, and `-Inf` differentiates to
#' `NaN` and stops the fit. A Royston-Parmar spline reaches `F == 1` in
#' double precision as soon as the cumulative hazard passes about 37,
#' which the optimizer walks through while it is still finding the
#' scale: measured on flexsurv's `bc` data, a slope of 30 in log time
#' puts 572 of 686 rows there.
#'
#' The squeeze is the same smooth positive part used on the hazard, once
#' at each end. With `eps2 = 1e-30` it is EXACT in double precision
#' wherever the distance to the boundary exceeds 6.7e-08, which is a
#' survival probability far below anything a fit reports, and it puts
#' `1 - F` at 5e-16 rather than 0 where the double runs out.
#'
#' @noRd
sp_squeeze <- function(p, eps2 = 1e-30) {
  lo <- 0.5 * (p + sqrt(p * p + eps2))
  1 - 0.5 * ((1 - lo) + sqrt((1 - lo) * (1 - lo) + eps2))
}

#' The Royston and Parmar natural cubic spline basis.
#'
#' Column 1 is the constant, column 2 is `x` itself, and column `j + 2`
#' is the natural cubic term at interior knot `j`, constrained to be
#' linear beyond both boundary knots. This is
#' `flexsurv:::basis(knots, x)` term for term, so a `gamma` vector
#' estimated by one package means the same curve in the other.
#'
#' @param knots The full knot vector, boundary knots included, on the
#'   `x` scale (log time), sorted.
#' @param x Points to evaluate at.
#' @return An `length(x)` by `length(knots)` matrix.
#' @noRd
sp_rp_basis <- function(knots, x) {
  nk <- length(knots)
  out <- list(x - x + 1, x)
  if (nk > 2L) {
    kmin <- knots[1L]
    kmax <- knots[nk]
    for (j in seq_len(nk - 2L)) {
      kj <- knots[j + 1L]
      lam <- (kmax - kj) / (kmax - kmin)
      out[[j + 2L]] <- sp_pospart(x - kj)^3 -
        lam * sp_pospart(x - kmin)^3 -
        (1 - lam) * sp_pospart(x - kmax)^3
    }
  }
  out
}

#' The derivative of that basis with respect to `x`.
#'
#' `flexsurv:::dbasis(knots, x)`.
#'
#' @noRd
sp_rp_dbasis <- function(knots, x) {
  nk <- length(knots)
  out <- list(x - x, (x - x) + 1)
  if (nk > 2L) {
    kmin <- knots[1L]
    kmax <- knots[nk]
    for (j in seq_len(nk - 2L)) {
      kj <- knots[j + 1L]
      lam <- (kmax - kj) / (kmax - kmin)
      out[[j + 2L]] <- 3 * sp_pospart(x - kj)^2 -
        3 * lam * sp_pospart(x - kmin)^2 -
        3 * (1 - lam) * sp_pospart(x - kmax)^2
    }
  }
  out
}

#' `sum_j basis_j(x) gamma_j`, with `gamma` a list of per-row vectors.
#'
#' The coefficients are DISTRIBUTIONAL PARAMETERS, so each one is a
#' length-n vector rather than a scalar: that is what lets a covariate
#' reach `gamma1` and become a time-varying effect. The sum is written
#' as a loop over the basis rather than as a matrix product for exactly
#' that reason.
#'
#' @noRd
sp_rp_eta <- function(bas, gam) {
  eta <- bas[[1L]] * gam[[1L]]
  for (j in seq_along(bas)[-1L]) eta <- eta + bas[[j]] * gam[[j]]
  eta
}
