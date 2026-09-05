# Pieces more than one family in this package needs, kept in one place so
# that the second family to want them does not grow a second copy.

#' `exp(la) - exp(lb)` from the two LOGS, without overflow and without
#' cancellation.
#'
#' Both differences this package forms in a tail are differences of two
#' exponentials whose logs are known exactly: the normal density
#' difference in the racing diffusion model's density, and the
#' difference of the two boundary terms in its survival. Written
#' directly as `exp(la) - exp(lb)` either overflows, when one log is
#' large and positive, or underflows one term to zero while the other
#' stays finite, which turns a small difference into the larger term
#' alone.
#'
#' Anchoring at the larger log fixes both at once. Subtract the anchor
#' from each: one argument becomes exactly zero and the other is at most
#' zero, so neither `exp()` can overflow. Then use `expm1` on both rather
#' than `exp`, because `expm1` of the anchor's own argument is exactly
#' zero: the difference is the non-anchored term alone, computed to full
#' relative precision, rather than `1 - something-near-1`.
#'
#' `ddm_floor(la, lb)` is the anchor. It is written as a floor and used
#' here as a maximum, which is the same function: the helper is
#' `max(x, lo)` for any `lo`, and nothing in it requires `lo` to be a
#' constant.
#'
#' @noRd
ddm_expdiff <- function(la, lb) {
  m <- ddm_floor(la, lb)
  exp(m) * (expm1(la - m) - expm1(lb - m))
}

#' Log of the standard normal density, as a plain expression.
#'
#' `log(dnorm(x))` underflows to `-Inf` for `|x|` past about 38, which is
#' inside the range a race's loser reaches. The log is an exact
#' polynomial and never does.
#'
#' @noRd
ddm_lphi <- function(x) -0.5 * log(2 * pi) - 0.5 * x * x

#' Refuse a `dec()` term on a race family.
#'
#' `dec()` carries a two-level boundary indicator. A race of `n`
#' accumulators needs a winner in `1..n`, so the term cannot mean
#' anything to one, and frmtmb has no seam for "this family does not
#' accept that addition term": `required_aterms` is a conjunction of
#' what the density NEEDS, with no complementary allow-list, so an extra
#' term is never refused by declaration and each family that wants to
#' refuse one writes the check itself.
#'
#' Left unrefused the term is not merely ignored, it is invisible.
#' Measured on `lba(3)` over 300 rows, `rt | dec(two) + vint(choice)`
#' and `rt | vint(choice)` give fixed effects whose largest absolute
#' difference is EXACTLY ZERO and an identical log likelihood, with no
#' warning at any point, while the fitted object carries both `dec` and
#' `vint1` in its aterm values. The user's column travels all the way
#' into the fit and changes nothing.
#'
#' Shared by [lba()] and [rdm()] so that the two sibling race families
#' cannot drift apart on it, and so that this is ONE condition-message
#' template with the family name filled in at run time rather than two
#' near-identical ones.
#'
#' @noRd
ddm_refuse_dec <- function(what, n, aterms) {
  if (!is.null(aterms[["dec"]])) {
    stop(what, ": dec() carries a two-level boundary indicator, and ",
         "this family needs to know which of ", n, " accumulators ",
         "reached the threshold first. The term cannot mean anything ",
         "here, so it is refused rather than silently dropped: a model ",
         "ported over from wiener(), where dec() IS the spelling, ",
         "would otherwise fit while quietly ignoring it. The winner ",
         "travels through vint() as a whole number from 1 to ", n, ".",
         call. = FALSE)
  }
  invisible(NULL)
}

#' Fit a bounded non-decision-time link to the observed response.
#'
#' The density of every family here is zero at and below `ndt`, so the
#' likelihood has a hard edge at `ndt = min(rt)`. A log link would let
#' the optimizer walk over it; a logit scaled onto `(0, ub)` makes the
#' constraint structural, which is what [wiener()] and [lba()] both do.
#'
#' Shared rather than copied because two families in this package want
#' the same link and the same refusal, and the refusal carries the
#' family's name at run time so that the message a user sees still names
#' the family they wrote.
#'
#' @noRd
ddm_ndt_finalize <- function(fam, y, max_ndt, what) {
  ub <- if (is.null(max_ndt)) min(y) else max_ndt
  if (!is.null(max_ndt) && ub > min(y)) {
    stop(what, ": max_ndt = ", format(ub), " is above the fastest ",
         "response (", format(min(y)), "). Nothing can be observed ",
         "before the non-decision time, so a bound above the fastest ",
         "response admits values at which that trial has no ",
         "likelihood.", call. = FALSE)
  }
  fam[["links"]][["ndt"]] <- list(
    name = paste0("scaled_logit(0, ", signif(ub, 4), ")"),
    linkfun = function(mu) log(mu / (ub - mu)),
    linkinv = function(eta) ub / (1 + exp(-eta)),
    mu_eta = function(eta) {
      p <- 1 / (1 + exp(-eta))
      ub * p * (1 - p)
    })
  fam
}
