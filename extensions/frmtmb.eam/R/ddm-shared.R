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

#' What each family in this package reads, in one place.
#'
#' The constructor passes its entry to `accepts_aterms` and the
#' compatibility rows derive the refusals from the same entry, so the
#' table cannot promise a term frame assembly refuses, nor refuse one it
#' takes. Written here rather than read back off the family objects
#' because the rules builder runs inside `.onLoad()`, where building
#' five families to ask them one question each is work nobody needs on a
#' package load.
#'
#' Spelled in TERMS, without parentheses, which is the vocabulary
#' `frmtmb_family(accepts_aterms =)` uses.
#'
#' @noRd
ddm_accepts <- list(
  wiener     = c("dec", "vint", "weights"),
  gddm       = c("dec", "vint", "vreal", "weights"),
  lba        = c("vint", "weights"),
  rdm        = c("vint", "weights", "cens", "trunc"),
  wiener_gng = c("dec", "vreal", "weights", "cens"))
