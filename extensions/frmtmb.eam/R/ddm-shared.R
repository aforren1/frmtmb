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

#' The fastest response, in seconds, above which the data is read as a
#' unit mistake rather than as a slow task.
#'
#' Not a plausible single response time. It is a ceiling on the fastest
#' response IN THE WHOLE DATA SET, which is a much rarer thing: one slow
#' trial does not reach it and cannot.
#'
#' Where 20 comes from, measured on simulated wiener data
#' (`dev/smallitems-findings.md` carries the tables). Over 192 cells at
#' the parameters a two-choice task actually produces, boundary 2 to 5,
#' non-decision time 0.5 to 5 seconds, drift 0.3 to 2, at 12 to 80
#' trials, the fastest response ran from 0.599 to 6.605 seconds and
#' NONE false alarmed; and over 81 cells at the standard published
#' range, all 81 are caught when the same data is read as milliseconds.
#' In the other direction the miss rate is zero at every threshold up to
#' 100 and rises at 200, so 20 sits a factor of five inside the band
#' where nothing is missed.
#'
#' It CAN fire on a correct model, and the rate is worth knowing rather
#' than assuming. The exposed design is a slow one with a short
#' session. At 20 trials the rate is 0.045 where the fastest response
#' averages 13.4 seconds, 0.155 at 15.7, 0.690 at 22.7 and 1.000 at
#' 41.7; sixty trials of the same tasks give 0.000, 0.000, 0.285 and
#' 0.995, because the fastest of more draws settles toward the floor.
#'
#' Quoted against the FASTEST response and not the median, because the
#' median does not determine the rate: at a median near 50 seconds the
#' rate runs from 0.185 to 0.935 depending on whether the task is slow
#' from a far boundary or from weak evidence, and the fastest response
#' tracks it in every one of those cells. The cost falls on
#' deliberation, insight and matrix-reasoning designs, whose
#' non-decision time is 1 to 2 seconds rather than 20.
#'
#' @noRd
ddm_seconds_ceiling <- 20

#' Warn when the response is not on the scale the defaults assume.
#'
#' Every default in this package reads the response as SECONDS. The
#' starting values, `gddm_control(dt = 0.01)` and the window `t_max`
#' takes from the data are all absolute times, and nothing in any of the
#' five likelihoods refuses milliseconds: the fit converges and reports
#' a boundary separation three orders of magnitude out. That is the
#' silent wrong answer this package ranks first, so it is worth a
#' warning even at the cost of an occasional false alarm.
#'
#' A warning and not a refusal, because a design with no trial under 20
#' seconds is possible, if rare. It carries a class so that such a
#' design can silence this one condition without also hiding the
#' convergence warnings beside it.
#'
#' Called from each family's `valid_y` rather than once from here, so
#' that the sentence a user sees names the family they wrote.
#'
#' @noRd
ddm_check_units <- function(y, what) {
  if (!length(y)) return(invisible(NULL))
  lo <- min(y)
  if (!is.finite(lo) || lo <= ddm_seconds_ceiling) return(invisible(NULL))
  warning(warningCondition(paste0(
    what, ": the fastest response in these data is ",
    format(lo, digits = 4), ", and every default in this package reads ",
    "the response as a time in SECONDS. Milliseconds is the usual ",
    "cause: a millisecond clock puts a typical response near 500, and a ",
    "task in which no trial at all finishes within ",
    ddm_seconds_ceiling, " seconds is rare. Divide the response by 1000 ",
    "if that is what happened, which puts these times between ",
    format(lo / 1000, digits = 4), " and ",
    format(max(y) / 1000, digits = 4), " seconds. If the task really is ",
    "this slow then the fit is correct and this warning is its only ",
    "cost; it carries the class frmtmb_eam_units_warning so that it can ",
    "be silenced on its own"),
    class = "frmtmb_eam_units_warning"))
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
