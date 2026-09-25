# The uniform contaminant of wiener(contaminant = TRUE), and the refusal
# for a censored model the distribution function does not cover.
#
# WHY INSIDE THE FAMILY. A contaminant is a second process that
# produced some of the trials: a guess, a lapse, a response to
# something other than the stimulus. Ratcliff and Tuerlinckx (2002)
# model it as a response time uniform over the observed range and a
# choice that is a coin flip. The uniform has NO free parameter, so it
# cannot be a brms family of its own (every brms family must estimate
# something), and in frmtmb it would have to be a mixture component
# with a response range fixed from data that mixture() never
# finalizes. Here it is one more term in the density and one more
# distributional parameter, `lambda`, the mixing proportion.
#
# THE JOINT DENSITY. The response is the pair (time, boundary), so the
# contaminant's density over the pair is
#
#   g(t, b) = 1 / (2 (hi - lo))   for lo <= t <= hi, either boundary
#
# and the row's density is (1 - lambda) f(t, b) + lambda g(t, b). The
# factor of one half is the coin flip. Without it the contaminant would
# put mass one on EACH boundary and the mixture would integrate to
# 1 + lambda.
#
# WHY THE RANGE IS THE OBSERVED ONE. It is the range the data could
# have come from, and it is fixed at frame assembly and kept, like the
# non-decision-time bound, so a refit on fewer rows scores against the
# range the fit used.

# The log density a row outside the contaminant's range is given, in
# place of -Inf. A difference against -Inf is NaN inside the smooth
# maximum below; against this it is finite, the exponential it feeds
# is exactly zero, and so is its gradient. A fitted row reaches it only
# under a `contaminant_range` narrower than the data; newdata can too.
ddm_cont_lout <- -1e300

#' Log density of the contaminant at each row.
#'
#' `y` is data here, so the range test is a comparison on data and
#' leaves no branch on the tape.
#'
#' @noRd
ddm_cont_lg <- function(y, crange) {
  if (anyNA(crange)) ddm_stop_cont_range()
  yy <- as.numeric(y)
  inside <- yy >= crange[1L] & yy <= crange[2L]
  ifelse(inside, -log(2 * (crange[2L] - crange[1L])), ddm_cont_lout)
}

#' Wrap a Wiener log density with the contaminant.
#'
#' `log((1 - lambda) f + lambda g)`, written as
#'
#'   lf + M + log((1 - lambda) exp(-M) + lambda exp(d - M))
#'
#' with `d = lg - lf` and `M = max(d, 0)`. Both exponents are at most
#' zero, so neither overflows, and at least one of the two terms is of
#' order one, so the log is finite for every row with `lambda` inside
#' its link. `M` is `(d + |d|) / 2`, which is EXACT in floating point:
#' it is `d + d` halved or `d - d`.
#'
#' The consequence that matters: at `lambda = 0` a row whose Wiener
#' density is at least the contaminant's (`d <= 0`) comes back as `lf`
#' BIT FOR BIT, because `M` is exactly zero and the log is of exactly
#' one. A row with `d > 0` comes back as `lf + M + log(exp(-M))`, which
#' is `lf` to within the rounding of `M`. That is the identity item 3.5
#' asks for, stated at the precision it holds.
#'
#' `(d - m)` is parenthesized on purpose. For a row below the
#' non-decision time under a stated `max_ndt`, `lf` is about -1.8e9, so
#' `d` is about 1.8e9 and `l + d` rounds `log(lambda)` to the ulp of
#' 1.8e9, 2.4e-7. The objective was then a staircase in `lambda` while
#' the tape gradient stayed smooth, and nlminb reported false
#' convergence on most recovery fits. `d - m` is exactly 0 there.
#'
#' The log of the sum is floored at 1e-300, which is inert on every
#' value a double can hold above that, and turns the one case that
#' could underflow (a row the Wiener density puts 745 log units below
#' the contaminant, at a `lambda` below exp(-745)) into a finite
#' barrier instead of a NaN gradient.
#'
#' @noRd
ddm_cont_lpdf <- function(lpdf, crange) {
  force(lpdf)
  force(crange)
  function(y, dpars, aterms) {
    lf <- lpdf(y, dpars, aterms)
    lc <- frmtmb::dpar_log_complement(dpars, "lambda", "logit")
    d <- ddm_cont_lg(y, crange) - lf
    m <- 0.5 * (d + abs(d))
    lf + m + log(ddm_floor(exp(lc[["l1m"]] - m) + exp(lc[["l"]] + (d - m)),
                           1e-300))
  }
}

#' The contaminant's distribution function at `q`, a probability.
#'
#' @noRd
ddm_cont_G <- function(q, crange) {
  if (anyNA(crange)) ddm_stop_cont_range()
  qq <- as.numeric(q)
  pmin(pmax((qq - crange[1L]) / (crange[2L] - crange[1L]), 0), 1)
}

#' The mixture's distribution function: `(1 - lambda) F + lambda G`.
#'
#' Formed on the probability scale because that is what frmtmb's `lcdf`
#' slot returns. The contaminant marginalizes over its coin flip, so
#' its distribution function over the response time alone is `G`, with
#' no factor of one half. A left- or interval-censored row is scored at
#' its known boundary, and there the contaminant's share is `G / 2`,
#' the coin flip landing on that boundary.
#'
#' @noRd
ddm_cont_lcdf <- function(lcdf, crange) {
  force(lcdf)
  force(crange)
  function(q, dpars, aterms) {
    lc <- frmtmb::dpar_log_complement(dpars, "lambda", "logit")
    cen <- aterms[["cens"]]
    gs <- if (is.null(cen)) 1 else 1 - 0.5 * as.numeric(cen == -1 | cen == 2)
    exp(lc[["l1m"]]) * lcdf(q, dpars, aterms) +
      exp(lc[["l"]]) * ddm_cont_G(q, crange) * gs
  }
}

#' The mixture's log survival: `log((1 - lambda) S + lambda (1 - G))`.
#'
#' On the log scale, anchored the way `ddm_cont_lpdf()` is and for the
#' same reason, so that a Wiener survival of 1e-30 at a deadline is
#' kept rather than lost against the contaminant's. `1 - G` is zero at
#' and past the top of the range, which is where a deadline usually
#' sits, so it is floored at 1e-300 before its log is taken: a
#' contaminant cannot outlast the range it is uniform over.
#'
#' @noRd
ddm_cont_lccdf <- function(lccdf, crange) {
  force(lccdf)
  force(crange)
  function(q, dpars, aterms) {
    ls <- lccdf(q, dpars, aterms)
    lc <- frmtmb::dpar_log_complement(dpars, "lambda", "logit")
    lsg <- log(pmax(1 - ddm_cont_G(q, crange), 1e-300))
    d <- lsg - ls
    m <- 0.5 * (d + abs(d))
    ls + m + log(ddm_floor(exp(lc[["l1m"]] - m) + exp(lc[["l"]] + (d - m)),
                           1e-300))
  }
}

#' The probability of the row's own boundary, off the tape.
#'
#' The contaminant's share of a row depends on it: given that a trial
#' ended at a boundary the Wiener process rarely reaches, it is more
#' likely to have been a guess. Plain model in closed form; under
#' `variability` the same average over the drift and the start point
#' that `ddm_mean_rt_var()` takes as its denominator.
#'
#' @noRd
ddm_boundary_prob_fn <- function(vv, nd) {
  if (!length(vv)) {
    return(function(dpars, aterms, n) {
      up <- rep_len(ddm_indicator(aterms), n)
      pu <- ddm_p_upper(rep_len(dpars[["mu"]], n), rep_len(dpars[["bs"]], n),
                        rep_len(dpars[["bias"]], n))
      ifelse(up == 1, pu, 1 - pu)
    })
  }
  gh <- ddm_gauss_hermite(21L)
  function(dpars, aterms, n) {
    g <- function(nm, default) {
      if (is.null(dpars[[nm]])) rep_len(default, n) else rep_len(dpars[[nm]], n)
    }
    v <- g("mu", 0); a <- g("bs", 1); w <- g("bias", 0.5)
    sv <- g("sv", 0); sz <- g("sz", 0)
    up <- rep_len(ddm_indicator(aterms), n)
    num <- numeric(n)
    den <- 0
    for (k in seq_along(gh[["x"]])) {
      nu <- v + sv * gh[["x"]][[k]]
      for (i in seq_along(nd[["sz"]][["x"]])) {
        om <- w + sz * (nd[["sz"]][["x"]][[i]] - 0.5)
        pu <- ddm_p_upper(nu, a, om)
        wt <- gh[["w"]][[k]] * nd[["sz"]][["w"]][[i]]
        num <- num + wt * ifelse(up == 1, pu, 1 - pu)
        den <- den + wt
      }
    }
    num / den
  }
}

#' The contaminant's posterior share of each row, given its boundary.
#'
#' `lambda / 2` of the mass is contaminant at each boundary and
#' `(1 - lambda) P(b)` is the diffusion's, so Bayes' rule is a ratio of
#' the two.
#'
#' @noRd
ddm_cont_share <- function(dpars, aterms, n, pb) {
  lam <- rep_len(dpars[["lambda"]], n)
  pw <- pb(dpars, aterms, n)
  0.5 * lam / (0.5 * lam + (1 - lam) * pw)
}

#' The conditional mean response time under the contaminant.
#'
#' @noRd
ddm_cont_mean <- function(mean_fn, pb, crange) {
  force(mean_fn)
  force(pb)
  force(crange)
  function(dpars, aterms) {
    base <- mean_fn(dpars, aterms)
    if (anyNA(crange)) ddm_stop_cont_range()
    pc <- ddm_cont_share(dpars, aterms, length(base), pb)
    pc * mean(crange) + (1 - pc) * base
  }
}

#' A draw under the contaminant, conditional on the row's boundary.
#'
#' Each row is a contaminant with its posterior share given the
#' boundary, and then draws from the uniform, or else from the
#' diffusion exactly as the plain family does.
#'
#' @noRd
ddm_cont_sim <- function(sim, pb, crange) {
  force(sim)
  force(pb)
  force(crange)
  function(dpars, aterms, n) {
    if (anyNA(crange)) ddm_stop_cont_range()
    out <- sim(dpars, aterms, n)
    pc <- ddm_cont_share(dpars, aterms, n, pb)
    hit <- stats::runif(n) < pc
    if (any(hit)) out[hit] <- stats::runif(sum(hit), crange[1L], crange[2L])
    out
  }
}

#' @noRd
ddm_stop_cont_range <- function() {
  frm_stop("wiener(contaminant = TRUE): the contaminant's range is not ",
           "set yet. It is the observed range of the response, which frm() ",
           "reads when it assembles the model frame; a family object used ",
           "outside frm() has not seen any data.", call. = FALSE)
}

#' @noRd
ddm_stop_cens_var <- function() {
  frm_stop("wiener: cens() and trunc() need the distribution function of ",
           "the response time, and this package has it for the plain model ",
           "only. Under `variability` it is the same series averaged over ",
           "the per-trial drift, start point and non-decision time, and the ",
           "drift average has no closed form there. Drop `variability`, or ",
           "the censoring.", call. = FALSE)
}

#' @noRd
ddm_stop_cens_trunc <- function() {
  frm_stop("wiener: left and interval censoring cannot be combined with ",
           "trunc(). A left- or interval-censored trial reached a known ",
           "boundary, so its probability is that boundary's defective ",
           "distribution function, while the truncation window is ",
           "normalized over both boundaries; frmtmb forms both from the ",
           "one distribution-function slot a family has, so it cannot ",
           "tell them apart. Right censoring with trunc() works.",
           call. = FALSE)
}

#' The contaminant's window when the user gives none.
#'
#' Only a model with a response DEADLINE has one, declared through
#' trunc(ub =): the deadline bounds the diffusion and the contaminant
#' alike, so the observed range is a window both processes filled.
#' Without a deadline the slowest diffusion trial sets the top, the
#' window grows with the sample, and lambda is estimated too low: on
#' one subject of 4000 trials with a true share of 0.049, the observed
#' range gave 0.0050 and covered on 5 of 25 seeds (review of
#' 2026-09-24). That is refused rather than fitted.
#'
#' Under a deadline the window runs from the fastest response to the
#' deadline itself, not to the slowest response. Both were measured on
#' one subject of 4000 trials, 25 seeds, contaminants uniform on `[0, 5]`
#' and a 5 s deadline (dev/phase3b-cont-window.R): lambda 0.0508 with
#' the deadline as the top against 0.0518 with the slowest response,
#' for a true recorded share of 0.0496, both covering on 25 of 25. The
#' deadline is also the one of the two that does not move with the
#' sample. A deadline that differs between rows gives no one window,
#' and is refused.
#'
#' @noRd
ddm_cont_default_range <- function(y, aterms) {
  ub <- aterms[["trunc_ub"]]
  if (!is.null(ub) && length(unique(as.numeric(ub))) > 1L) {
    frm_stop("wiener(contaminant = TRUE): trunc(ub = ) differs between ",
             "rows, so there is no one deadline to end the contaminant's ",
             "window at. Give `contaminant_range =`.", call. = FALSE,
             class = "frmtmb_eam_contaminant_range_error")
  }
  if (is.null(ub)) {
    frm_stop("wiener(contaminant = TRUE) needs `contaminant_range =`, the ",
             "response times the contaminant is uniform over, because this ",
             "model declares no response deadline. Without a deadline ",
             "bounding both processes, the slowest diffusion trial would set ",
             "the window, the window would grow with the sample, and lambda ",
             "would collapse toward zero. Give the task's response window, ",
             "for example contaminant_range = c(0, 5), or declare the ",
             "deadline with trunc(ub = ).", call. = FALSE,
             class = "frmtmb_eam_contaminant_range_error")
  }
  c(min(y), as.numeric(ub)[1L])
}

#' Refuse a response outside the contaminant's window, by name.
#'
#' The review of punch round 1 found a response past the window's end
#' fitted without a word: the contaminant gave that row density zero and
#' only the Wiener part was left to explain it. The interval-censored
#' row's upper edge is a response time too and is checked with it.
#'
#' @noRd
ddm_cont_check_range <- function(y, aterms, crange) {
  yy <- c(as.numeric(y), as.numeric(aterms[["cens_y2"]]))
  yy <- yy[is.finite(yy)]
  out <- yy[yy < crange[1L] | yy > crange[2L]]
  if (length(out)) {
    frm_stop("wiener(contaminant = TRUE): ", length(out), " response ",
             "time", if (length(out) > 1L) "s lie" else " lies",
             " outside contaminant_range = c(", format(crange[1L]), ", ",
             format(crange[2L]), "), for example ", format(out[1L]),
             ". The contaminant has density zero there. Widen ",
             "`contaminant_range =` to the task's response window.",
             call. = FALSE, class = "frmtmb_eam_contaminant_range_error")
  }
  invisible(NULL)
}

#' The refusals for what item 3.4 and 3.5 built for wiener() only.
#'
#' Without these, lba(contaminant = TRUE) failed with R's "unused
#' argument", and lba() or gddm() with cens() or trunc() failed with
#' frmtmb's generic message, which names neither the family nor the
#' reason. The slots below are filled with the refusal for the reason
#' ddm_family() gives: frmtmb checks that a censored model's family has
#' a distribution function before it reads the family's validator.
#'
#' @noRd
ddm_stop_not_built <- function(what, feature) {
  frm_stop(what, "(): ", feature, " is not built for this family. It is ",
           "built for wiener() only; see NEWS.md for what is planned.",
           call. = FALSE)
}

#' @noRd
ddm_refuse_contaminant <- function(contaminant, what) {
  # anything but TRUE or FALSE is a typo, not a request for the feature
  if (!is.logical(contaminant) || length(contaminant) != 1L ||
      is.na(contaminant)) {
    frm_stop(what, "(): `contaminant` must be FALSE; contaminant = TRUE ",
             "is not built for this family.", call. = FALSE)
  }
  if (contaminant) ddm_stop_not_built(what, "contaminant = TRUE")
  invisible(NULL)
}

#' @noRd
ddm_refuse_cens_slot <- function(what) {
  force(what)
  function(q, dpars, aterms) {
    ddm_stop_not_built(what, "censoring with cens() or trunc()")
  }
}

#' @noRd
ddm_refuse_cens_aterms <- function(aterms, what) {
  if (!is.null(aterms[["cens"]]) || !is.null(aterms[["trunc_lb"]]) ||
      !is.null(aterms[["trunc_ub"]])) {
    ddm_stop_not_built(what, "censoring with cens() or trunc()")
  }
  invisible(NULL)
}
