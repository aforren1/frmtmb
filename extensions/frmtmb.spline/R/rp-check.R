# The post-fit check that stands in for a core `lccdf` slot.
#
# frmtmb forms a right-censored contribution on the PROBABILITY scale,
# `log(Fub - F(y))` (`R/objective.R:100`), which without truncation is
# `log(1 - F)`. A survival family therefore cannot hand back an accurate
# `log S` once `1 - F` stops being representable, and no amount of care
# inside the family changes that: the number core asks for is `F`, and
# the complement of a double near 1 carries absolute error about
# `.Machine$double.eps`.
#
# So the error in the scored `log S` is about `eps / S`, and that single
# expression governs all three scales. Measured, by forming `F` for a
# given `-log S` on each scale and reading `log(1 - F)` back:
#
#   -log S    computed        abs error    eps / S
#   10        -10             1.3e-13      4.9e-12
#   15        -15             9.0e-11      7.3e-10
#   19.2      -19.2           2.4e-10      4.8e-08
#   25        -25             3.8e-06      1.6e-05
#   30        -29.99983       1.7e-04      2.4e-03
#   36        -34.94504       1.05         9.6e-01
#   40        -35.12736       4.87         5.2e+01
#
# The three scales agree to every printed digit at every row, because
# they differ only in how `eta` maps to `S`. The reviewer's two
# thresholds are the same threshold said twice: `eta = 6` on the normal
# scale is `-log S = 20.74`, and `-log S = 19.2` is `eta = 5.745`.
#
# Past `-log S = 30` the value is not merely inaccurate, it is FLAT: the
# gradient of the scored term with respect to the coefficients is exactly
# zero, so the optimizer prices such a row at a constant and fits the
# others as if it were free. That is what produces a converged,
# warning-free fit whose reported log likelihood is wrong by thousands
# and whose treatment coefficient is wrong by tens of percent.
#
# The fix belongs in core (`dev/spline-seam-proposal.md`, Part 1, the
# `lccdf` slot). Until it lands, this package REFUSES rather than
# floors: a fit in that region is not one this family may report numbers
# for.

#' The finalized `royston_parmar()` family of a fit, or `NULL`.
#'
#' @noRd
sp_rp_family_of <- function(object) {
  if (!inherits(object, "frmtmb_fit")) return(NULL)
  fam <- try(stats::family(object), silent = TRUE)
  if (inherits(fam, "try-error") || is.null(fam)) return(NULL)
  if (!identical(fam[["family"]], "royston_parmar")) return(NULL)
  fam
}

#' The knot vector the fit's family actually carries.
#'
#' `family_finalize()` rebuilds the family with its knots baked into the
#' densities' enclosing environment, so this is where they live. It is
#' this package's own closure, not a reach into core.
#'
#' @noRd
sp_rp_knots_of <- function(fam) environment(fam[["lpdf"]])$allknots

#' The fitted spline, evaluated at the observed rows.
#'
#' Every spline coefficient is a distributional parameter, so each one
#' is read off with `predict(type = "link", dpar = )` and the basis is
#' rebuilt at the observed log times. This is the same arithmetic the
#' objective does, on the same knots.
#'
#' @noRd
sp_rp_fitted <- function(object, fam) {
  rspec <- frmtmb::single_response(object, "rp_floored()")
  rnm <- rspec$resp_name
  y <- as.numeric(object$frame[["y"]][[rnm]])
  kn <- sp_rp_knots_of(fam)
  dp <- lapply(fam[["dpars"]], function(p) {
    as.numeric(stats::predict(object, type = "link", dpar = p))
  })
  x <- log(y)
  eta <- sp_rp_eta(sp_rp_basis(kn, x), dp)
  detadx <- sp_rp_eta(sp_rp_dbasis(kn, x), dp)
  cens <- object$frame[["aterm_values"]][[rnm]][["cens"]]
  if (is.null(cens)) cens <- rep(0, length(y))
  scale <- environment(fam[["lpdf"]])$cfg$scale
  # -log S, the quantity whose size governs the accuracy of the scored
  # censored term on every scale
  nlogS <- switch(scale,
    hazard = exp(eta),
    odds = log1p(exp(eta)),
    normal = -stats::pnorm(-eta, log.p = TRUE))
  list(eta = eta, detadx = detadx, nlogS = nlogS, cens = as.numeric(cens),
       scale = scale, n = length(y))
}

#' Rows this family answered with a floor rather than with a density
#'
#' Two things in `royston_parmar()` are floors rather than answers, and
#' both are silent in the fitted object: `logLik()` and `AIC()` report
#' the floored value with nothing to say it is floored. This function is
#' where that goes to be read, and by default it REFUSES rather than
#' reports, because a fit in either region is one whose numbers are not
#' the model's.
#'
#' @section The censored-row floor:
#' frmtmb forms a right-censored contribution as `log(1 - F(y))` on the
#' probability scale (`R/objective.R:100`), so the scored `log S` carries
#' absolute error about `.Machine$double.eps / S`. Past `-log S` of about
#' 19.2 that error passes 1e-8; past 30 the term is FLAT, its gradient
#' exactly zero, and the optimizer prices the row at a constant.
#'
#' The size of the error is a property of the data rather than of the
#' family: a floored row contributes -35.127363 instead of its own
#' `-log S`, so the reported log likelihood is short by about
#' `-log S - 35` per floored row. Two runs of one 600-subject design
#' differing only in seed give 2.4e+03 and 2.166e+04, both converged
#' without a warning and both with the treatment coefficient out by tens
#' of percent.
#'
#' The quantity checked is `-log S` at the fitted parameters, on every
#' censored row. It is one quantity for all three scales: `exp(eta)` on
#' `"hazard"`, `log1p(exp(eta))` on `"odds"` and `-log(Phi(-eta))` on
#' `"normal"`. On the hazard scale it is the cumulative hazard `H`.
#'
#' The real fix is a complementary log-CDF slot in core, so that a family
#' can hand back `log S` instead of `F`. See `dev/spline-seam-proposal.md`
#' in the package sources.
#'
#' @section The monotonicity floor:
#' The cumulative hazard has to increase, so the spline's derivative in
#' log time has to stay positive; nothing enforces it and flexsurv does
#' not enforce it either. Where it goes non-positive there is no hazard
#' and the true log density is `-Inf`, and this family replaces it with a
#' large finite number so that the optimizer has something to work with.
#' That keeps the fit alive and makes `logLik()` a pseudo-likelihood: a
#' 60 percent cure-fraction dataset has been measured converging with 6
#' such rows and a reported log likelihood 3952 units away from the
#' density's.
#'
#' @section What this cannot do:
#' The refusal is POST-FIT. `logLik()` reads `object$opt$objective`
#' directly (`R/methods-fit.R:233-240`) and the family protocol has no
#' hook that runs when a fit finishes, so nothing in this package can
#' make `logLik()` or `AIC()` refuse on their own. The optimizer may
#' therefore have walked through, or stopped inside, the flat region
#' before this function is ever called. Call it on every
#' `royston_parmar()` fit whose data carry censoring; `frm_curve()` and
#' its two companions call it for you.
#'
#' @param object A `frmtmb_fit` with a [royston_parmar()] family.
#' @param action `"error"`, the default, refuses when either floor was
#'   used. `"report"` returns the same numbers without refusing.
#' @param max_nlogS The largest `-log S` on a censored row that is still
#'   scored accurately. The default 19.2 is where `eps / S` passes 1e-8.
#'
#' @return A list with `n_censored_floored`, `max_nlogS`, `threshold`,
#'   `n_nonmonotone`, `scale` and `n_obs`, returned invisibly when
#'   nothing was floored. The offending row indices are the `"rows"`
#'   attribute, a list with elements `censored` and `nonmonotone`.
#'
#' @seealso [royston_parmar()]
#' @examples
#' set.seed(1)
#' n <- 300
#' dd <- data.frame(trt = rep(0:1, each = n / 2))
#' dd$t <- rweibull(n, shape = 1.4, scale = exp(1 - 0.5 * dd$trt))
#' dd$censored <- as.integer(dd$t > 3)
#' dd$t <- pmin(dd$t, 3)
#' fit <- frmtmb::frm(frmtmb::bf(t | cens(censored) ~ trt),
#'                    family = royston_parmar(df = 2), data = dd)
#' str(rp_floored(fit, action = "report"))
#' @export
rp_floored <- function(object, action = c("error", "report"),
                       max_nlogS = 19.2) {
  action <- match.arg(action)
  if (!is.numeric(max_nlogS) || length(max_nlogS) != 1L ||
      !is.finite(max_nlogS) || max_nlogS <= 0) {
    stop("`max_nlogS` must be one positive finite number: it is the ",
         "largest -log S on a censored row that this family still scores ",
         "accurately", call. = FALSE)
  }
  fam <- sp_rp_family_of(object)
  if (is.null(fam)) {
    stop("rp_floored() reads the floors of a royston_parmar() fit, and ",
         "this object's family is '",
         if (inherits(object, "frmtmb_fit")) {
           stats::family(object)[["family"]]
         } else {
           paste0("not a fit at all (", class(object)[1L], ")")
         },
         "'. Nothing else in this package floors anything",
         call. = FALSE)
  }
  f <- sp_rp_fitted(object, fam)
  cens_rows <- which(f$cens != 0 & f$nlogS > max_nlogS)
  mono_rows <- which(f$cens == 0 & f$detadx <= 0)
  mx <- if (any(f$cens != 0)) max(f$nlogS[f$cens != 0]) else 0
  out <- list(n_censored_floored = length(cens_rows),
              max_nlogS = mx,
              threshold = max_nlogS,
              n_nonmonotone = length(mono_rows),
              scale = f$scale,
              n_obs = f$n)
  attr(out, "rows") <- list(censored = cens_rows, nonmonotone = mono_rows)
  if (identical(action, "report")) return(out)
  if (length(cens_rows) || length(mono_rows)) {
    stop(sp_rp_refusal(out, f), call. = FALSE)
  }
  invisible(out)
}

#' The refusal text. One template, so that both floors read the same way
#' and the message-uniqueness property holds.
#'
#' @noRd
sp_rp_refusal <- function(out, f) {
  qty <- switch(out$scale, hazard = "the cumulative hazard H",
                odds = "log(1 + exp(eta))", "-log(Phi(-eta))")
  parts <- character(0)
  if (out$n_censored_floored) {
    parts <- c(parts, paste0(
      out$n_censored_floored, " of ", sum(f$cens != 0),
      " censored rows are scored past the accurate region: -log S, which ",
      "on the ", out$scale, " scale is ", qty, ", reaches ",
      format(out$max_nlogS, digits = 6), " where this family stays ",
      "accurate only to ", format(out$threshold),
      ". frmtmb forms a right-censored term as log(1 - F) on the ",
      "probability scale and core has no complementary log-CDF (lccdf) ",
      "slot a family could use instead, so the scored log S is floored ",
      "at -35.127363 and its gradient is exactly zero past -log S of 30"))
  }
  if (out$n_nonmonotone) {
    parts <- c(parts, paste0(
      out$n_nonmonotone, " of ", out$n_obs,
      " observed rows have a non-positive d(eta)/d(log t) at the fitted ",
      "parameters, so no hazard exists there and their true log density ",
      "is -Inf. This family floors them to keep the optimizer alive, ",
      "which makes logLik() and AIC() a pseudo-likelihood rather than a ",
      "density"))
  }
  paste0("rp_floored(): this fit's reported likelihood is not the ",
         "model's. ", paste(parts, collapse = ". Separately, "),
         ". The row indices are in the \"rows\" attribute of ",
         "rp_floored(action = \"report\"). The remedy for the censored ",
         "term is the lccdf slot proposed in ",
         "dev/spline-seam-proposal.md; until it exists, shorten the ",
         "censoring horizon or refit without the offending rows. The ",
         "remedy for a non-monotone spline is fewer knots. Note that ",
         "this check is POST-FIT: logLik() reads the optimizer's own ",
         "value and no family hook runs when a fit finishes, so nothing ",
         "here could have refused earlier")
}

#' Refuse before reporting a curve off a royston_parmar fit whose
#' likelihood was floored.
#'
#' The curve functions are the documented way to inspect this family, so
#' they are the entry points this package owns that a user reaches after
#' fitting. A fit whose likelihood is a floor artifact is one whose
#' fitted curve is too.
#'
#' @noRd
sp_rp_gate <- function(object) {
  fam <- sp_rp_family_of(object)
  if (is.null(fam)) return(invisible(NULL))
  rp_floored(object, action = "error")
  invisible(NULL)
}
