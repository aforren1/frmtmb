# The post-fit report on the two floors this family used to have.
#
# Up to frmtmb 0.51.0 core formed a right-censored contribution on the
# PROBABILITY scale, `log(Fub - F(y))` (`R/objective.R:100`), which
# without truncation is `log(1 - F)`. A survival family could not hand
# back an accurate `log S` once `1 - F` stopped being representable: the
# scored error was about `eps / S`, the term was FLAT past `-log S` of
# 30 with a gradient of exactly zero, and a converged, warning-free fit
# could report a log likelihood wrong by tens of thousands.
#
# frmtmb 0.52.0 added the `lccdf` slot the proposal asked for, and
# `royston_parmar()` supplies it (`R/royston-parmar.R`). A right-censored
# row is now scored from `log S` directly, in closed form on all three
# scales, so THE CENSORED FLOOR NO LONGER EXISTS. Measured on the
# hazard scale, `-log S` at 40 was scored as -35.127363 before and is
# scored as -40 exactly now, and its gradient in `eta` is right rather
# than zero.
#
# So the censored count below is a DIAGNOSTIC rather than a refusal.
# `-log S` past about 19.2 still says something worth hearing: a fit
# with a censored row whose survival probability is `exp(-40)` is a fit
# whose data barely constrain that row, whatever the arithmetic does.
# It is reported, and it does not stop anything.
#
# The monotonicity floor is UNCHANGED and still refuses. The cumulative
# hazard has to increase, nothing enforces it, and where the spline's
# derivative in log time goes non-positive there is no hazard and the
# true log density is `-Inf`. `sp_floor_pos()` replaces it with a large
# finite number, which keeps the optimizer alive and makes `logLik()` a
# pseudo-likelihood. No core seam addresses that; it is a property of
# the model, not of the arithmetic.
#
# What this still cannot do: `logLik()` reads `object$opt$objective`
# directly, so this is post-fit either way. frmtmb 0.52.0 does now run a
# family's `post$fit_check` when a fit finishes, which is how the
# monotonicity report reaches a user who never calls this function.

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

#' Report the deep censored rows and the non-monotone rows of a fit
#'
#' Two things in `royston_parmar()` used to be floors rather than
#' answers. One of them is gone; the other refuses.
#'
#' @section The censored rows: a report, no longer a floor:
#' Up to frmtmb 0.51.0 core formed a right-censored contribution as
#' `log(1 - F(y))` on the probability scale, so the scored `log S`
#' carried absolute error about `.Machine$double.eps / S`: past `-log S`
#' of about 19.2 that error passed 1e-8, and past 30 the term was FLAT
#' with a gradient of exactly zero.
#'
#' frmtmb 0.52.0 added the `lccdf` slot and this family supplies it, in
#' closed form on all three scales, so a right-censored row is scored
#' from `log S` directly and no complement is formed. Measured on the
#' hazard scale, `-log S = 40` was scored as -35.127363 and is now
#' scored as -40 exactly.
#'
#' The count is therefore a DIAGNOSTIC and never refuses. It still says
#' something: a censored row whose fitted survival probability is
#' `exp(-40)` is one the data barely constrain, whatever the arithmetic
#' does. The quantity is `-log S` at the fitted parameters, and it is
#' one quantity for all three scales: `exp(eta)` on `"hazard"`,
#' `log1p(exp(eta))` on `"odds"` and `-log(Phi(-eta))` on `"normal"`.
#' On the hazard scale it is the cumulative hazard `H`.
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
#' `logLik()` reads `object$opt$objective` directly, so a check that
#' runs after the fit cannot make `logLik()` or `AIC()` refuse on their
#' own. What frmtmb 0.52.0 does provide is a fit-end hook: this family
#' declares `post$fit_check`, so a fit with a non-monotone row warns as
#' it is returned rather than only when someone calls this function.
#' `frm_curve()` and its two companions still call it for you.
#'
#' @param object A `frmtmb_fit` with a [royston_parmar()] family.
#' @param action `"error"`, the default, refuses when the MONOTONICITY
#'   floor was used. The censored count never refuses under either
#'   value; since frmtmb 0.52.0 it is a diagnostic. `"report"` returns
#'   the same numbers without refusing.
#' @param max_nlogS The `-log S` on a censored row above which the row
#'   is reported as barely constrained. The default 19.2 is where the
#'   OLD probability-scale arithmetic passed 1e-8 of error; it is kept
#'   as the threshold so that the two versions report the same rows.
#'
#' @return A list with `n_censored_deep`, `max_nlogS`, `threshold`,
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
  # named for what it counts. It was `n_censored_floored` while the
  # censored term WAS floored; with lccdf exact, a row past the
  # threshold is deep rather than floored, and a field that says
  # otherwise is the same defect as a help page that contradicts the
  # code.
  out <- list(n_censored_deep = length(cens_rows),
              max_nlogS = mx,
              threshold = max_nlogS,
              n_nonmonotone = length(mono_rows),
              scale = f$scale,
              n_obs = f$n)
  attr(out, "rows") <- list(censored = cens_rows, nonmonotone = mono_rows)
  if (identical(action, "report")) return(out)
  # Only the monotonicity floor refuses now. The censored rows are
  # scored exactly through the family's lccdf slot, so their count is a
  # diagnostic and stopping on it would refuse a correct fit.
  if (length(mono_rows)) stop(sp_rp_refusal(out, f), call. = FALSE)
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
         "rp_floored(action = \"report\"). The remedy for a ",
         "non-monotone spline is fewer knots. The censored term is no ",
         "longer part of this refusal: since frmtmb 0.52.0 this family ",
         "supplies lccdf and log S is scored exactly, so the censored ",
         "count that rp_floored() also returns is a diagnostic")
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
