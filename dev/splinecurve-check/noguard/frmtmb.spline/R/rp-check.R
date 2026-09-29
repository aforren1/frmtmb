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
#' is read off with `frm_linpred(type = "link", dpar = )` and the basis is
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
    as.numeric(frmtmb::frm_linpred(object, type = "link", dpar = p))
  })
  x <- log(y)
  eta <- sp_rp_eta(sp_rp_basis(kn, x), dp)
  detadx <- sp_rp_eta(sp_rp_dbasis(kn, x), dp)
  av <- object$frame[["aterm_values"]][[rnm]]
  cens <- av[["cens"]]
  if (is.null(cens)) cens <- rep(0, length(y))
  # An interval row's likelihood reads the survival function at BOTH of
  # its ends, so a spline that turns over at the upper end is as much a
  # non-hazard as one that turns over at the lower end.
  detadx_hi <- rep(NA_real_, length(y))
  i_int <- which(cens == 2)
  if (length(i_int) && !is.null(av[["cens_y2"]])) {
    x2 <- log(as.numeric(av[["cens_y2"]])[i_int])
    detadx_hi[i_int] <- sp_rp_eta(sp_rp_dbasis(kn, x2),
                                  lapply(dp, function(v) v[i_int]))
  }
  scale <- environment(fam[["lpdf"]])$cfg$scale
  # -log S, the quantity whose size governs the accuracy of the scored
  # censored term on every scale
  nlogS <- switch(scale,
    hazard = exp(eta),
    odds = log1p(exp(eta)),
    normal = -stats::pnorm(-eta, log.p = TRUE))
  x_hi <- rep(NA_real_, length(y))
  if (length(i_int) && !is.null(av[["cens_y2"]])) x_hi[i_int] <- x2
  list(eta = eta, detadx = detadx, detadx_hi = detadx_hi, nlogS = nlogS,
       cens = as.numeric(cens), scale = scale, n = length(y), x = x,
       x_hi = x_hi, dp = dp, knots = kn)
}

#' How far the fitted survival rises, over the given rows.
#'
#' For each row, S is evaluated with that row's own coefficients on a
#' grid from the smallest observed log time to the row's own (upper)
#' time, and the rise is the largest amount S climbs above its running
#' minimum there. A user judges a rising survival function by that
#' size, which the row count does not give.
#'
#' @noRd
sp_rp_rise <- function(f, rows, n_grid = 200L) {
  if (!length(rows)) return(0)
  lo <- min(f$x)
  surv <- function(eta) switch(f$scale,
    hazard = exp(-exp(eta)),
    odds = 1 / (1 + exp(eta)),
    normal = stats::pnorm(-eta))
  max(vapply(rows, function(i) {
    hi <- max(f$x[i], f$x_hi[i], na.rm = TRUE)
    if (hi <= lo) return(0)
    g <- seq(lo, hi, length.out = n_grid)
    s <- surv(sp_rp_eta(sp_rp_basis(f$knots, g),
                        lapply(f$dp, function(v) v[i])))
    max(s - cummin(s))
  }, 0))
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
#' @section Censored rows where the survival function rises:
#' The floor above is only ever used on an EVENT row, because only an
#' event row has a density. A censored row contributes `log S`, which
#' the family scores exactly whatever the sign of the derivative. So a
#' censored row with a non-positive `d(eta)/d(log t)` does not make
#' `logLik()` wrong. It makes the MODEL wrong: the fitted survival
#' function rises with time there, and a survival function that rises
#' is not one.
#'
#' This can happen with no event row affected at all. A random effect or
#' a covariate on `gamma1` gives each group its own slope, the
#' `log(gamma1 + u)` barrier that holds a slope positive lives in the
#' density, and a group whose rows are ALL censored contributes no
#' density. Measured on 40 centres of 10 with five centres followed to
#' a common administrative time with no deaths, seeds 20260910 to
#' 20260915: on 4 of 6 seeds those five centres converge at slopes of
#' -0.18 to -0.31, and one centre's fitted survival goes from 4.2e-50 at
#' `t = 1e-12` to 0.774 at `t = 2.8`.
#'
#' `n_nonmonotone_censored` counts those rows. They are counted apart
#' from `n_nonmonotone` because the two answer different questions:
#' `n_nonmonotone` says the reported likelihood is not the model's, and
#' `n_nonmonotone_censored` says the fitted model is not a survival
#' distribution where the data are. An interval-censored row is tested at
#' both ends.
#'
#' `action = "error"` REFUSES on `n_nonmonotone` and WARNS on
#' `n_nonmonotone_censored`, and [frm_curve()] and its two companions do
#' the same, so on a fit with censored rows only they answer with the
#' warning. The warning gives `max_survival_rise`, the largest amount the
#' fitted survival climbs above its running minimum on any flagged row's
#' own coefficients, because the size is what a user needs to judge it.
#' The line is drawn where a negative hazard contradicts an observed
#' event, as flexsurv draws it: `flexsurv::dsurvspline()` sets the
#' density to 0 where the derivative is not positive, so an event row
#' there makes the likelihood `-Inf` and the fit avoids it, while a
#' censored row is never checked and a rise there passes silently.
#' frmtmb draws the same line and says so. brms's `cox()` family cannot
#' produce a rising survival at all, because its baseline hazard is an
#' M-spline basis times non-negative weights; rstpm2 penalizes a
#' negative hazard during the fit.
#'
#' It also fires on an ordinary design, and a user should expect it
#' there: a cure fraction with a time-varying effect. With
#' `gamma1 ~ arm` and part of each arm never having the event, the
#' spline past an arm's last event is fitted to censored rows only, and
#' it can turn over there. Measured on two arms of 250 with 40 and 55
#' percent cured, `df` 3 to 6, 20 seeds each: 27 of 80 fits fire, every
#' flagged row lies past its own arm's last event, and the fitted
#' survival rises across them by 5.2e-04 to 0.13. The count is right,
#' since that fitted survival function does rise, and the rise is where
#' the data carry no event to prevent it. With no time-varying effect
#' this cannot happen past the last event: every row then shares one
#' spline, and past the boundary knot its derivative is the one at the
#' last event row. A fit whose flagged rows are ALL censored warns and
#' does not refuse. A fit with a flagged EVENT row refuses, as in 0.7.0,
#' whatever else it flags, and a cure design fitted with proportional
#' hazards can: its shared spline can turn over at an event row.
#' Fewer knots or a single `gamma1` remove a censored-row rise.
#'
#' The test is at each row's own time. For `df = 1` the derivative is
#' the same at every time, so a row test is exact. For `df >= 2` the
#' derivative can dip between two observed times of one group and
#' recover at both; this check does not search for that.
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
#'   floor was used on an event row, and warns when a censored row sits
#'   where the fitted survival function rises. The deep censored count
#'   never refuses under either value; since frmtmb 0.52.0 it is a
#'   diagnostic. `"report"` returns the same numbers silently.
#' @param max_nlogS The `-log S` on a censored row above which the row
#'   is reported as barely constrained. The default 19.2 is where the
#'   OLD probability-scale arithmetic passed 1e-8 of error; it is kept
#'   as the threshold so that the two versions report the same rows.
#'
#' @return A list with `n_censored_deep`, `max_nlogS`, `threshold`,
#'   `n_nonmonotone`, `n_nonmonotone_censored`, `max_survival_rise`,
#'   `scale` and `n_obs`,
#'   returned invisibly when nothing refuses. The offending row indices
#'   are the `"rows"` attribute, a list with elements `censored`,
#'   `nonmonotone` and `nonmonotone_censored`.
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
  action <- frm_match_arg(action)
  if (!is.numeric(max_nlogS) || length(max_nlogS) != 1L ||
      !is.finite(max_nlogS) || max_nlogS <= 0) {
    frm_stop("`max_nlogS` must be one positive finite number: it is the ",
             "largest -log S on a censored row that this family still scores ",
             "accurately", call. = FALSE)
  }
  fam <- sp_rp_family_of(object)
  if (is.null(fam)) {
    frm_stop("rp_floored() reads the floors of a royston_parmar() fit, and ",
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
  # A censored row is not floored, so it cannot join mono_rows without
  # changing what n_nonmonotone means. It is the only row a group with
  # no events has, so without this count such a group is invisible.
  rise_rows <- which(f$cens != 0 &
                       (f$detadx <= 0 | (!is.na(f$detadx_hi) &
                                           f$detadx_hi <= 0)))
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
              n_nonmonotone_censored = length(rise_rows),
              max_survival_rise = sp_rp_rise(f, rise_rows),
              scale = f$scale,
              n_obs = f$n)
  attr(out, "rows") <- list(censored = cens_rows, nonmonotone = mono_rows,
                            nonmonotone_censored = rise_rows)
  if (identical(action, "report")) return(out)
  # The deep censored rows are scored exactly through the family's
  # lccdf slot, so their count is a diagnostic. An EVENT row on a
  # non-positive slope makes logLik() a pseudo-likelihood, which
  # refuses. A CENSORED row there is scored exactly too, and the rise it
  # marks is usually extrapolation past a group's last event: that
  # warns, with its size (user decision, 2026-09-24; flexsurv draws the
  # same line, silently).
  if (length(mono_rows)) frm_stop(sp_rp_refusal(out, f), call. = FALSE)
  if (length(rise_rows)) frm_warning(sp_rp_rise_text(out), call. = FALSE)
  invisible(out)
}

#' The warning for censored rows where the fitted survival rises. One
#' template, shared with the fit-end hook, so the two read the same.
#'
#' @noRd
sp_rp_rise_text <- function(out) {
  paste0("the fitted survival rises by up to ",
         signif(out[["max_survival_rise"]], 2), " across ",
         out[["n_nonmonotone_censored"]], " of ", out[["n_obs"]],
         " censored rows, where d(eta)/d(log t) is not positive. Those rows ",
         "are scored exactly, so logLik() is the model's, but the fitted ",
         "function is not a survival function there. It is usually past a ",
         "group's last event, where the spline extrapolates: a cure ",
         "fraction with a covariate or random effect on gamma1, or a group ",
         "with no events. Fewer knots, or no slope term for such a group, ",
         "removes it. rp_floored(action = \"report\") names the rows")
}

#' The refusal text. One template, so that both floors read the same way
#' and the message-uniqueness property holds.
#'
#' @noRd
sp_rp_refusal <- function(out, f) {
  parts <- character(0)
  if (out$n_nonmonotone) {
    parts <- c(parts, paste0(
      "this fit's reported likelihood is not the model's. ",
      out$n_nonmonotone, " of ", out$n_obs,
      " observed rows have a non-positive d(eta)/d(log t) at the fitted ",
      "parameters, so no hazard exists there and their true log density ",
      "is -Inf. This family floors them to keep the optimizer alive, ",
      "which makes logLik() and AIC() a pseudo-likelihood rather than a ",
      "density"))
  }
  paste0("rp_floored(): ", paste(parts, collapse = ". Separately, "),
         ". The row indices are in the \"rows\" attribute of ",
         "rp_floored(action = \"report\"). The remedy for a ",
         "non-monotone spline is fewer knots. The deep censored count is ",
         "not part of this refusal: since frmtmb 0.52.0 this family ",
         "supplies lccdf and log S is scored exactly, so that count is a ",
         "diagnostic")
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
