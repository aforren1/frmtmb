# Checks that need the finished fit.
#
# `family_finalize()` runs at frame assembly, before any parameter has a
# value, so a family that has to look at where the optimizer LANDED has
# had nowhere to do it. That gap is not cosmetic: a family whose
# likelihood is floored in some region of the parameter space converges
# there without a warning, `logLik()` reads `object$opt$objective`
# directly (R/methods-fit.R), and no extension can gate either that or
# `AIC()`. This is the hook that closes it.

#' Run every fit-end check on a finished `frmtmb_fit`.
#'
#' Called once from `fit_assembled()`, after `check_convergence()` and
#' before the object is returned, so a check sees the estimates, the
#' frame and the objective together.
#'
#' @noRd
fit_end_checks <- function(fit) {
  # Every check here runs AFTER frm() has done all the work, so a check
  # that throws must not take the fit with it. A third-party family's
  # hook is arbitrary code from outside this package, and even core's
  # own coverage report evaluates a nonlinear body to find its
  # arguments. A fit that finished is worth more than a diagnostic that
  # did not: the failure is reported and the fit is returned.
  tryCatch(ps_coverage_warning(fit), error = keep_missing(function(e) {
    frm_warning("The ps() knot-span coverage report failed after the fit ",
                "finished, so nothing was checked about it: ",
                conditionMessage(e),
                ". The fit itself is complete and unaffected", call. = FALSE)
  }))
  tryCatch(discrete_loglik_warning(fit), error = keep_missing(function(e) {
    frm_warning("The check of the discrete log-likelihood failed after the ",
                "fit finished, so nothing was checked about it: ",
                conditionMessage(e),
                ". The fit itself is complete and unaffected", call. = FALSE)
  }))
  for (resp in fit$spec$responses) {
    fc <- resp$family[["post"]][["fit_check"]]
    if (!is.function(fc)) next
    fam_nm <- resp$family[["family"]]
    # A family's own verdict on where its fit landed (a threshold no row
    # places, an unidentified disc intercept or coi, a degenerate
    # mixture component) names the parameters and the remedy, so the
    # standard-error warning that would follow it is not raised.
    explained <- function(w) {
      if (is.environment(fit$cache)) fit$cache$se_explained <- fam_nm
    }
    tryCatch(withCallingHandlers(fc(fit, resp$resp_name),
                                 warning = explained), error = function(e) {
      frm_warning("The '", fam_nm, "' family's post$fit_check hook failed ",
                  "after the fit finished, so whatever it checks was not ",
                  "checked: ", conditionMessage(e),
                  ". The fit itself is complete and unaffected",
                  call. = FALSE)
    })
  }
  invisible(NULL)
}

#' A warning for a fit of discrete responses with a log-likelihood above 0.
#'
#' A probability mass is at most 1, so that value cannot be a likelihood. It
#' marks a density that lost its digits (`beta_binomial()` read +94208 at code
#' 0 before its density was rewritten) or a custom family that dropped a
#' constant. It is one-sided: a value below 0 can be just as wrong, so this does
#' not vouch for a fit. Under a prior `logLik()` holds the log prior as well, and
#' REML integrates over a flat prior, so neither is a mass there.
#'
#' @noRd
discrete_loglik_warning <- function(fit) {
  if (isTRUE(fit$REML) || NROW(fit$prior) > 0L) return(invisible(NULL))
  ll <- -fit$opt$objective
  if (length(ll) != 1L || !is.finite(ll)) return(invisible(NULL))
  # a sum of near-zero log masses can pass 0 by its rounding
  if (ll <= 1e-8 * max(1, fit$frame[["n_obs"]])) return(invisible(NULL))
  types <- vapply(fit$spec$responses, function(resp) {
    as.character(resp$family[["type"]] %||% "")[1L]
  }, "")
  if (!all(types == "discrete")) return(invisible(NULL))
  frm_warning(
    "The log-likelihood is ", format(ll, digits = 7), ", above 0, but every ",
    "response is discrete and a probability mass is at most 1. The density ",
    "has lost its digits at these estimates, or a custom family has dropped ",
    "a constant. Do not read this fit",
    call. = FALSE, class = "frmtmb_discrete_loglik")
  invisible(NULL)
}
