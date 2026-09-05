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
  tryCatch(ps_coverage_warning(fit), error = function(e) {
    warning("The ps() knot-span coverage report failed after the fit ",
            "finished, so nothing was checked about it: ",
            conditionMessage(e),
            ". The fit itself is complete and unaffected", call. = FALSE)
  })
  for (resp in fit$spec$responses) {
    fc <- resp$family[["post"]][["fit_check"]]
    if (!is.function(fc)) next
    fam_nm <- resp$family[["family"]]
    tryCatch(fc(fit, resp$resp_name), error = function(e) {
      warning("The '", fam_nm, "' family's post$fit_check hook failed ",
              "after the fit finished, so whatever it checks was not ",
              "checked: ", conditionMessage(e),
              ". The fit itself is complete and unaffected",
              call. = FALSE)
    })
  }
  invisible(NULL)
}
