# What this package tells frmtmb about itself at load time.
#
# The registration uses a seam frmtmb exports for the purpose
# (?frmtmb::`frmtmb-sampling-api`, section "The compatibility
# registry"). Registering from .onLoad() rather than at top level is
# what a contributor outside the package must do: by then every
# namespace is sealed, so the collation-order question that governs
# frmtmb's own in-package contributors does not arise.

#' @noRd
.onLoad <- function(libname, pkgname) {
  frmtmb_register_compat(features = c(wiener = "family"),
                         rules = ddm_compat_rules)
  invisible()
}

#' The compatibility rules that name the `wiener` family.
#'
#' Every row here was run, not reasoned about. `untested` is used where
#' the pair was not exercised, which is the honest third state: an
#' absent guard and a passing guard look the same from outside.
#'
#' @noRd
ddm_compat_rules <- function() {
  b <- compat_rule_builder()
  r <- b$r
  r("wiener", "vint()", "works",
    "Required, not merely supported: the boundary a trial ended at is data and reaches the density through vint(). Omitting it is refused, because the density would otherwise read a NULL and the log likelihood would silently collapse to zero terms.")
  r("wiener", "cens()", "refused",
    "The family declares no lcdf. The Wiener first-passage distribution function is a third series with its own truncation problem and none of it is written here.")
  r("wiener", "trunc()", "refused",
    "Same reason as cens(): no lcdf, so frmtmb has no normalizing constant to divide by.")
  r("wiener", "weights()", "works",
    "Verified: the weighted log likelihood is the unweighted one at unit weights and scales as it should.")
  r("wiener", "simulate", "works",
    "Verified against the fitted parameters. Draws are conditional on each row's boundary, by inverse transform through RWiener's defective quantile function; without RWiener a discretized forward simulation stands in.")
  r("wiener", "fitted", "works",
    "The mean is the conditional mean response time for the row's own boundary, in closed form.")
  r("wiener", "predict", "works",
    "vint() is mandatory on newdata, so the boundary must be supplied there as well.")
  r("wiener", "residuals", "conditional",
    "type = \"response\" works. \"pearson\" and \"deviance\" are both refused: the family declares no variance function and no unit deviance, because the conditional variance of a first-passage time was not worth writing for a residual nobody reads on response times.")
  r("wiener", "residuals_osa", "untested",
    "One-step-ahead residuals re-tape the objective with the response promoted to a parameter. Nothing here exercises that path.")
  r("wiener", "REML", "untested",
    "The drift rate is the primary dpar and would be integrated out. Not exercised.")
  b$rules()
}
