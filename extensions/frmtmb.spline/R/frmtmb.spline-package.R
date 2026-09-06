#' @keywords internal
#'
#' @section What this package reads that frmtmb does not promise:
#' Nothing, since frmtmb 0.52.0. Everything read off a fitted object is
#' now a documented seam.
#'
#' * `fit$estimates` and `fit$obj` both have precedent: the
#'   `frmtmb::frmtmb-extension-api` example reads `fit$estimates`, and
#'   `frmtmb.sample::frm_sample()` reads `fit$obj`.
#' * The joint covariance of a grid prediction comes from
#'   [frmtmb::frm_lp_basis()], which is exported and documented. Up to
#'   frmtmb 0.51.0 there was no such seam and `frm_curve()` read
#'   `fit$cache$Vjoint`, the memo written by an unexported
#'   `get_joint_cov()`. That reach is gone, and with it the whole
#'   section this one replaces.
#'
#' The covariance check survives the change and still runs on every
#' call: what it now verifies is that this package reads the seam
#' correctly, rather than that a reconstruction reproduced core's
#' number. [frm_curve()], section "The route to the covariance", sets it
#' out.
#'
#' @section What core has since supplied, and what is still missing:
#' Three things this package used to work around are seams in frmtmb
#' 0.52.0.
#'
#' * `lccdf`. A right-censored row is scored from `log S` directly, in
#'   closed form on all three [royston_parmar()] scales, so the floor at
#'   -35.127363 and the flat region past `-log S = 30` are both gone.
#' * `post$fit_check`, the fit-end family hook. A non-monotone fit warns
#'   as it is returned rather than only when someone calls
#'   [rp_floored()].
#' * [frmtmb::frm_lp_basis()], above.
#'
#' What is still missing is one thing, and it is the protocol's rather
#' than this family's. `frmtmb_structure(loglik =)` returns one AD
#' scalar, so core never sees the individual factors of a likelihood
#' that factorizes. A **per-row log-likelihood slot** would give
#' `loo()` and `waic()` the pointwise matrix they need, and a
#' **per-group log-likelihood slot** would give `frm(importance =)` the
#' one value per group it corrects with. Neither exists, so `logLik()`
#' and `AIC()` on a [royston_parmar()] fit still report whatever the
#' optimizer reached with no way to say which rows carried it;
#' [rp_floored()] is where that is read instead.
#'
#' @section One limit that is core's to fix:
#' * A mapped random-effect block is untested here, because
#'   `frmtmb::frmtmb_control()` takes no `map` argument and there is
#'   therefore no supported route to one. The nearest reachable
#'   analogue is a distributional parameter held fixed, which sets
#'   `betad_fixed_idx` and takes the same index-remapping path:
#'   `bf(y ~ s(x, k = 8), sigma = 0.5)` works at a measured
#'   `cov_rel_error` of 1.55e-15.
"_PACKAGE"

# frmtmb is a Depends, so that one library(frmtmb.spline) call gives a
# user the formula grammar, the frm() a curve is read off, and the
# family constructor. The seams this package builds on are imported by
# name as well, because a namespace that is loaded and not attached
# reaches nothing through the search path.
#' @importFrom frmtmb custom_family frmtmb_register_compat
#'   compat_rule_builder
#' @importFrom stats predict qnorm quantile rnorm setNames
NULL
