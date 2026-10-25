#' frmtmb.survival: survival families for frmtmb
#'
#' Two survival families for [frmtmb::frm()].
#'
#' * [cox()] is brms's proportional hazards model with an M-spline
#'   baseline hazard. [cox_baseline()] returns the fitted baseline
#'   weights.
#' * [royston_parmar()] is the flexible parametric model of Royston and
#'   Parmar (2002), with the log cumulative hazard a natural cubic
#'   spline in log time, parameterized as `flexsurv::flexsurvspline()`
#'   parameterizes it. [rp_floored()] reports where a fit's likelihood
#'   is a floor rather than a density.
#'
#' Both moved here from elsewhere: `cox()` from frmtmb itself and
#' `royston_parmar()` from frmtmb.spline. Load this package with
#' `library(frmtmb.survival)`; it attaches frmtmb as well. Without it,
#' `family = "cox"` and `brmsfamily("cox")` stop with a message that
#' names this package.
#'
#' @section What this package reads from frmtmb:
#' Only exported functions. The families are built with
#' [frmtmb::frmtmb_family()] and [frmtmb::custom_family()]. The `cox()`
#' baseline is a simplex declared through `simplex_pars =`, which is
#' what makes frmtmb chart it, search it off a face, give it brms's
#' `sbhaz` prior class and draw it on frmtmb.sample's sampler, and its
#' spline bases are data built through `aterm_data =`. brms's `bhaz()`
#' addition term, which sets that basis on the response side, is
#' registered with [frmtmb::frmtmb_register_aterm()] as an option of the
#' family (`option =`). The Royston and Parmar refusal of a floored fit
#' reaches [frmtmb::frm_curve()] and its two companions through the
#' family's `post$curve_check` hook.
#' `test-bracket-access.R` checks that no function here reads a hazard
#' container, and `test-no-reach.R` that no source file reaches into
#' frmtmb's namespace.
#'
#' @references
#' Royston, P. and Parmar, M. K. B. (2002). Flexible parametric
#' proportional-hazards and proportional-odds models for censored
#' survival data. *Statistics in Medicine* 21, 2175-2197.
#' @keywords internal
"_PACKAGE"

# frmtmb is a Depends, so that one library(frmtmb.survival) call gives a
# user the formula grammar and frm() as well as the families. The seams
# this package builds on are imported by name too, because a namespace
# that is loaded and not attached reaches nothing through the search
# path.
#' @importFrom frmtmb custom_family frmtmb_family frmtmb_register_compat
#' @importFrom frmtmb frmtmb_register_aterm
#' @importFrom frmtmb frm_match_arg frm_message frm_stop frm_warning
#'   compat_rule_builder compat_aterm_rules mo_simplex single_response
#' @importFrom stats predict qnorm quantile rnorm setNames
NULL
