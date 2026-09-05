#' @keywords internal
#'
#' @section What this package reads that frmtmb does not promise:
#' Two of the three things this package reads off a fitted object are
#' documented seams and one is not. It is named here so that a reader
#' does not have to discover it from the sources.
#'
#' * `fit$estimates` and `fit$obj` both have precedent: the
#'   `frmtmb::frmtmb-extension-api` example reads `fit$estimates`, and
#'   `frmtmb.sample::frm_sample()` reads `fit$obj`.
#' * **`fit$cache$Vjoint` has no precedent and is an internal.** It is
#'   the memo written by frmtmb's `get_joint_cov()`, which is not
#'   exported and not documented, and neither the slot nor its
#'   `list(V =, names =)` shape appears anywhere in frmtmb's own
#'   documentation. `frm_curve()` reads it because the alternative,
#'   recomputing the joint covariance here, was both far slower and
#'   WRONG on an autoscaled fit.
#'
#' The consequences of that reach are bounded and are set out in full
#' under [frm_curve()], section "The one internal this reaches into": a
#' change to the name or the shape costs speed, a change to the meaning
#' is caught by the covariance check that every call makes, and absence
#' is the ordinary case on a fresh fit and is handled by warming the
#' cache first. None of the three can produce a wrong number.
#'
#' `dev/spline-seam-proposal.md` Part 1a asks core for an exported
#' accessor, which would retire the reach entirely.
#'
#' @section Two limits that are core's to fix:
#' * `logLik()` and `AIC()` on a [royston_parmar()] fit report whatever
#'   the optimizer reached, and neither can be gated from an extension:
#'   `logLik()` reads `object$opt$objective` and the family protocol has
#'   no hook that runs when a fit finishes. Where the likelihood was
#'   floored, those two numbers are wrong and say nothing.
#'   [rp_floored()] refuses, and [frm_curve()] calls it, but a user who
#'   reads `AIC()` and nothing else gets no signal. The fix is core's:
#'   an `lccdf` slot, or a post-fit family hook.
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
