# What this package tells frmtmb about itself at load time.
#
# The registration uses the seam frmtmb exports for the purpose
# (?frmtmb::`frmtmb-sampling-api`, section "The compatibility registry").
# Registering from .onLoad() rather than at top level is what a
# contributor outside the package must do: by then every namespace is
# sealed, so the collation-order question that governs frmtmb's own
# in-package contributors does not arise.
#
# The rules below name a feature of THIS package on one side and a
# feature of core frmtmb on the other, and every core name used here is
# one that exists on the base this package was written against
# (frmtmb 0.50.0, commit 5dfdd84): cens(), trunc(), weights(), s(),
# t2(), gp, smooth, rr, mixture, mvbf, nl, REML, quadrature, fitted,
# predict, simulate, residuals and residuals_osa are all in
# frm_compat_features() there, spelled as the `name` column spells them
# (a callable feature carries its parentheses; a covariance structure
# and a method do not). Only the two names this package brings are new,
# and both are declared through `features =` rather than assumed.
#
# A sibling lane is changing frmtmb_register_compat() this round so that
# an unknown feature is refused at registration rather than accepted and
# never resolved. Nothing here should be affected: every name on the
# core side of a rule below resolves on 5dfdd84 today, which was checked
# against frm_compat_features() rather than remembered.

#' @noRd
.onLoad <- function(libname, pkgname) {
  frmtmb_register_compat(
    features = c(royston_parmar = "family", frm_curve = "method",
                 rp_floored = "method"),
    rules = sp_compat_rules)
  invisible()
}

#' The compatibility rules for what this package supplies.
#'
#' Every row was run. `untested` is used where the pair was not
#' exercised, which is the honest third state: an absent guard and a
#' passing guard look the same from outside.
#'
#' @noRd
sp_compat_rules <- function() {
  b <- compat_rule_builder()
  r <- b$r

  ## ---- the family -------------------------------------------------
  r("royston_parmar", "cens()", "works",
    "Exact since frmtmb 0.52.0. Core used to form a right-censored term as log(1 - F(y)) on the probability scale and offered a family no complementary log-CDF slot, so the scored log S was floored at -35.127363 and its GRADIENT was exactly zero past -log S of 30: a converged, warning-free fit could report a log likelihood short by tens of thousands and a treatment coefficient out by tens of percent. Core now takes an lccdf slot and this family supplies it, in closed form on all three scales, so no complement is ever formed. Measured on the hazard scale: -log S of 40 was scored as -35.127363 and is scored as -40 exactly, and -log S of 100 as -100. The identity against flexsurv on bc still reaches only -log S = 2.01, so the evidence for the deep region is the direct arithmetic check rather than that test.")
  r("royston_parmar", "trunc()", "conditional",
    "Unchanged, and the one the new slot does NOT close. A truncation window is divided out as log(Fub - Flb) on the probability scale (R/objective.R), so an upper bound deep in the tail loses accuracy exactly as a censored row used to, and a LEFT truncation bound (delayed entry, which is routine in survival) meets the same representability problem from the other side. lccdf is a log SURVIVOR function and closing truncation needs a windowed log difference. Verified for the case that must change nothing: a bound below every observed time reproduces the untruncated coefficients to 0.05.")
  r("royston_parmar", "weights()", "works",
    "Verified: the weighted log likelihood equals the unweighted one at unit weights, and doubling every weight doubles it.")
  r("royston_parmar", "s()", "works",
    "A smooth on any of the spline coefficients fits. s() on mu is a smooth covariate effect under proportional hazards; s() on gamma1 or higher is a smooth time-varying effect. Both fit and both are readable with frm_curve().")
  r("royston_parmar", "smooth", "works",
    "The random-effect block a penalized smooth becomes, for the same reason. The family imposes nothing on the block structure: it consumes one linear predictor per spline coefficient and does not care how each was built.")
  r("royston_parmar", "predict", "works",
    "type = \"link\" gives the linear predictor of any spline coefficient, which is what frm_curve() reads. type = \"response\" is refused, and correctly: see the fitted row.")
  r("royston_parmar", "fitted", "refused",
    "The family declares no post$mean_fn, because the mean of a Royston-Parmar survival time has no closed form: it is an integral of exp(-exp(spline(log t))) that no reparameterization removes. core::cox() refuses for the same reason. Read the curve instead, with frm_curve() on the log cumulative hazard.")
  r("royston_parmar", "simulate", "works",
    "The sim slot inverts the fitted survival function by vectorized bisection on log time, so the simulator and the density are one statement of the model rather than two. Verified: five draws over 686 rows are all finite and positive, and refitting one of them as fully observed recovers the coefficient it was drawn from to within 0.35, which is what one draw of that size supports.")
  r("royston_parmar", "residuals", "refused",
    "All three types refuse, and all three for one reason: every residual needs a fitted mean first and this family declares none, so \"pearson\" never reaches the point of missing a variance function. \"deviance\" refuses one step earlier still, on the missing unit deviance. What a Royston-Parmar fit is checked with is its fitted log cumulative hazard against the Nelson-Aalen estimate, which is a plot rather than a residual.")
  r("royston_parmar", "residuals_osa", "untested",
    "One-step-ahead residuals re-tape the objective with the response promoted to a parameter. The spline basis is written branch-free with 0.5*(e + abs(e)) so that it WOULD tape, which is the reason it is written that way, but nothing here exercises the path.")
  r("royston_parmar", "REML", "untested",
    "gamma0 is the primary dpar and would be integrated out. Not exercised.")
  r("royston_parmar", "mixture", "untested",
    "Not exercised. The density is floored rather than allowed to reach NaN, which is the property a mixture needs, so this is the more likely of the untested rows to work.")
  r("royston_parmar", "quadrature", "untested",
    "Not exercised. Marginalizing a random effect by Gauss-Kronrod is orthogonal to what this family does, and nothing about the density argues against it.")

  ## ---- curve inference --------------------------------------------
  r("frm_curve", "s()", "works",
    "The case this package is for. The penalized smooth's wiggly part is a random-effect block, so the curve covariance needs the joint covariance of the fixed and random coefficients; frm_curve() assembles it and checks it against predict(se.fit = TRUE) on every call.")
  r("frm_curve", "smooth", "works",
    "Same thing named as a covariance structure rather than as a formula term.")
  r("frm_curve", "t2()", "works",
    "A t2() tensor smooth is several random-effect blocks rather than one, and the assembly is per COEFFICIENT rather than per block, so nothing changes. Verified on a two-dimensional t2() fit.")
  r("frm_curve", "gp", "works",
    "An exact gp() at an unseen position contributes a kriging variance that is not coefficient uncertainty at all. frm_lp_basis() returns it separately, as extra_var, rather than folding it into A V A', and this package adds it to the diagonal where it belongs. Measured on y ~ gp(x) over 90 points with a grid offset off the observed positions: agreement 0.0e+00, extra_var 7.3e-07 to 8.5e-07.")
  r("frm_curve", "rr", "works",
    "Both re.form settings, since this package reads frmtmb::frm_lp_basis() instead of rebuilding the design. A reduced-rank block's loadings live in theta, so a design over (beta, b) alone is INCOMPLETE; the unit-perturbation rebuild this package used could not see the derivative with respect to the loadings and its standard errors came out 27 percent away from predict(se.fit = TRUE)'s at re.form = NULL, which the covariance check refused. frm_lp_basis() carries the loading columns through rr_jacobians(). Measured on a 5-trait rr(d = 2) fit over 60 groups: agreement 0.0e+00 at re.form = NA (5 columns) and at re.form = NULL (134 columns, 9 of them theta).")
  r("frm_curve", "predict", "works",
    "frm_curve() is built out of predict(): the design is the difference between predictions one coefficient apart, and the answer is checked against predict(se.fit = TRUE).")
  r("frm_curve", "autoscale", "works",
    "Measured, and it is a consequence of reading core's cached covariance rather than recomputing one. An autoscaled fit carries par_units and its covariance comes from autoscale_sdreport(), which reparameterizes; a fresh RTMB::sdreport() would not, and this package used to call one. On y ~ s(x, k = 10) + z with z at scale 1e6, par_units spans 9.7e-07 to 1 and frm_curve() agrees with predict(se.fit = TRUE) to 3.7e-16 absolute, cov_rel_error 6.7e-15.")
  r("frm_curve", "REML", "untested",
    "Under REML the fixed effects are integrated out and beta joins the joint precision as a random block. The assembly reads the joint precision by row name and should follow it there, and the built-in check would catch it if it did not. Not exercised.")
  r("frm_curve", "mvbf", "untested",
    "A multivariate fit reaches the right linear predictor through resp =, which is passed through to predict(). Not exercised.")
  r("frm_curve", "nl", "works",
    "The consumer dev/spline-seam-proposal.md was written for. A nonlinear body is not linear in its coefficients, so the linearity probe this package used to run refused it; frm_lp_basis() tapes the body instead and returns d eta / d coef as a Jacobian. predict(se.fit = TRUE) is still refused for a nonlinear predictor, so there is NO second route to check against: cov_rel_error comes back NA and print() says the check did not run. Measured on the warped-growth model of D'Alessandro, Thoresen and Sorensen (2026) fitted to brokenstick::smocc_200, and on a smaller ps() fit whose Jacobian agrees with a central difference of predict() to 1e-5.")
  r("rp_floored", "cens()", "works",
    "It still reports the censored rows whose fitted -log S passes 19.2, and since frmtmb 0.52.0 it does not REFUSE for them: the family supplies lccdf and the term is exact there. The count is kept because a censored row whose fitted survival probability is exp(-40) is one the data barely constrain, whatever the arithmetic does. What refuses is the monotonicity floor.")
  r("rp_floored", "royston_parmar", "works",
    "The only family it applies to; it refuses any other by name. Both floors are counted: the censored-row one above, and the rows whose fitted d(eta)/d(log t) is non-positive, where the density is a floor rather than a density.")
  r("rp_floored", "frm_curve", "works",
    "frm_curve() and its two companions call rp_floored() on a royston_parmar fit before they assemble anything, so the documented way to inspect this family refuses a fit whose likelihood is a floor artifact. Every other family passes straight through.")
  b$rules()
}
