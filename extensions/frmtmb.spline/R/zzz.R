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
                 frm_curve_contrast = "method", rp_floored = "method"),
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
  ## ---- the difference curve ---------------------------------------
  r("frm_curve_contrast", "s()", "works",
    "The case it is for. On y ~ fac + s(x, by = fac, k = 8) at n = 400 the difference agrees with gratia::difference_smooths(group_means = TRUE) on the same mgcv ML fit to 1.95e-08, which is 1.5 times the gap between the two packages' own fitted curves (1.29e-08) and so is the fit difference rather than a method difference. The attribution is provable rather than argued: the difference's disagreement with gratia equals the difference of the two curves' own disagreements to 5.55e-16. The standard errors agree to 1.56 percent, which sits BETWEEN curve A's 1.49 and curve B's 3.11 rather than below both, and does the same under unconditional = TRUE (1.32 between 1.15 and 1.58). mgcv's Vp is conditional on the smoothing parameter and the joint precision this package inverts is not; nothing measurable cancels, and the difference tracks the better curve. Simultaneous-band coverage of the true difference over the whole grid is 0.970 over 200 seeds (binomial mcse 0.0154) against a pointwise band's 0.685.")
  r("frm_curve_contrast", "smooth", "works",
    "Same thing named as a covariance structure. The two designs must sit at the same rows of the joint covariance for the subtraction to mean anything; coef_pos is compared on every call and a mismatch refuses rather than pairing the wrong columns.")
  r("frm_curve_contrast", "nl", "works",
    "Measured on a two-ps() nonlinear body, y ~ lev + ps(t) + b1 w + b2 w t + ps(z), differencing across w. frm_lp_basis() tapes the body and returns a Jacobian per grid, and the difference of two Jacobians is the Jacobian of the difference because the subtraction is linear. predict(se.fit = TRUE) is refused for a nonlinear predictor, so cov_rel_error is NA here for the same reason it is NA for one curve, and print() says so.")
  r("frm_curve_contrast", "ps()", "works",
    "Both grids are checked against the frozen knot span, and each message says which grid raised it. frm_curve_feature()'s second span gate asks whether a column other than the search variable varies DOWN each grid; it does NOT ask whether the two grids differ from one another, which they always do. Measured on the fixture above: a clean difference curve is accepted with no warning, and an out-of-span row in row 10 of either grid refuses by name.")
  r("frm_curve_contrast", "predict", "works",
    "Two calls rather than one, one per grid, and both are the covariance check. What predict(se.fit = TRUE) cannot supply is the covariance BETWEEN the grids, which is the whole content of a difference, so the check licenses the two designs and not the difference's own standard error.")
  r("frm_curve_contrast", "gp", "conditional",
    "Works when the two grids sit at the SAME gp() positions, which is the ordinary case: both then load the same kriging residual of the same field and it cancels exactly, so the difference is (A1 - A2) V (A1 - A2)' with nothing left over. Measured on y ~ fac + gp(x) over 90 points with the grid offset off the observed positions: the difference across fac is flat at 0.0618114535, which is sqrt(vcov(fit)[\"facB\", \"facB\"]) to 12 digits. Sameness is tested on the DESIGN, not on the variances: the two grids must agree bit for bit on EVERY column of A outside the fixed effects, where equality of extra_var alone would pass two different levels of one grouping block, which have identical marginal variances by construction and are different draws. That test is stricter than the mathematics needs and the refusal says so: only the block carrying the residual has to match, so a contrast across fac on y ~ fac + s(x, by = fac, k = 6) + gp(x) is refused even though its gp() columns are bit-identical, because the by-factor smooth's own columns differ. It fails closed, so that is a missing answer and not a wrong one. Refused too when the two grids sit at different gp() positions, which moving the second grid's x by 1e-10 already triggers. That refusal is a SEAM limit and not a mathematical one: the conditional cross-covariance is k(x1, x2) - Xr1 K Xr2' and core forms every piece of it at R/predict.R:383-405, but reduces the result to one variance per row before frm_lp_basis() returns. Filed as a core seam in dev/diffcurve-findings.md.")
  r("frm_curve_contrast", "frm_lp_basis", "works",
    "The whole implementation is two reads of that seam and one subtraction. Where core has a second route to the answer, the two agree bit for bit: on y ~ fac * x with no random effect the difference standard error is identical() to the one formed from vcov() and the same contrast matrix.")
  r("frm_curve_contrast", "t2()", "untested",
    "Nothing about a tensor smooth argues against it, since the assembly is per coefficient rather than per block, but no difference was taken across one.")
  r("frm_curve_contrast", "rr", "untested",
    "A reduced-rank block's loading columns reach A through rr_jacobians() and would subtract like any other column. Not exercised.")
  r("frm_curve_contrast", "mvbf", "untested",
    "resp = is passed to both grids alike. Not exercised.")

  r("rp_floored", "cens()", "works",
    "It still reports the censored rows whose fitted -log S passes 19.2, and since frmtmb 0.52.0 it does not REFUSE for them: the family supplies lccdf and the term is exact there. The count is kept because a censored row whose fitted survival probability is exp(-40) is one the data barely constrain, whatever the arithmetic does. What refuses is the monotonicity floor.")
  r("rp_floored", "royston_parmar", "works",
    "The only family it applies to; it refuses any other by name. Both floors are counted: the censored-row one above, and the rows whose fitted d(eta)/d(log t) is non-positive, where the density is a floor rather than a density.")
  r("rp_floored", "frm_curve", "works",
    "frm_curve() and its two companions call rp_floored() on a royston_parmar fit before they assemble anything, so the documented way to inspect this family refuses a fit whose likelihood is a floor artifact. Every other family passes straight through.")
  b$rules()
}
