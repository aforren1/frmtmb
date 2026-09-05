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
  r("royston_parmar", "cens()", "conditional",
    "Works only while every censored row stays inside the accurate region, and the package cannot enforce that during the fit. frmtmb forms a right-censored term as log(1 - F(y)) on the probability scale (R/objective.R:100) and offers a family no complementary log-CDF slot, so the scored log S carries absolute error about eps/S: it is exact to 1.3e-13 at -log S of 10, wrong by 1.7e-04 at 30, floored at -35.127363 past 36, and its GRADIENT is exactly zero past 30. The size of the resulting error is not a property of the family but of the data: a censored row scored at the floor contributes -35.127363 instead of its own -log S, so the reported log likelihood is short by about (-log S - 35) for every such row, which is thousands to tens of thousands as soon as one censored time sits well past the event times. Two independent runs of the same 600-subject design differ only in seed and give 2.4e+03 and 2.166e+04. In both the fit converged without a warning and the treatment coefficient was out by tens of percent, on data flexsurv declines to fit at all. rp_floored() refuses such a fit post-hoc and frm_curve() calls it; logLik() and AIC() cannot be gated, because they read the optimizer value and no family hook runs at fit end. The identity test against flexsurv reaches only -log S = 2.01, so it cannot see this. Unconditional again once core has an lccdf slot: see dev/spline-seam-proposal.md.")
  r("royston_parmar", "trunc()", "conditional",
    "Same lcdf and the same limit, from both ends. A truncation window is divided out as log(Fub - Flb) (R/objective.R:110), so an upper bound deep in the tail loses accuracy exactly as a censored row does, and a LEFT truncation bound (delayed entry) runs into the same representability problem from the other side, which an lccdf slot alone would not close. Verified for the case that must change nothing: a bound below every observed time reproduces the untruncated coefficients to 0.05.")
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
  r("frm_curve", "gp", "untested",
    "An exact gp() term adds a kriging variance to predict(se.fit = TRUE) that is not a coefficient covariance at all, so the check inside frm_curve() would fail at an unseen position and refuse. Not exercised, and the refusal is the right behavior if it does.")
  r("frm_curve", "rr", "conditional",
    "re.form = NA works and was measured at 1.2e-14: the reduced-rank block drops out of a population prediction, so it contributes nothing to rebuild. re.form = NULL is REFUSED, and by the covariance check rather than by the linearity probe. A reduced-rank block reaches the linear predictor through loadings that live in theta, so eta is linear in b at fixed theta and the probe passes; what the perturbation cannot see is the derivative with respect to the loadings, which core's own delta method carries as rr_jacobians(). Measured: the assembled standard errors are 27 percent away from predict(se.fit = TRUE)'s, and the call refuses rather than reporting them. This row is why the check exists.")
  r("frm_curve", "predict", "works",
    "frm_curve() is built out of predict(): the design is the difference between predictions one coefficient apart, and the answer is checked against predict(se.fit = TRUE).")
  r("frm_curve", "autoscale", "works",
    "Measured, and it is a consequence of reading core's cached covariance rather than recomputing one. An autoscaled fit carries par_units and its covariance comes from autoscale_sdreport(), which reparameterizes; a fresh RTMB::sdreport() would not, and this package used to call one. On y ~ s(x, k = 10) + z with z at scale 1e6, par_units spans 9.7e-07 to 1 and frm_curve() agrees with predict(se.fit = TRUE) to 3.7e-16 absolute, cov_rel_error 6.7e-15.")
  r("frm_curve", "REML", "untested",
    "Under REML the fixed effects are integrated out and beta joins the joint precision as a random block. The assembly reads the joint precision by row name and should follow it there, and the built-in check would catch it if it did not. Not exercised.")
  r("frm_curve", "mvbf", "untested",
    "A multivariate fit reaches the right linear predictor through resp =, which is passed through to predict(). Not exercised.")
  r("frm_curve", "nl", "untested",
    "A nonlinear body is not linear in its coefficients, so the linearity probe should refuse it. That is the right answer for the construction as it stands and the wrong answer for the model: a warped spline inside a nonlinear body is exactly the consumer dev/spline-seam-proposal.md is written for.")
  r("rp_floored", "cens()", "works",
    "The reason it exists. It recomputes -log S on every censored row at the fitted parameters and refuses past 19.2, naming the row count, the maximum, the threshold, the reason and the remedy. Pinned by a 600-subject fit with one subject censored far beyond the event times, which converges without a warning and which this function refuses.")
  r("rp_floored", "royston_parmar", "works",
    "The only family it applies to; it refuses any other by name. Both floors are counted: the censored-row one above, and the rows whose fitted d(eta)/d(log t) is non-positive, where the density is a floor rather than a density.")
  r("rp_floored", "frm_curve", "works",
    "frm_curve() and its two companions call rp_floored() on a royston_parmar fit before they assemble anything, so the documented way to inspect this family refuses a fit whose likelihood is a floor artifact. Every other family passes straight through.")
  b$rules()
}
