# What this package tells frmtmb about itself at load time.
#
# The registration uses the seam frmtmb exports for the purpose
# (?frmtmb::`frmtmb-sampling-api`, section "The compatibility registry").
# Registering from .onLoad() rather than at top level is what a
# contributor outside the package must do: by then every namespace is
# sealed, so the collation-order question that governs frmtmb's own
# in-package contributors does not arise.
#
# The rows are the ones core and frmtmb.spline registered for these
# features before the two families moved here, word for word, so that
# frm_compat() reports the same pairs with this package loaded as it
# did before (dev/splitsurv-compat-dump.R compares the two sessions).

#' @noRd
.onLoad <- function(libname, pkgname) {
  # brms's spelling of the cox() basis, bhaz(df = 6) on the response
  # side. It is registered before the rules, so that cox()'s rules see
  # it as a term cox() takes
  frmtmb_register_aterm("bhaz", option = cox_bhaz_option)
  frmtmb_register_compat(
    features = c(cox = "family", royston_parmar = "family",
                 rp_floored = "method"),
    rules = sv_compat_rules)
  invisible()
}

#' The compatibility rules for what this package supplies.
#'
#' Every row was run. `untested` is used where the pair was not
#' exercised, which is the honest third state: an absent guard and a
#' passing guard look the same from outside.
#'
#' @noRd
sv_compat_rules <- function() {
  b <- compat_rule_builder()
  r <- b$r

  ## ---- cox() ------------------------------------------------------
  r("cens()", "cox", "works",
    "The reason the family exists. An event contributes the density log h0(t) + eta - H0(t) exp(eta) and a censored row the survivor function, which is what cox()'s log-density and log-CDF already are, so right, left and interval censoring all run through the ordinary cens() machinery. Verified against a hand-rolled M-spline PH likelihood exactly and against survival::coxph() to 2e-2 on the log hazard ratio.")
  r("trunc()", "cox", "conditional",
    "Runs through the same log-CDF the censoring uses, so trunc(lb = ) is delayed entry. The truncation bound is evaluated against the SAME spline basis the response is, which means a bound outside the boundary knots is clamped to them rather than extrapolated. Untested against an external left-truncated reference.")
  r("bhaz()", "cox", "works",
    "brms's spelling of the baseline basis: bhaz(df =, degree =, intercept =, knots =, Boundary.knots =) sets the cox() argument of the same name, and the two spellings give the same fit. The basis equals brms 2.23.0's standata() Zbhaz and Zcbhaz to rounding, and the log likelihood equals brms's log_lik() at brms's own draws. An argument set in both places with two values is refused, and so are df and knots that disagree, a pair brms resolves by ignoring df. bhaz(gr = ), brms's stratified baseline, is refused by name, and so are boundary knots that do not cover the observed times, where brms extrapolates the basis to a cumulative hazard that can be negative. A shortened argument name brms reads (d = for df, kn =, deg =) is read the same way; one brms drops without a message (int =, Bound =) is refused with the name that was probably meant.")
  r("bhaz()", "royston_parmar", "refused",
    "Refused: bhaz() sets the M-spline baseline of cox(), and a Royston-Parmar model has no such baseline; its spline is set by royston_parmar(df =, knots =, bknots =). frm() stops with \"family 'royston_parmar' has none\".")
  r("mvbf", "cox", "refused",
    "Refused: the Cox baseline is an extra parameter block that the post-fit methods read without a response in hand. Fit the survival response in a model of its own.")
  r("simulate", "cox", "refused",
    "Refused: drawing a survival time means inverting the cumulative baseline hazard, which this family does not carry a quantile function for. simulate(), posterior_predict() and frm_simulate() each say so in their own words and then repeat the family's reason.")
  r("predict", "cox", "conditional",
    "type = \"response\" and fitted() are refused: a survival time has no mean the censored rows identify. type = \"link\" gives the log hazard ratio, and cox_baseline() the fitted baseline weights.")
  r("fitted", "cox", "refused",
    "Refused: a survival time has no mean on the response scale here. Use frm_linpred(type = \"link\") for the log hazard ratio.")

  ## ---- royston_parmar() --------------------------------------
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
  r("royston_parmar", "us", "works",
    "A (1 | centre) block on mu is a shared log-normal frailty: mu is gamma0, so the centre deviation adds to the whole log cumulative hazard and multiplies the cumulative hazard by exp(b). Checked against rstpm2's own log-normal frailty, stpm2(cluster =, RandDist = \"LogN\"), which is the same model in another basis: the change of basis from nsx() coordinates to Royston and Parmar's is exact to 3.0e-14 of the linear predictor's scale. Over 200 replicates at 2000 subjects in 40 centres, 40 percent censored, the treatment coefficient agrees to 3.2e-04 of the run's own standard error on average and 8.8e-04 at worst, the frailty standard deviation to 1.4e-03 relative, and the standard errors on that coefficient to 2.6e-06 relative. The two report DIFFERENT log likelihoods and the difference is the integration rule, not the model: frmtmb's Laplace value is 0.058 units under an exact per-cluster integral where rstpm2's 9-node adaptive Gauss-Hermite is 4.9e-06 under it. That offset costs the ESTIMATE 5.5e-05 units of exact log likelihood, which is 0.1 percent of the offset itself, and it grows as the cluster shrinks, reaching 1.38 units at 4 subjects per centre. A block on gamma1 or higher is a per-centre TIME-VARYING effect and not a frailty; ?royston_parmar says what it means and what it recovers.")
  r("royston_parmar", "|ID|", "works",
    "(1 | c | centre) on mu and gamma1 together fits a correlated centre intercept and centre slope in log time, and it is the spelling that makes a random effect on gamma1 mean the same thing in every time unit. A block on gamma1 ALONE is anchored at t = 1, because that is where log t is zero, and a time unit is not a fact about the data. Measured on one 2000-row dataset read in years, months and days, comparing log likelihoods as logLik + n_event log c, which is what a density owes a rescale: the paired block moves by 4.3e-08 units and holds sd(gamma1 | centre) at 0.2259 in all three, while the slope-only block moves by 8.19 units and its sd(gamma1 | centre) falls from 0.2201 to 0.0370.")
  r("royston_parmar", "predict", "works",
    "type = \"link\" gives the linear predictor of any spline coefficient, which is what frm_curve() reads. type = \"response\" is refused, and correctly: see the fitted row.")
  r("royston_parmar", "fitted", "refused",
    "The family declares no post$mean_fn, because the mean of a Royston-Parmar survival time has no closed form: it is an integral of exp(-exp(spline(log t))) that no reparameterization removes. cox(), in this package too, refuses for the same reason. Read the curve instead, with frm_curve() on the log cumulative hazard.")
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

  r("rp_floored", "cens()", "works",
    "It still reports the censored rows whose fitted -log S passes 19.2, and since frmtmb 0.52.0 it does not REFUSE for them: the family supplies lccdf and the term is exact there. The count is kept because a censored row whose fitted survival probability is exp(-40) is one the data barely constrain, whatever the arithmetic does. What refuses is the monotonicity floor.")
  r("rp_floored", "royston_parmar", "works",
    "The only family it applies to; it refuses any other by name. Three counts: the deep censored rows above; the EVENT rows whose fitted d(eta)/d(log t) is non-positive, where the density is a floor; and, after frmtmb.spline 0.7.0, the CENSORED rows where it is non-positive, where the fitted survival function rises. The last is the only row a group with no events has. It refuses on the event rows and warns, with the size of the rise, on the censored ones.")
  r("rp_floored", "frm_curve", "works",
    "frm_curve() and its two companions call rp_floored() on a royston_parmar fit before they assemble anything, so the documented way to inspect this family refuses a fit whose likelihood is a floor artifact. Every other family passes straight through.")
  hand <- b$rules()
  # cox() declares the addition terms it reads, and each one it does not
  # read is refused by name at frame assembly; core derives those rows
  # from the declaration, as it did while cox() was a core family
  rbind(hand, compat_aterm_rules(list(cox = cox()), existing = hand))
}
