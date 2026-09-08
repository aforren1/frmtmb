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
  frmtmb_register_compat(
    features = c(wiener = "family", gddm = "family", lba = "family",
                 "dec()" = "aterm"),
    rules = ddm_compat_rules)
  # The spelling every reference on the drift-diffusion model uses, and
  # the one brms takes. Before frmtmb had an addition-term registry the
  # indicator had to travel as vint(), which meant the user hand-coded
  # a factor to 0/1 for a family that could have done it. This is that
  # coercion, contributed once.
  frmtmb_register_aterm("dec", arity = 1L, coerce = ddm_coerce_dec)
  # The two families added after the first three, registered in their own
  # call rather than folded into the one above. Two reasons: the rows
  # below were measured on this worktree and the ones above were not, and
  # a separate call keeps the two sets of rows from having to be merged
  # by hand when the lanes meet.
  frmtmb_register_compat(
    features = c(rdm = "family", wiener_gng = "family"),
    rules = rdm_gng_compat_rules)
  invisible()
}

#' Coerce a decision indicator to the 0/1 the density reads.
#'
#' brms's rule, kept: a factor or character vector is read on its
#' levels, and the SECOND level is the upper boundary, so
#' `c("lower", "upper")` and `c(FALSE, TRUE)` both read the way they
#' look. A numeric column is passed through and validated by the family,
#' which is where a 2 rather than a 1 gets its refusal.
#'
#' @noRd
ddm_coerce_dec <- function(x) {
  if (is.factor(x) || is.character(x) || is.logical(x)) {
    f <- if (is.factor(x)) x else factor(x)
    if (nlevels(f) != 2L) {
      stop("dec(): a decision indicator has two levels, one per ",
           "boundary, and this one has ", nlevels(f), ": ",
           paste(levels(f), collapse = ", "),
           ". Give it as a factor carrying both levels, or as a 0/1 ",
           "column.", call. = FALSE)
    }
    return(as.numeric(as.integer(f) - 1L))
  }
  as.numeric(x)
}

#' The compatibility rules for the families this package supplies.
#'
#' Every `wiener` row here was run, not reasoned about. `untested` is
#' used where the pair was not exercised, which is the honest third
#' state: an absent guard and a passing guard look the same from
#' outside.
#'
#' The `gddm` and `lba` rows are the exception and are marked as such
#' where they sit: those two families were written in sibling worktrees
#' and cannot be loaded from this one, so their rows are transcribed
#' from the review that read all three lanes together. They become
#' testable at the merge, and the merged suite is what confirms them.
#'
#' @noRd
ddm_compat_rules <- function() {
  b <- compat_rule_builder()
  r <- b$r
  r("wiener", "dec()", "works",
    "The spelling to use, contributed to frmtmb's addition-term registry when this package loads. It takes a factor, a character vector or a logical the way brms does, reading the second level as the upper boundary, and it is required: the boundary a trial ended at is data and the density is meaningless without it.")
  r("wiener", "vint()", "works",
    "The other spelling for the same thing, carrying the indicator as a plain 0/1 integer. It was the only route before frmtmb had an addition-term registry and it still works unchanged. Supplying neither is refused, because the density would otherwise read a NULL and the log likelihood would silently collapse to zero terms.")
  r("wiener", "cens()", "refused",
    "The family declares no lcdf. The Wiener first-passage distribution function is a third series with its own truncation problem and none of it is written here.")
  r("wiener", "trunc()", "refused",
    "Same reason as cens(): no lcdf, so frmtmb has no normalizing constant to divide by.")
  r("wiener", "weights()", "works",
    "Verified: the weighted log likelihood is the unweighted one at unit weights and scales as it should.")
  r("wiener", "simulate", "works",
    "Verified against the fitted parameters. Draws are conditional on each row's boundary, by inverse transform through RWiener's defective quantile function; without RWiener a discretized forward simulation stands in. Under variability = the per-trial parameters are drawn first and accepted with the boundary probability they imply, because conditioning on the boundary reweights which of them the trial could have had.")
  r("wiener", "fitted", "works",
    "The mean is the conditional mean response time for the row's own boundary, in closed form. Under variability = it is a ratio of two quadratures instead, for the same reason simulate() rejects: the closed form is the mean of a model whose parameters do not vary, and at an unbiased start point it cannot even tell the two boundaries apart.")
  r("wiener", "predict", "works",
    "The decision indicator is mandatory on newdata, under whichever of dec() and vint() the model was written with, so the boundary must be supplied there as well.")
  r("wiener", "residuals", "conditional",
    "type = \"response\" works. \"pearson\" and \"deviance\" are both refused: the family declares no variance function and no unit deviance, because the conditional variance of a first-passage time was not worth writing for a residual nobody reads on response times.")
  r("wiener", "residuals_osa", "untested",
    "One-step-ahead residuals re-tape the objective with the response promoted to a parameter. Nothing here exercises that path.")
  r("wiener", "REML", "untested",
    "The drift rate is the primary dpar and would be integrated out. Not exercised.")
  r("wiener", "mixture", "works",
    "Verified with a lognormal contaminant. The Wiener component needs wiener(max_ndt =, allow_unreachable = TRUE) so that its non-decision time may sit above the fastest response time, which is the whole point of putting a contaminant there.")
  r("wiener", "quadrature", "refused",
    "By frmtmb, and correctly. quadrature = TRUE marginalizes RANDOM EFFECTS by Gauss-Kronrod and refuses a model with no random-effect block. The across-trial variability of wiener(variability =) is a different integral entirely, over per-row parameter distributions that no level is shared across, and it is done inside the density.")

  # The other two families in this package. Their rows are written from
  # the review's enumeration rather than from this file's own runs,
  # because gddm() and lba() were built in sibling worktrees and are not
  # loadable here; the merge is where they are exercised. Every row
  # below is a statement about the merged package, and any of them that
  # the merged suite contradicts is a defect in this table, not in the
  # family.
  r("gddm", "dec()", "works",
    "The boundary is read from dec() when it is there and from vint() otherwise, as wiener() does. vint() numbers positionally: alongside dec() the condition index is the first vint() value, and inside vint(upper, cond) it is the second.")
  r("gddm", "vint()", "works",
    "Required, and twice over: vint1 is the boundary a trial ended at, coded 0/1, and vint2 is the condition index the solver groups on. Both are declared, so omitting either is refused by name rather than silently summed over no rows.")
  r("gddm", "vreal()", "works",
    "Carries the per-condition covariate the drift nonlinearity reads.")
  r("gddm", "cens()", "refused",
    "No lcdf. The generalized model's distribution function would be a second pass over the Fokker-Planck solve and is not written.")
  r("gddm", "trunc()", "refused",
    "No lcdf, so there is no normalizing constant for frmtmb to divide the window by.")
  r("gddm", "weights()", "untested",
    "Not exercised. Nothing in the family works against case weights, but nothing has run them either.")
  r("gddm", "simulate", "works",
    "gddm_simulate() draws from the same solved density the likelihood scores, so the simulator and the density are one statement of the model.")
  r("gddm", "fitted", "works",
    "The family defines post$mean_fn, so the mean response time comes back on the response scale.")
  r("gddm", "predict", "works",
    "Both vint() columns are mandatory on newdata, so the boundary and the condition must be supplied there as well.")
  r("gddm", "residuals", "conditional",
    "type = \"response\" works. \"pearson\" and \"deviance\" are refused for the same reason as wiener: no variance function and no unit deviance.")
  r("gddm", "residuals_osa", "untested",
    "One-step-ahead residuals re-tape the objective with the response promoted to a parameter. Nothing exercises that path.")
  r("gddm", "REML", "untested",
    "Not exercised.")
  r("gddm", "mixture", "untested",
    "Not exercised, and worth exercising before it is relied on: gddm does not floor its density where wiener does, so a below-support row is a NaN rather than a finite zero, and a NaN inside a log-sum-exp takes every other component with it.")
  r("gddm", "quadrature", "refused",
    "By frmtmb, for the same reason it refuses wiener: quadrature = TRUE integrates random effects, and this family has no random effect to integrate.")

  r("lba", "dec()", "refused",
    "Refused by the family, and by declaration rather than by a check of its own: lba() names the addition terms it takes in frmtmb_family(accepts_aterms =), dec() is not among them, and frame assembly refuses it by name. Until frmtmb 0.53.0 there was no such declaration and the refusal was a hand-written check shared with rdm(); before that check this row was wrong in a way worth recording, because a model written rt | dec(two) + vint(choice) FITTED, dropping the term with no warning and with fixed effects bit-identical to the model without it. dec() IS the spelling under wiener(), so a ported model quietly ignored half of what its author wrote.")
  r("lba", "vint()", "works",
    "Required: vint1 is which accumulator won, counted from 1. Note that this is 1-based where wiener's boundary indicator is 0-based, which is a difference between the two families and not a typo.")
  r("lba", "cens()", "refused",
    "No lcdf. The race distribution function is a product of survivals with no closed form written here.")
  r("lba", "trunc()", "refused",
    "Same reason as cens(): no lcdf, so there is no normalizer.")
  r("lba", "weights()", "untested",
    "Not exercised.")
  r("lba", "simulate", "works",
    "The family supplies a sim slot that races the accumulators and returns the winner's time.")
  r("lba", "fitted", "refused",
    "The family declares no post$mean_fn, because the mean of the winning accumulator's time has no closed form: it is an expectation over the minimum of n truncated-normal-rate arrivals.")
  r("lba", "predict", "works",
    "vint() is mandatory on newdata, so the winning accumulator must be supplied there as well.")
  r("lba", "residuals", "conditional",
    "type = \"response\" works; the two standardized types are refused for want of a variance function.")
  r("lba", "residuals_osa", "untested",
    "Not exercised.")
  r("lba", "REML", "untested",
    "Not exercised.")
  r("lba", "mixture", "untested",
    "Not exercised.")
  r("lba", "quadrature", "refused",
    "By frmtmb. No random effect, nothing to marginalize.")
  # The refusals each family has already declared, appended and handed
  # the rows above so that it defers to them: a pair written by hand
  # keeps its own note, and what is added is the cells nobody wrote.
  # Those are the terms outside a family's accepts_aterms allow-list,
  # which frame assembly refuses by name, so they stop reading
  # `untested` when the guard has been there all along.
  hand <- b$rules()
  rbind(hand,
        compat_aterm_rules(ddm_accepts[c("wiener", "gddm", "lba")], hand))
}

#' The compatibility rules for the two families added at 0.3.0.
#'
#' Unlike the `gddm` and `lba` rows above, every row here was RUN on the
#' worktree that wrote it, by a probe that calls the feature and records
#' what came back. `untested` below therefore means the pair was
#' deliberately not exercised, and says why, rather than meaning nobody
#' looked.
#'
#' @noRd
rdm_gng_compat_rules <- function() {
  b <- compat_rule_builder()
  r <- b$r

  r("rdm", "vint()", "works",
    "Required: vint1 is which accumulator reached the threshold, counted from 1, exactly as it is for lba(). Omitting it is refused by name, because the density indexes it.")
  r("rdm", "dec()", "refused",
    "Refused by declaration: rdm() names the terms it takes in frmtmb_family(accepts_aterms =) and dec() is not one of them, because dec() carries a two-level boundary indicator and a race of n accumulators needs a winner in 1..n. The winner travels through vint(). lba() refuses it through the same seam, so the two race families cannot drift apart on it.")
  r("rdm", "vreal()", "refused",
    "Refused by name since frmtmb 0.53.0. Nothing in the density reads a real-valued addition term, and this row used to read works, on the ground that a model supplying one fitted and gave the same answer as one that did not. That IS the defect: the column travelled into the fit and changed nothing, silently. rdm() now declares the terms it takes and frame assembly refuses the rest.")
  r("rdm", "cens()", "works",
    "All four codes. The family declares both an lccdf and an lcdf, and neither needed new algebra: the race is unfinished exactly when every accumulator is, so log S is the sum of the same per-accumulator survivals the density already forms for the losers of an observed trial, over all n instead of n - 1. Verified against the likelihood written by hand at the fitted parameters: right censoring alone on 400 rows of which 100 are censored agrees to 2.2e-15 relative, and all four codes together on 300 rows to 5.0e-16. A censored row still needs a vint() winner, which the likelihood does not read; that is the price of a declaration that cannot be conditional on a censoring code.")
  r("rdm", "trunc()", "conditional",
    "Works, and the range matters. The family declares an lcdf, written as -expm1(log S) so that the distribution function of a race whose survival is within a rounding of one does not come back as exactly zero; a fit left-truncated at 0.30 on 1146 surviving rows reproduces the hand-written normalized likelihood to 9.4e-16, and one right-truncated at 1.20 to 4.1e-16. What limits it is not this family: core forms the left-truncation normalizer as 1 - F(lb) on the PROBABILITY scale (R/objective.R, then ll - log(Fub - Flb)), so the declared lccdf is not used for it and the cancellation the lccdf exists to remove comes back at the bound. Measured on rdm(2), v = (3, 2), A = 0.8, k = 0.5, ndt = 0.15: exact at log S(lb) = -0.406, 5.0e-11 at -14.0, 1.1e-07 at -21.5, 8.7e-05 at -28.7 and 20 percent wrong at -35.8. A left bound in the ordinary range is exact; a bound out in the tail, past about log S = -25, is not. Core documents the limitation at R/families.R:278-285 and calls closing it a windowed log-difference slot it does not yet have.")
  r("rdm", "weights()", "works",
    "Verified: the weighted log likelihood equals the unweighted one at unit weights.")
  r("rdm", "simulate", "works",
    "The family supplies a sim slot that races the accumulators from inverse-Gaussian draws and returns the winner's time, conditioned on each row's observed vint() winner by rejection. rdm_simulate() is the unconditional joint draw of choice and time.")
  r("rdm", "fitted", "refused",
    "The family declares a post$mean_fn that stops, because the mean of the winning accumulator's arrival is an expectation over the minimum of several inverse-Gaussian first passages and has no closed form. The refusal is deliberate and replaced a silent wrong answer: with that slot empty frmtmb returns the first primary dpar on the response scale, which here is a DRIFT RATE, and predict(type = \"response\") gave 3.51 for data whose response times average 0.36.")
  r("rdm", "predict", "conditional",
    "type = \"link\" works, on the training data and on newdata. type = \"response\" is refused with the same message fitted() gives. Note what is NOT true: vint() is not mandatory on newdata for a link-scale prediction, because no linear predictor reads it.")
  r("rdm", "residuals", "refused",
    "All three types. \"response\" needs the mean the family refuses; \"pearson\" reports that same refusal, because it asks for the mean before it asks for a variance function; \"deviance\" is refused by frmtmb for want of a unit deviance. Before the refusing mean was added, \"response\" returned a length-zero vector rather than an error.")
  r("rdm", "residuals_osa", "refused",
    "Reached and refused, but not gracefully: it fails inside RTMB with a type error about S4 and double rather than with a sentence naming the family. One-step-ahead residuals re-tape the objective with the response promoted to a parameter, which this density does not survive.")
  r("rdm", "REML", "works",
    "Verified: REML = TRUE fits. The drifts are the primary dpars and are integrated out.")
  r("rdm", "quadrature", "works",
    "Verified on a model with a random effect. This differs from the wiener row above, which records a refusal: frmtmb refuses quadrature only when there is no random-effect block to marginalize, and a racing-diffusion model with a grouping factor has one.")
  r("rdm", "mixture", "refused",
    "By frmtmb: mixture() components need a dpar called mu, and this family's primary dpars are v1..vn. A contaminant on a race would need drifts named the way mixture() expects, which is a change to mixture() rather than to this family.")

  r("wiener_gng", "dec()", "works",
    "Required, and the only spelling: dec() says whether a trial produced a response, coded 1, or did not, coded 0. With one observable boundary that is the same 0/1 wiener() reads and it means the same thing. A logical column gives the same fit as a 0/1 one, verified.")
  r("wiener_gng", "vint()", "refused",
    "Refused by name, so that a model written against wiener()'s vint(upper) spelling fails loudly rather than fitting with the indicator ignored. The response indicator travels through dec() here.")
  r("wiener_gng", "vreal()", "works",
    "Carries the per-row deadline when it varies between trials, and is then required. A constant deadline goes on the family as wiener_gng(deadline =) instead; supplying neither is refused by name, and the two spellings give the same log likelihood to 1e-10, verified.")
  r("wiener_gng", "cens()", "conditional",
    "RIGHT censoring works, and is the same statement the family already makes: a trial whose clock stopped before it responded is a no-go trial with the deadline moved, so the declared lccdf is this family's own no-go probability at the censoring time. Verified to the LAST BIT: the same 400 rows scored as no-go trials at the deadline and as trials right-censored at the deadline give log likelihoods that differ by exactly zero. Left and interval censoring are refused by frmtmb, because the family declares no lcdf, and that is deliberate rather than unwritten; see the trunc() row.")
  r("wiener_gng", "trunc()", "refused",
    "By frmtmb, because the family declares no lcdf, and the family declares none on purpose. This likelihood is a defective density plus a point mass at no-response, and a truncation window on the response scale renormalizes the density while saying nothing about the mass, so the two halves of every row would be divided by different things. Right censoring has no such problem, because it replaces a row rather than reweighting it.")
  r("wiener_gng", "weights()", "works",
    "Verified: the weighted log likelihood equals the unweighted one at unit weights.")
  r("wiener_gng", "simulate", "works",
    "The sim slot redraws a go trial's time by rejection from the unconditional process and gives a no-go trial its deadline, which is the placeholder the family documents. wiener_gng_simulate() is the unconditional draw of outcome and time together.")
  r("wiener_gng", "fitted", "refused",
    "The family declares a post$mean_fn that stops, and here the refusal is a statement about the model rather than a missing integral: a go/no-go trial produces a PAIR, whether a response happened and when, and the no-go rows have no response time to average at all. It also replaced a silent wrong answer, and a worse one than rdm's because this family's primary dpar is literally called mu: fitted() returned a constant drift of 1.05 for data whose go response times average 0.6.")
  r("wiener_gng", "predict", "conditional",
    "type = \"link\" works, on the training data and on newdata. type = \"response\" is refused with the same message fitted() gives. dec() is not mandatory on newdata for a link-scale prediction, because no linear predictor reads it.")
  r("wiener_gng", "residuals", "refused",
    "All three types, and for the reasons fitted() gives: \"response\" and \"pearson\" both report the mean's refusal, and \"deviance\" is refused by frmtmb for want of a unit deviance.")
  r("wiener_gng", "residuals_osa", "refused",
    "Reached and refused inside RTMB rather than by a sentence naming the family, exactly as for rdm.")
  r("wiener_gng", "REML", "works",
    "Verified: REML = TRUE fits.")
  r("wiener_gng", "quadrature", "works",
    "Verified on a model with a random effect, as for rdm and unlike wiener().")
  r("wiener_gng", "mixture", "conditional",
    "It assembles and runs, where mixture(rdm(3), ...) is refused outright, because this family has a dpar called mu. But the one case tried did not converge, reporting false convergence and a maximum absolute gradient of 7.7e13, so nothing here supports relying on it. What a contaminant should do with the no-go rows is a modelling question this package has not answered.")
  hand <- b$rules()
  rbind(hand,
        compat_aterm_rules(ddm_accepts[c("rdm", "wiener_gng")], hand))
}
