# What this package tells frmtmb about itself at load time: the
# addition terms its families read, and the compatibility rules for
# what those families can and cannot do.
#
# Registering from .onLoad() rather than at top level is what a
# contributor outside the core must do: by then every namespace is
# sealed, so the collation-order question that governs frmtmb's own
# in-package contributors does not arise. The aterms go first, because
# registering a term puts it in the compatibility vocabulary and the
# rules below name it.

#' The addition terms this package brings, and their arities.
#'
#' Three rather than one, because an addition term's arity is fixed when
#' it is registered and a K-armed task needs a K-column term. `reward`
#' is spelled exactly as `frmtmb`'s own reinforcement-learning example
#' registers it, so a session that has sourced that example and then
#' loads this package finds one term rather than two spellings of one.
#'
#' @noRd
ln_aterms <- c(reward = 2L, payoff = 4L, stage2 = 2L)

#' Register each term once, and never over the top of another package.
#'
#' `frmtmb_register_aterm()` REPLACES an existing entry of the same
#' name, which is what makes reloading a package harmless and what makes
#' a foreign registration at a different arity silently destructive. So
#' this checks first. What it cannot check is the arity of a term that
#' is already there: `frm_compat_features()` carries a term's name, key
#' and kind, and no exported accessor returns its arity (recorded as a
#' seam in dev/learn-findings.md). That case is caught one step later
#' instead, and loudly: a family declares `required_aterms` by their
#' indexed names, so a `reward` registered at arity 1 makes `frm()`
#' refuse the fit for a missing `reward2` rather than fit something
#' wrong.
#'
#' @noRd
ln_register_aterms <- function() {
  known <- frmtmb::frm_compat_features()
  known <- known[["key"]][known[["kind"]] == "aterm"]
  for (nm in names(ln_aterms)) {
    if (nm %in% known) next
    frmtmb::frmtmb_register_aterm(nm, arity = ln_aterms[[nm]])
  }
  invisible(NULL)
}

#' @noRd
.onLoad <- function(libname, pkgname) {
  ln_register_aterms()
  frmtmb_register_compat(
    features = c(bandit2arm_delta = "family", bandit2arm_dual = "family",
                 prl_fictitious = "family",
                 bandit4arm2_kalman_filter = "family",
                 ts_par7 = "family", igt_pvl_delta = "family",
                 igt_orl = "family", rlddm = "family",
                 frm_value_trace = "method", frm_task_simulate = "method"),
    # `dec()` is frmtmb.eam's term, not this package's. rlddm() reads it
    # for the boundary the response reached, so a rule here names a
    # feature another package supplies, which is what expects = is for.
    #
    # It IS already registered by the time this call runs, and the
    # reason is the `importFrom(frmtmb.eam, ...)` in NAMESPACE rather
    # than the DESCRIPTION Imports line: R loads a package's imports
    # before running its .onLoad(), and frmtmb.eam registers `dec()` in
    # its own. An earlier version of this comment credited the
    # DESCRIPTION entry, which loads nothing, and `rlddm()` was
    # unusable from a clean session for exactly that reason. `expects =`
    # tolerates a name already present either way, and is what stops
    # this call failing if the ordering ever changes.
    # `ndt_group()` is the same story: frmtmb.eam's term, registered in
    # the same .onLoad(), and read by rlddm() through that package's
    # bound seam.
    expects = c("dec()", "ndt_group()"),
    rules = ln_compat_rules)
  invisible(NULL)
}

#' The rows every family in this package shares.
#'
#' They are shared because the FACT is shared: each one is a consequence
#' of the likelihood not factorizing over rows, or of the recursion
#' reading a whole trial history, and neither depends on which learning
#' rule the family carries. Family-specific rows are added on top of
#' these in `ln_compat_rules()`.
#'
#' @noRd
ln_common_rules <- function(r, nm, nominal = TRUE) {
  r(nm, "weights()", "refused",
    paste0("Refused in check_spec, by name. A row weight has no meaning ",
           "here: a trial's factor depends on every earlier trial of the ",
           "same subject, so there is nothing separable to reweight. The ",
           "structure's loglik slot does take `weights` and ignores them, ",
           "which frmtmb_structure() permits only because the term is ",
           "refused before it can arrive."))
  r(nm, "cens()", "refused",
    paste0("Refused in check_spec, with se() and trunc(), for one ",
           "reason: all three reshape a per-row likelihood contribution ",
           "and this family's contributions are not per row. A censored ",
           "CHOICE is also not a thing the literature defines."))
  r(nm, "trunc()", "refused",
    paste0("Refused in check_spec. Same reason as cens(): a truncation ",
           "window divides a row's density by a probability, and a ",
           "trial's density here is conditional on its subject's whole ",
           "history."))
  r(nm, "se()", "refused",
    "Refused in check_spec. A choice has no measurement standard error.")
  r(nm, "importance", "works",
    paste0("ADMITTED, and it was refused by name until this family ",
           "declared how its likelihood factorizes. The correction ",
           "reweights draws from the Laplace Gaussian and needs one ",
           "log-likelihood value per GROUP; the family now supplies ",
           "frmtmb_structure(loglik_group = ), one value per subject, ",
           "and frmtmb_structure(loglik_row = ), one per trial, both ",
           "off the same recursion the objective tapes. The core checks ",
           "the family's units against the model's grouping factor and ",
           "refuses a model grouped on anything else, and imp_verify() ",
           "checks the family's handling of the stacked design against ",
           "the plain objective at the first freeze. MEASURED, and the ",
           "answer depends on the design rather than on the family: on ",
           "40 subjects by 100 trials the fixed effects move by less ",
           "than 0.02 and the subject-level standard deviation moves up ",
           "substantially where the Laplace fit had not collapsed it; ",
           "on 20-trial data the correction still runs but the quantity ",
           "it corrects is already gone. dev/learn2-findings.md carries ",
           "the per-dataset table with its Monte Carlo error and its ",
           "effective sample sizes, and vignette('learning') carries ",
           "the summary."))
  r(nm, "quadrature", "refused",
    paste0("Refused by the protocol's conservative default, and a real ",
           "refusal rather than a missing slot. Gauss-Kronrod ",
           "integration of one scalar random effect works against a ",
           "PRODUCT of per-row densities, and this likelihood is not a ",
           "product of per-row densities."))
  r(nm, "REML", "refused",
    paste0("Refused by the protocol's conservative default. Not ",
           "exercised, and not obviously meaningful: REML integrates the ",
           "fixed effects out of a linear predictor, and every ",
           "predictor here feeds a nonlinear recursion."))
  r(nm, "profile", "refused",
    "Refused by the protocol's conservative default. Not exercised.")
  r(nm, "mvbf", "refused",
    paste0("Refused in check_spec, in the family's own words. The ",
           "recursion is a likelihood over one response's trial ",
           "sequences."))
  r(nm, "rescor", "refused", "Refused with mvbf(), by the same check.")
  r(nm, "mixture", "refused",
    paste0("Refused by the protocol's conservative default. A mixture ",
           "over learning strategies is a real model and a wanted one, ",
           "but core's mixture() combines per-row densities and would ",
           "need per-SEQUENCE ones here. The family now declares its ",
           "per-sequence values, which is what admitted the importance ",
           "correction, so what is left is a mixture() that reads a ",
           "declared factorization rather than a rowwise density. That ",
           "is a change under core's R/ and not a declaration this ",
           "package can make."))
  r(nm, "residuals_osa", "refused",
    paste0("Refused by name. One-step-ahead residuals re-tape the ",
           "objective with the response promoted to a parameter, and the ",
           "response here is a category the recursion selects with a ",
           "zero-one indicator rather than a continuous quantity the ",
           "tape could differentiate."))
  if (nominal) {
    r(nm, "fitted", "refused",
      paste0("The family declares no fitted_mean, and that is a ",
             "decision rather than an omission. The response is the ",
             "option a subject took, coded 1 to K, and it is NOMINAL: ",
             "arm 2 is not twice arm 1, and the four decks of the Iowa ",
             "gambling task have no order at all. Core forms a fitted ",
             "value as the conditional mean of the response and a ",
             "residual as y - mean, so any mean this family supplied ",
             "would be arithmetic on a category code. core::cox() and ",
             "frmtmb.spline::royston_parmar() decline to invent a mean ",
             "for the same kind of reason. WHAT REPLACES IT: ",
             "frm_value_trace(), which returns more than a mean could. ",
             "Its `p` column is the per-trial factor of the likelihood, ",
             "so sum(log(p)) is the CONDITIONAL data log-likelihood: ",
             "equal to logLik(fit) to machine precision with no random ",
             "effects, and different from it in a hierarchical fit, ",
             "where logLik() is the Laplace MARGINAL (measured: ",
             "-1507.8 against -1542.5 at sd(id) = 1.07). NOTE the ",
             "difference from frmtmb's own rw_delta example, which does ",
             "have fitted(): it codes a two-armed choice 0 and 1 and ",
             "returns P(arm 1), so y - p is the ordinary binary ",
             "residual. That does not survive a fourth arm, and one ",
             "engine over two and four options was judged worth more ",
             "than fitted() on the two-option half."))
  } else {
    r(nm, "fitted", "refused",
      paste0("Refused for a different reason from the option-code ",
             "families, and the difference is worth stating: this ",
             "response is a response TIME, which is ordered and does ",
             "have a conditional mean, so nothing about the model ",
             "forbids fitted(). What is missing is a route to it. The ",
             "mean of a Wiener first-passage time conditional on the ",
             "boundary reached is a closed form that belongs to ",
             "frmtmb.eam, and what that package exports to this one is ",
             "the DENSITY, so declaring fitted_mean here would mean ",
             "re-deriving the mean rather than importing it. That is a ",
             "second seam and it is recorded as one in ",
             "dev/learn2-findings.md rather than guessed at. WHAT IS ",
             "AVAILABLE MEANWHILE: frm_value_trace(), whose `dens` ",
             "column is the per-trial joint density of the boundary and ",
             "the time, and whose `drift_t` column is the drift rate ",
             "the value difference produced on that trial."))
  }
  r(nm, "predict", "conditional",
    paste0("type = 'link' works everywhere, on new data too, and is how ",
           "a fitted learning parameter is read: ",
           "predict(type = 'link', dpar = 'alpha'). type = 'response' is ",
           "refused, for the reason in the fitted row, and separately ",
           "would be refused on newdata because a trial's choice ",
           "probability is conditional on a history newdata does not ",
           "carry."))
  r(nm, "residuals", "refused",
    paste0("All four types refuse and no type is available. 'response' ",
           "and 'pearson' refuse first, on the missing mean; see the ",
           "fitted row. 'osa' refuses in this package's own words for a ",
           "reason that would apply even if a mean existed. 'deviance' ",
           "is the one that changed this round: the family now declares ",
           "loglik_row(), so each trial's own log-likelihood IS ",
           "available, and what is still missing is ",
           if (nominal) {
             "the SIGN, which needs a conditional mean to depart from"
           } else {
             paste0("the saturated comparison: a row here contributes a ",
                    "DENSITY, whose supremum over the parameters at a ",
                    "fixed response time is unbounded, so there is no ",
                    "constant to subtract")
           },
           ". What a fit of one of these families is checked with is ",
           "its fitted value trajectory against the observed choices, ",
           "which is frm_value_trace() and a plot rather than a ",
           "residual."))
  r(nm, "s()", "works",
    paste0("A smooth reaches any parameter of the family, because every ",
           "parameter is an ordinary distributional parameter with its ",
           "own linear predictor. A learning rate that drifts smoothly ",
           "over a session is s(trial) on alpha."))
  r(nm, "smooth", "works",
    paste0("The random-effect block a penalized smooth becomes, for the ",
           "same reason. The family imposes nothing on the block ",
           "structure."))
  r(nm, "us", "works",
    paste0("The correlated per-subject parameters these models are ",
           "usually fitted with, spelled (1 | p | id) across two or more ",
           "parameters. This is the standard hierarchical ",
           "parameterization in the literature and it is ordinary ",
           "frmtmb grammar here."))
  r(nm, "|ID|", "works", "Same thing named as a grammar feature.")
  r(nm, "prior", "works",
    paste0("set_prior() reaches every parameter of the family. A prior ",
           "on a learning rate's intercept is the usual remedy for a ",
           "subject whose choices are perfectly consistent, which is ",
           "where an unpenalized inverse temperature runs away."))
  r(nm, "autoscale", "untested",
    "Not exercised. The parameters are all of order 1 by construction.")
  r(nm, "nl", "untested",
    paste0("A nonlinear body on one of these parameters is not ",
           "exercised. Nothing about the recursion argues against it: ",
           "the family reads a linear predictor's VALUE and does not ",
           "care how it was built."))
  r(nm, "mi()", "untested",
    "Not exercised. A missing choice drops its row; see keep_na = FALSE.")
  invisible(NULL)
}

#' The compatibility rules for what this package supplies.
#'
#' Every row was run unless it says `untested`, which is the honest
#' third state: an absent guard and a passing guard look the same from
#' outside.
#'
#' @noRd
ln_compat_rules <- function() {
  b <- compat_rule_builder()
  r <- b$r
  fams <- c("bandit2arm_delta", "bandit2arm_dual", "prl_fictitious",
            "bandit4arm2_kalman_filter", "ts_par7", "igt_pvl_delta",
            "igt_orl", "rlddm")
  # rlddm() is the one family whose response is not an option code, so
  # it is the one whose missing mean has a different reason. Everything
  # else about the eight is the same, which is what one engine buys.
  for (nm in fams) ln_common_rules(r, nm, nominal = nm != "rlddm")
  # ndt_group() is frmtmb.eam's term and only the one family here with a
  # non-decision time reads it. The other seven are refused at frame
  # assembly rather than by their own declaration, because the refusal
  # belongs to the term's owner: a grouping no family read leaves no
  # per-group table behind, and frmtmb.eam's frame check refuses exactly
  # that. The row says so rather than leaving the pair untested.
  for (nm in setdiff(fams, "rlddm")) {
    r(nm, "ndt_group()", "refused",
      paste0("Refused at frame assembly, by frmtmb.eam's own check. The ",
             "term says which trials share a non-decision-time bound ",
             "and this family has no non-decision time: its trial ",
             "contributes a choice probability, not a response-time ",
             "density. rlddm() is the family in this package that reads ",
             "it."))
  }

  ## ---- what differs, family by family ------------------------------
  for (nm in setdiff(fams, c("ts_par7", "rlddm"))) {
    r(nm, "simulate", "works",
      paste0("The structured simulator walks each subject forward, ",
             "drawing a choice from the current value store and learning ",
             "from the payoff that choice earns. It is the SAME ",
             "recursion the likelihood tapes, at mode = 'simulate', so a ",
             "draw and the density that scores it cannot drift apart. ",
             "simulate(), posterior_predict() and frm_simulate() all ",
             "reach it. It is coherent because the payoff schedule of ",
             "every arm is fixed in advance, which is what the reward() ",
             "and payoff() terms carry. A schedule that does NOT carry ",
             "it is refused by name, on all eight families: where every ",
             "column of the term holds the same value on every trial, ",
             "the data records the received outcome alone, the fit is ",
             "still exactly right and the draw is not, because the ",
             "option a simulated subject picks has no payoff to read. ",
             "Measured on the two-armed design, a duplicated schedule ",
             "draws at exact chance where the real one reproduces the ",
             "data. frm_task_simulate() is checked the same way, which ",
             "is what covers rlddm() and ts_par7()."))
  }
  r("ts_par7", "simulate", "refused",
    paste0("REFUSED, and the only capability that differs between this ",
           "family and its siblings. One two-step trial's draw is three ",
           "numbers, the stage-one choice, the stage-two state the ",
           "environment answers with, and the stage-two choice; ",
           "simulate() returns one response vector and there is nowhere ",
           "to put the other two. Use frm_task_simulate(), which returns ",
           "whole data frames and which every other family in this ",
           "package supports as well."))
  r("bandit2arm_delta", "reward()", "works",
    paste0("Carries what each arm would have paid on this trial, in arm ",
           "order. The likelihood reads only the chosen arm's entry, so ",
           "data recording the received outcome alone can pass it twice ",
           "and the fit is exactly right; measured, the two spellings ",
           "give one log likelihood to printed precision. The second ",
           "column is what makes the SIMULATOR coherent, and passing it ",
           "twice is refused there rather than drawn from. The ",
           "arity-one spelling reward(pay) is not available: an ",
           "addition term's arity is fixed at registration and frmtmb's ",
           "parser refuses any other argument count, so one column has ",
           "to be written twice."))
  r("bandit2arm_dual", "reward()", "works", "As bandit2arm_delta().")
  r("prl_fictitious", "reward()", "works",
    paste0("As bandit2arm_delta(), and this row said the opposite from ",
           "0.1.0 until it was measured. It claimed both columns were ",
           "READ here because counterfactual updating moves the ",
           "unchosen option's value. The update forms the outcome as ",
           "c1 * reward1 + c2 * reward2 with the CHOSEN indicators and ",
           "then FLIPS ITS SIGN, so what moves the unchosen value is ",
           "the negative of the realized outcome and the second column ",
           "never enters. Measured: replacing the unchosen entries with ",
           "N(100, 50) noise leaves the log likelihood bitwise ",
           "unchanged. ?prl_fictitious said this correctly all along. ",
           "So the second column is carried for the SIMULATOR here too, ",
           "and a duplicated schedule is refused there."))
  r("bandit4arm2_kalman_filter", "payoff()", "works",
    "Four columns, one per arm, for the same reason reward() has two.")
  r("igt_pvl_delta", "payoff()", "works",
    "Four columns, one per deck.")
  r("ts_par7", "payoff()", "works",
    paste0("Four columns, one per stage-two option, indexed as ",
           "2 * (state - 1) + choice. The likelihood reads the realized ",
           "one; frm_task_simulate() reads whichever the drawn path ",
           "reaches, so four identical columns are refused there. ",
           "stage2() is NOT checked that way: its two columns are the ",
           "observed state and choice rather than a payoff schedule."))
  r("ts_par7", "stage2()", "works",
    paste0("Carries the observed stage-two state and choice, which are ",
           "data the likelihood conditions on. The transition ",
           "probability is a known constant of the task rather than a ",
           "parameter, so it contributes a constant to the log ",
           "likelihood and is dropped."))

  r("rlddm", "dec()", "works",
    paste0("The boundary the response reached, 0 for the lower and 1 ",
           "for the upper, which is frmtmb.eam's term and frmtmb.eam's ",
           "coding. Arm 1 is the lower boundary and arm 2 the upper, so ",
           "reward(pay1, pay2) is in the same arm order as every other ",
           "family here. The term is registered by frmtmb.eam, which ",
           "this package imports for the density, so it is present ",
           "whenever rlddm() is."))
  r("rlddm", "ndt_group()", "works",
    paste0("Which trials share a non-decision-time BOUND, and this is ",
           "the one family in this package that reads it. Without it ",
           "the bound is the whole data set's fastest response, which ",
           "for a hierarchical fit is the fastest of every learner: a ",
           "subject deviation on `ndt` is then a deviation on a ",
           "fraction of somebody else's floor, and the scale tier's 100 ",
           "by 200 design reached a maximum gradient of 6.36e+09 with a ",
           "Hessian that was not positive definite and four NaN ",
           "standard errors. With ndt_group(id) each row is bounded by ",
           "its own group's fastest response, `ndt` is a fraction of ",
           "that bound on a plain logit, and the density multiplies it ",
           "back out. The term, its coercion and the bound are ",
           "frmtmb.eam's, reached through ndt_bound() rather than ",
           "copied. max_ndt and ndt_group() together are refused: they ",
           "set one bound to two different things. ",
           "frmtmb.eam::ndt_time() reports the fitted non-decision time ",
           "in seconds under either parameterization, and ",
           "predict(dpar = 'ndt', type = 'response') reports the ",
           "FRACTION under this one."))
  r("rlddm", "reward()", "works",
    paste0("As bandit2arm_delta(). The learning rule is that family's ",
           "exactly; what differs is that the value difference drives a ",
           "drift rate rather than a softmax. simulate() refuses on ",
           "this family for its own reason, but frm_task_simulate() is ",
           "the route the refusal names, so the duplicated-schedule ",
           "check runs there."))
  r("rlddm", "simulate", "refused",
    paste0("REFUSED, for the reason ts_par7() refuses it. One trial's ",
           "draw is two numbers, the boundary reached and the time it ",
           "took, and simulate() returns one response vector. Drawing ",
           "the time alone against the observed choice is a draw from a ",
           "different model, not a cheaper version of this one. Use ",
           "frm_task_simulate(), which returns whole data frames with ",
           "both columns written."))
  r("igt_orl", "payoff()", "works",
    "Four columns, one per deck, as igt_pvl_delta().")

  ## ---- the two methods this package adds ---------------------------
  r("frm_value_trace", "fitted", "works",
    paste0("THE REPLACEMENT for fitted() in this package, and the ",
           "reason refusing fitted() costs nothing: the per-option value ",
           "estimates each choice was made on, the prediction error the ",
           "outcome produced, and the choice probability. These are the ",
           "quantities papers plot. It runs at the estimates, off the ",
           "tape, through the exported accessors eval_dpars() and ",
           "frame_block_of()."))
  r("frm_value_trace", "predict", "conditional",
    paste0("The two answer different questions and neither replaces the ",
           "other. predict(type = 'link', dpar = ) gives a fitted ",
           "PARAMETER, one number per row, and works on new data. The ",
           "trace gives the fitted STATE, which exists only where a ",
           "history exists, so it runs on the training data alone. ",
           "predict(type = 'response') is refused, and the trace's `p` ",
           "column is what it would have returned."))
  r("frm_task_simulate", "simulate", "works",
    paste0("The other simulation route. It takes parameters directly, ",
           "one value per subject on their natural scales, and returns ",
           "whole data frames; frm_simulate() goes through the fitted ",
           "grammar and returns response vectors. The two share the ",
           "recursion and are checked against each other."))
  b$rules()
}
