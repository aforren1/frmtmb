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
                 frm_value_trace = "method", frm_task_simulate = "method"),
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
ln_common_rules <- function(r, nm) {
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
  r(nm, "importance", "refused",
    paste0("REFUSED BY NAME, and the refusal is wider than the ",
           "mathematics. frm(importance =) reweights draws from the ",
           "Laplace Gaussian and needs one log-likelihood value per ",
           "GROUP; this family HAS that value, because its likelihood ",
           "factorizes over subjects and again over trials, but ",
           "frmtmb_structure(loglik =) returns one AD scalar for the ",
           "whole response and there is no slot to put the factors in. ",
           "THE SEAM: a structure slot carrying the finest factorization ",
           "the family has, per row where one exists and per group ",
           "otherwise, with `unit` left as the separate declaration of ",
           "the leave-one-out granularity. It is written up in frmtmb's ",
           "dev/rl-findings.md under 'Protocol seams'. This package ",
           "computes a per-row conditional likelihood already, for ",
           "frm_value_trace(), so it would fill a per-row slot the day ",
           "exists. Until then, the Laplace error is measured by ",
           "simulation instead: see vignette('learning')."))
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
           "need per-sequence ones here. Same seam as importance."))
  r(nm, "residuals_osa", "refused",
    paste0("Refused by name. One-step-ahead residuals re-tape the ",
           "objective with the response promoted to a parameter, and the ",
           "response here is a category the recursion selects with a ",
           "zero-one indicator rather than a continuous quantity the ",
           "tape could differentiate."))
  r(nm, "fitted", "refused",
    paste0("The family declares no fitted_mean, and that is a decision ",
           "rather than an omission. The response is the option a ",
           "subject took, coded 1 to K, and it is NOMINAL: arm 2 is not ",
           "twice arm 1, and the four decks of the Iowa gambling task ",
           "have no order at all. Core forms a fitted value as the ",
           "conditional mean of the response and a residual as ",
           "y - mean, so any mean this family supplied would be ",
           "arithmetic on a category code. core::cox() and ",
           "frmtmb.spline::royston_parmar() decline to invent a mean for ",
           "the same kind of reason. WHAT REPLACES IT: ",
           "frm_value_trace(), which returns more than a mean could. ",
           "Its `p` column is the per-trial factor of the likelihood, so ",
           "sum(log(p)) is the CONDITIONAL data log-likelihood: equal to ",
           "logLik(fit) to machine precision with no random effects, and ",
           "different from it in a hierarchical fit, where logLik() is ",
           "the Laplace MARGINAL (measured: -1507.8 against -1542.5 at ",
           "sd(id) = 1.07). NOTE the difference from frmtmb's own rw_delta ",
           "example, which does have fitted(): it codes a two-armed ",
           "choice 0 and 1 and returns P(arm 1), so y - p is the ",
           "ordinary binary residual. That does not survive a fourth ",
           "arm, and one engine over two and four options was judged ",
           "worth more than fitted() on the two-option half."))
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
           "fitted row. 'deviance' and 'osa' refuse in this package's own ",
           "words for reasons that would apply even if a mean existed. ",
           "What a fit of one of these families is checked with is its ",
           "fitted value trajectory against the observed choices, which ",
           "is frm_value_trace() and a plot rather than a residual."))
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
            "bandit4arm2_kalman_filter", "ts_par7", "igt_pvl_delta")
  for (nm in fams) ln_common_rules(r, nm)

  ## ---- what differs, family by family ------------------------------
  for (nm in setdiff(fams, "ts_par7")) {
    r(nm, "simulate", "works",
      paste0("The structured simulator walks each subject forward, ",
             "drawing a choice from the current value store and learning ",
             "from the payoff that choice earns. It is the SAME ",
             "recursion the likelihood tapes, at mode = 'simulate', so a ",
             "draw and the density that scores it cannot drift apart. ",
             "simulate(), posterior_predict() and frm_simulate() all ",
             "reach it. It is coherent because the payoff schedule of ",
             "every arm is fixed in advance, which is what the reward() ",
             "and payoff() terms carry."))
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
           "data recording the received outcome alone can pass it twice; ",
           "the second column is what makes the simulator coherent."))
  r("bandit2arm_dual", "reward()", "works", "As bandit2arm_delta().")
  r("prl_fictitious", "reward()", "works",
    paste0("Both columns are READ here rather than only carried: ",
           "counterfactual updating moves the unchosen option's value ",
           "too, so the second column enters the likelihood and not just ",
           "the simulator."))
  r("bandit4arm2_kalman_filter", "payoff()", "works",
    "Four columns, one per arm, for the same reason reward() has two.")
  r("igt_pvl_delta", "payoff()", "works",
    "Four columns, one per deck.")
  r("ts_par7", "payoff()", "works",
    paste0("Four columns, one per stage-two option, indexed as ",
           "2 * (state - 1) + choice. The likelihood reads the realized ",
           "one; frm_task_simulate() reads whichever the drawn path ",
           "reaches."))
  r("ts_par7", "stage2()", "works",
    paste0("Carries the observed stage-two state and choice, which are ",
           "data the likelihood conditions on. The transition ",
           "probability is a known constant of the task rather than a ",
           "parameter, so it contributes a constant to the log ",
           "likelihood and is dropped."))

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
