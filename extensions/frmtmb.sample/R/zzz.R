# What this package tells frmtmb about itself at load time.
#
# Both registrations use seams frmtmb exports for the purpose
# (?frmtmb::`frmtmb-sampling-api`). Registering from .onLoad() rather
# than at top level is what a contributor outside the package must do:
# by then every namespace is sealed, so the collation-order question
# that governs frmtmb's own in-package contributors does not arise.

#' @noRd
.onLoad <- function(libname, pkgname) {
  # expects =: hmm and lca belong to frmtmb.latent, which this package
  # only suggests, so the two rules naming them are a forward reference
  # rather than a misspelling. Declaring them is what says so; without
  # it frmtmb refuses the whole registration, and before frmtmb checked,
  # the two rules were dropped in silence whenever frmtmb.latent was not
  # loaded.
  frmtmb_register_compat(features = c(frm_sample = "method"),
                         rules = sample_compat_rules,
                         expects = c("hmm", "lca"))
  # This answers get_prior(route = "sample"), and nothing else: the
  # default route reports what frm() applies and reads no registry, so
  # loading this package cannot change the table a user was already
  # getting. Without a registration that route refuses rather than
  # reporting flat, which is why this is a load-time job and not a
  # first-call one.
  frmtmb_register_prior_defaults(function(spec, frame) {
    # default_priors_for() reads only $spec and $frame of the object it
    # is given, so the two the registry hands over are the whole of it
    default_priors_for(list(spec = spec, frame = frame))
  })
  invisible()
}

#' The compatibility rules that name `frm_sample()`.
#'
#' They travel with the feature rather than staying behind, and that is
#' the general rule: a pair rule belongs to whichever package makes the
#' pair possible. `hmm x frm_sample` and `lca x frm_sample` are here for
#' that reason even though `hmm()` and `lca()` are frmtmb.latent's - the
#' pair only exists when this package is loaded. They are declared in
#' `frmtmb_register_compat(expects =)`, because a rule naming a feature
#' the table does not have would otherwise dangle.
#'
#' @noRd
sample_compat_rules <- function() {
  b <- compat_rule_builder()
  r <- b$r
  r("rescor", "frm_sample", "works",
    "Verified: tmbstan samples the multivariate outer parameters, the residual correlation included.")
  r("mvbf", "frm_sample", "works",
    "Verified: tmbstan samples the multivariate outer parameters. A boundary variance component still warns about mode initialization, as it does for a univariate fit.")
  r("mixture", "frm_sample", "conditional",
    "Mixture posteriors are multimodal. Sample with init = \"random\" rather than the mode-anchored default.")
  r("hmm", "frm_sample", "works",
    "Verified by a short run. The posterior is multimodal in the state labels, as a mixture's is.")
  r("lca", "frm_sample", "works",
    "Verified by a tiny fit: the objective has no random effects, so tmbstan samples the gating coefficients and item logits directly. The posterior is label-invariant, so read it with the same caution as any mixture posterior.")
  # override: the jitter condition holds wherever frm_sample() is used,
  # so it outranks the permissive and untested blanket defaults.
  r("frm_sample", "*", "conditional",
    "Chains start jittered around the fitted mode. Use init_jitter to widen the spread, or init = \"random\" for a multimodal posterior.",
    override = TRUE)
  b$rules()
}
