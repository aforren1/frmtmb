#' Load hook.
#'
#' Two things this package must tell frmtmb about itself: a frame check
#' that refuses a dynamics input which is not constant inside a solve
#' group, and the one pair the package knows to be unusable. Registering
#' at load time is what keeps frmtmb free of any mention of ODEs.
#'
#' @noRd
.onLoad <- function(libname, pkgname) {
  frmtmb::frmtmb_register_frame_check(check_ode_constancy)
  # expects =: frm_sample is frmtmb.sample's feature, and this package
  # neither depends on nor suggests it. Declaring the forward reference
  # is what keeps the rule from dangling; without it frmtmb refuses the
  # whole registration.
  frmtmb::frmtmb_register_compat(features = c("frm_ode()" = "special"),
                                 rules = ode_compat_rules,
                                 expects = "frm_sample")
  invisible()
}

#' The one compatibility rule this package carries.
#'
#' The pair belongs here rather than in frmtmb.sample, by the rule that
#' a pair rule lives with whichever package makes the pair possible.
#'
#' The evidence is in `dev/upstream/rtmbode-issues.md`, section 1, and
#' reproduces without frmtmb: `RTMBode::ode()` calls `deSolve::ode()`
#' with no guard, and deSolve answers an extreme parameter either with
#' an R error or with fewer rows than were asked for, which the adjoint
#' node reports as `Wrong output length`. Stan reaches such parameters
#' by construction, since its first warmup step is of size 1 on the
#' unconstrained scale, so the chain aborts deterministically even when
#' it starts at the fitted optimum. The abort then crosses Stan's C++
#' boundary and leaves rstan's nested autodiff arena unbalanced for the
#' rest of the session, so the cost of trying is not one failed call
#' but a session that has to be restarted. That is what makes a refusal
#' worth more here than an attempt.
#'
#' @noRd
ode_compat_rules <- function() {
  b <- frmtmb::compat_rule_builder()
  b$r("frm_ode()", "frm_sample", "refused",
      paste0("The chain aborts at warmup iteration 1 and takes the R ",
             "session's rstan state with it. RTMBode calls deSolve ",
             "unguarded, and a failed solve escapes as an R error ",
             "rather than as NaN, which is what an optimizer shortens ",
             "a step on and a sampler rejects a proposal on; Stan's ",
             "first unconstrained step of size 1 reaches such a ",
             "parameter every time. Fit by maximum likelihood and read ",
             "the Wald intervals, or apply the patch series in this ",
             "package's dev/upstream/ to RTMBode. See ?frm_ode, ",
             "section 'Sampling an ODE fit'."))
  b$rules()
}
