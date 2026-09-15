# BLOCKER E, step 2: the fix. The arm/disarm pair joins core's internal
# export list beside the hyp_* helpers frmtmb.sample already uses, and
# hypothesis.frmtmb_draws() arms the note the way core's methods do.
ROOT = "C:/Users/adf44/source/r/frmtmb-wt-generics/"


def edit(rel, pairs):
    p = ROOT + rel
    s = open(p, encoding="utf-8", newline="").read()
    for old, new in pairs:
        n = s.count(old)
        if n != 1:
            raise SystemExit("%s: pattern found %d times: %r"
                             % (rel, n, old[:70]))
        s = s.replace(old, new)
    open(p, "w", encoding="utf-8", newline="\n").write(s)
    print("edited", rel)


# R/sampling-api.R was applied by the first run of this script, which
# then stopped on a wrong pattern below; it is not re-applied.

edit("R/confint.R", [
("""#' Arm the shadowing note for one user-level call and return the state
#' to restore afterwards (nested calls therefore stay one-shot too).
#'
#' @noRd
hyp_shadow_arm <- function() {""",
"""#' Arm the shadowing note for one user-level call and return the state
#' to restore afterwards (nested calls therefore stay one-shot too).
#'
#' Exported for the methods that live in other packages, and documented
#' on `?frmtmb-sampling-api`: every `hypothesis()` method has to arm the
#' note itself, because the generic that dispatched to it may be
#' brms's. `frmtmb.sample`'s draws method lost the note when this
#' arming left core's generic and before this pair was exported.
#'
#' @noRd
hyp_shadow_arm <- function() {"""),
("""  # UseMethod(): anything this package puts in its own generic is
  # simply not run then. Measured, when it was in the generic:""",
"""  # UseMethod(): anything this package puts in its own generic is
  # simply not run then. Every hypothesis() method arms it the same
  # way, including frmtmb.sample's draws method, which is why the pair
  # is on the ?frmtmb-sampling-api export list. Measured, when it was
  # in the generic:"""),
])

edit("extensions/frmtmb.sample/R/methods-draws.R", [
("""hypothesis.frmtmb_draws <- function(x, hypothesis, alpha = 0.05,
                                    class = NULL, group = NULL, ...) {
  fit <- x$fit
""",
"""hypothesis.frmtmb_draws <- function(x, hypothesis, alpha = 0.05,
                                    class = NULL, group = NULL, ...) {
  # Arm core's reserved-name shadowing note for this call. It used to
  # be armed by core's GENERIC, and this method relied on that. The
  # exported `hypothesis` now resolves to brms's generic whenever brms
  # is loaded, so core moved the arming into its methods, and this one
  # lost the note in every session until it armed it too. Measured on
  # a covariate named `sigma`: 1 note on the fit, 0 on its draws.
  old <- hyp_shadow_arm()
  on.exit(hyp_shadow_disarm(old), add = TRUE)
  fit <- x$fit
"""),
])
