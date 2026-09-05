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
    "dec() coerces to a two-level 0/1 boundary and an lba response is a 1..n accumulator index, so a three-alternative model is refused by the coercion. The indicator arrives through vint() instead.")
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
  b$rules()
}
