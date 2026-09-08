# What this package tells frmtmb about itself at load time.
#
# The registration uses the seam frmtmb exports for the purpose
# (?frmtmb::`frmtmb-sampling-api`, section "The compatibility
# registry"). Registering from .onLoad() rather than at top level is
# what a contributor outside the package must do: by then every
# namespace is sealed.
#
# Every name on the CORE side of a rule below was read off
# frm_compat_features() on frmtmb 0.53.0 rather than remembered, and
# spelled as the `name` column spells it: a callable feature carries its
# parentheses, a covariance structure and a method do not. The three
# names this package brings are declared through `features =`.

#' @noRd
.onLoad <- function(libname, pkgname) {
  frmtmb_register_compat(
    features = c(cross_wishart = "family", frm_coherence = "method",
                 frm_cross_simulate = "method"),
    rules = cp_compat_rules)
  invisible()
}

#' The compatibility rules for what this package supplies.
#'
#' Every row was run unless it says `untested`, which is the honest
#' third state: an absent guard and a passing guard look the same from
#' outside.
#'
#' @noRd
cp_compat_rules <- function() {
  b <- compat_rule_builder()
  r <- b$r

  ## ---- the family -------------------------------------------------
  r("cross_wishart", "vreal()", "works",
    "Required, and declared so. Three vreal columns carry the rest of the Hermitian matrix past the response slot: the second auto-spectrum and the two parts of the cross-spectrum. Without required_aterms an absent vreal reaches the density as NULL, the arithmetic gives numeric(0), and the fit returns a log-likelihood of zero; the declaration turns that into a refusal at frame assembly.")
  r("cross_wishart", "vint()", "works",
    "Required, and declared so. One vint column carries the degrees of freedom, which are per-row data rather than a parameter: they are decided by how the cross-spectrum was built and enter both the exponent on det(W) and the normalizing constant.")
  r("cross_wishart", "weights()", "works",
    "Verified: the weighted log likelihood equals the unweighted one at unit weights and doubles when every weight doubles. A weight of k on a row is not the same as k times the degrees of freedom, and nothing here pretends otherwise; use vint() to say a row rests on more draws.")
  r("cross_wishart", "cens()", "refused",
    "The family declares neither lcdf nor lccdf, and accepts_aterms does not list cens, so the term is refused by name at frame assembly. A censored cross-spectral matrix is not a thing anyone has: the CDF of a complex Wishart in two dimensions has no closed form worth taping, and there is no applied case asking for it.")
  r("cross_wishart", "trunc()", "refused",
    "Same reason, same mechanism.")
  r("cross_wishart", "s()", "works",
    "The case the package is for. s() on coh is a coherence spectrum, a smooth in frequency of the logit of coherence, and it cannot leave (0, 1) at any point on the smooth because the link is the constraint. s() on mu or pow2 is an ordinary log-spectrum smooth. Measured on a 40-subject fit: a coherence smooth over 60 frequencies recovers a two-peak truth.")
  r("cross_wishart", "smooth", "works",
    "The random-effect block a penalized smooth becomes, for the same reason. The family consumes one linear predictor per dpar and does not care how each was built.")
  r("cross_wishart", "predict", "works",
    "type = \"link\" on any of the four dpars, with se.fit, is what frm_coherence() and frm_phase() are built out of. type = \"response\" returns the mean of the RESPONSE, which is n * S11, the first auto-spectrum and not a coherence; read coherence with frm_coherence() instead.")
  r("cross_wishart", "fitted", "works",
    "post$mean_fn is supplied and is exact: W11 is Gamma(n, scale S11), so its mean is n * S11. It is the mean of the response column alone, which is the least interesting of the four numbers in the row.")
  r("cross_wishart", "residuals", "works",
    "All three types. \"pearson\" divides by sqrt(n) * S11, the exact standard deviation of W11. \"deviance\" uses the closed-form unit deviance against the saturated fit S = W/n, where tr(S^-1 W) is exactly 2n and log det S is log det W - 2 log n. Both speak about w11 only, which is a real limit rather than an approximation: the residual for a row's coherence has no place to be returned.")
  r("cross_wishart", "simulate", "refused",
    "Deliberate, through sim_refusal. A draw is a whole Hermitian matrix and the response slot carries one of its four numbers, so a simulated w11 sitting beside the OBSERVED w22, w12r and w12i is not a draw from anything, and it would usually not even be positive definite. frm_cross_simulate() returns all four columns.")
  r("cross_wishart", "residuals_osa", "untested",
    "One-step-ahead residuals re-tape the objective with the response promoted to a parameter. The density is branch-free and would tape, but promoting w11 alone while its three vreal companions stay fixed would break the positive definiteness the density needs, so this is the row most likely to be wrong.")
  r("cross_wishart", "REML", "untested",
    "mu is the primary dpar and would be integrated out. Not exercised.")
  r("cross_wishart", "mixture", "untested",
    "Not exercised. A two-component mixture on coherence is a plausible model of a channel pair that is coupled on some trials and not others.")
  r("cross_wishart", "mvbf", "untested",
    "Not exercised, and the natural way to reach more than two channels would be through it. See the family's refusal of p > 2.")
  r("cross_wishart", "quadrature", "untested", "Not exercised.")
  r("cross_wishart", "autoscale", "untested",
    "Not exercised. The four dpars are already on comparable scales because three of the four are on a log or logit link.")

  ## ---- the extractors ---------------------------------------------
  r("frm_coherence", "predict", "works",
    "It IS predict(type = \"link\", se.fit = TRUE) on the coh dpar, with the interval formed on the logit scale and pushed through plogis so that it cannot leave (0, 1). frm_phase() is the same call on phase, where no transform is needed and none is applied.")
  r("frm_coherence", "smooth", "works",
    "re.form is passed through, so re.form = NULL gives the per-group coherence with its shrinkage and re.form = NA gives the population one. Measured coverage of the population interval over 150 replicates: 0.912 at 4 segments per subject.")
  r("frm_coherence", "s()", "works",
    "A coherence spectrum is read by passing a frequency grid as newdata; the band comes from the same predict() call.")
  r("frm_coherence", "cross_wishart", "works",
    "The only family it applies to; every other is refused by name, because coh and phase are this family's dpars and a dpar of that name in another family would mean something else.")
  r("frm_cross_simulate", "cross_wishart", "works",
    "Verified against the density's own first moment: over 20000 draws the mean simulated matrix divided by n matches the fitted matrix to 0.0011 at n = 16 and 0.0112 at n = 2.")
  r("frm_cross_simulate", "simulate", "works",
    "It is the replacement for the refused method rather than a wrapper around it, and it takes nsim and seed with the same meaning.")
  b$rules()
}
