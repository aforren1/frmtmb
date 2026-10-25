# What this package tells frmtmb at load time: the structures it
# supplies, through frmtmb_register_covstruct(), the residual
# correlation term sar(), through frmtmb_register_autocor(), and their
# rows in frm_compat(), through frmtmb_register_compat(). The structures
# are registered FIRST, because registering one adds it to the
# compatibility vocabulary that the rules then name.
#
# Every note below is core's note for the same pair before the move,
# copied from frmtmb's R/compat.R (dev/spatial-gen-zzz.R in the frmtmb
# repository writes this file), so frm_compat() reports the same rows
# with this package loaded as core reported with the structures in
# core. The rows for importance = and quadrature are not here:
# frmtmb_register_covstruct() writes them for every structure it
# registers, with the notes core's groups gave these two.

#' @noRd
.onLoad <- function(libname, pkgname) {
  register_car()
  register_spde()
  register_sar()
  frmtmb_register_compat(rules = spatial_compat_rules)
  invisible()
}

#' The compatibility rules for the structures this package supplies.
#'
#' @noRd
spatial_compat_rules <- function() {
  b <- compat_rule_builder()
  r <- b$r
  r("car", "*", "conditional",
    "car(M, gr = g, type = ) is a predictor special, not a bar term. M is a symmetric adjacency matrix with dimnames (rownames, colnames, or both, which then have to agree) covering every location; entries must be present and non-negative, and non-zero weights are binarized. type = \"escar\" is the proper CAR, \"esicar\" the intrinsic one under an EXACT sum-to-zero constraint, \"icar\" the same field under brms's soft constraint (con_sd), \"bym2\" the scaled mixture; escar needs every location to have a neighbor. con_sd changes nothing an esicar term reports: the constraint is exact and the delta method uses the centering projection, so the estimates, the prediction standard errors and the ranef() conditional SDs are all invariant to it. Every constrained type constrains each connected component, where brms constrains the global sum; the two agree on a connected graph. M belongs in data2 = list(M = M).",
    override = TRUE)
  r("car", "simulate", "works",
    "Draws come from the block's own fitted precision.")
  r("car", "predict", "conditional",
    "Locations outside the fitted set have no marginal variance of their own; allow_new_levels predicts them at the population level with no block contribution.")
  r("spde", "*", "conditional",
    "spde(fem, gr = node) is a predictor special taking a mesh's finite-element matrices (M0/M1/M2 or c0/g1/g2) as fixed data; gr maps observations onto mesh nodes BY ROW NUMBER (whole numbers in 1..nrow(M0), as integers or as a factor/character spelling of them), because the matrices carry no dimnames to match labels against. Unobserved nodes keep their column; a general projector matrix is not supported yet. The matrices belong in data2 = list(fem = fem).",
    override = TRUE)
  r("spde", "simulate", "works",
    "Draws come from the block's own fitted precision.")
  r("spde", "predict", "conditional",
    "Locations outside the fitted set have no marginal variance of their own; allow_new_levels predicts them at the population level with no block contribution.")
  # sar(): frmtmb_register_autocor() writes the guards every residual
  # correlation term shares and the families row, with this package's
  # note (register_sar()); these rows are the term's own
  r("sar()", "fitted", "conditional",
    "Under type = \"lag\" fitted() is (I - rho W)^-1 mu, brms's posterior_epred(); under \"error\" it is mu. newdata with one row per row of the term's matrix (the rows of the data, the ones the fit dropped for missing values among them) is read in their order, as brms reads it, so a counterfactual (x + 1 everywhere) works and a row with a missing predictor gives NA; other newdata is refused, as brms refuses it without a matrix over the new rows (its newdata2). On such newdata the lag mean is A^-1 mu_new. frm_linpred() gives mu on any newdata.")
  r("sar()", "predict", "conditional",
    "Draws are one draw of the whole response through the model. newdata with one row per row of the term's matrix (the rows of the data, the ones the fit dropped for missing values among them) is read in their order, as brms reads it, so a counterfactual (x + 1 everywhere) works and a row with a missing predictor gives NA; other newdata is refused, as brms refuses it without a matrix over the new rows (its newdata2).")
  r("sar()", "simulate", "conditional",
    "simulate() draws the whole response, A^-1 (mu + e) or mu + A^-1 e. newdata with one row per row of the term's matrix (the rows of the data, the ones the fit dropped for missing values among them) is read in their order, as brms reads it, so a counterfactual (x + 1 everywhere) works and a row with a missing predictor gives NA; other newdata is refused, as brms refuses it without a matrix over the new rows (its newdata2).")
  r("sar()", "residuals", "conditional",
    "\"response\" is y minus the mean above. \"pearson\" divides by the marginal sd, sigma * sqrt(diag((A'A)^-1)). They are not decorrelated. \"osa\" is refused.")
  b$rules()
}
