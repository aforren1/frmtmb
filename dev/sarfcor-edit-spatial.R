# One-off source edit for lane sarfcor: register sar() from
# frmtmb.spatial's .onLoad() and add its compat rows.
setwd("C:/Users/adf44/source/r/frmtmb-wt-sarfcor/extensions/frmtmb.spatial")
rep1 <- function(p, old, new) {
  txt <- paste(readLines(p), collapse = "\n")
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop(p, " matches: ", n, " for ", substr(old, 1, 60))
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], p)
}
rep1("R/zzz.R", "  register_car()
  register_spde()
", "  register_car()
  register_spde()
  register_sar()
")
rep1("R/zzz.R", "# What this package tells frmtmb at load time: the structures it
# supplies, through frmtmb_register_covstruct(), and their rows in
# frm_compat(), through frmtmb_register_compat().",
"# What this package tells frmtmb at load time: the structures it
# supplies, through frmtmb_register_covstruct(), the residual
# correlation term sar(), through frmtmb_register_autocor(), and their
# rows in frm_compat(), through frmtmb_register_compat().")
rep1("R/zzz.R", "  r(\"spde\", \"predict\", \"conditional\",
    \"Locations outside the fitted set have no marginal variance of their own; allow_new_levels predicts them at the population level with no block contribution.\")
",
"  r(\"spde\", \"predict\", \"conditional\",
    \"Locations outside the fitted set have no marginal variance of their own; allow_new_levels predicts them at the population level with no block contribution.\")
  # sar(): frmtmb_register_autocor() writes the guards every residual
  # correlation term shares; these rows are the term's own
  r(\"sar()\", \"group:autocor_families\", \"works\",
    \"sar(M, type = \\\"lag\\\") is brms's SAR of the response, A y - mu ~ N(0, sigma^2 I) with A = I - rho W, whose mean is A^-1 mu; type = \\\"error\\\" is SAR of the residuals, A (y - mu) ~ N(0, sigma^2 I). student() gets the multivariate-t analog. The log-likelihood carries log|A| from the eigenvalues of W. rho is estimated on (1 / min eigenvalue, 1 / max eigenvalue), where |A| > 0; brms bounds it by (min eigenvalue, max eigenvalue). M is the weights matrix over every row of the data (a matrix, a sparse Matrix, or spdep's listw or nb), and the rows the model drops are dropped from it, as brms drops them. sigma and nu must be constant, as in brms.\",
    override = TRUE)
  r(\"sar()\", \"fitted\", \"conditional\",
    \"Under type = \\\"lag\\\" fitted() is (I - rho W)^-1 mu, brms's posterior_epred(); under \\\"error\\\" it is mu. newdata is refused, as brms refuses it without a matrix over the new rows (its newdata2); frm_linpred() gives mu on newdata.\")
  r(\"sar()\", \"predict\", \"conditional\",
    \"Draws in sample are one draw of the whole response through the model. newdata is refused, as in brms without newdata2.\")
  r(\"sar()\", \"simulate\", \"conditional\",
    \"simulate() draws the whole response, A^-1 (mu + e) or mu + A^-1 e. newdata is refused: M covers the fitted rows only.\")
  r(\"sar()\", \"residuals\", \"conditional\",
    \"\\\"response\\\" is y minus the mean above. \\\"pearson\\\" divides by the marginal sd, sigma * sqrt(diag((A'A)^-1)). They are not decorrelated. \\\"osa\\\" is refused.\")
")
rep1("R/frmtmb.spatial-package.R", "#' Through exported functions only. The structures are registered with
#' [frmtmb::frmtmb_register_covstruct()] and their rows in",
"#' Through exported functions only. The structures are registered with
#' [frmtmb::frmtmb_register_covstruct()], the residual correlation term
#' `sar()` with [frmtmb::frmtmb_register_autocor()], and their rows in")
rep1("R/frmtmb.spatial-package.R", "#' Without this package, a formula that writes `car()` or `spde()` is
#' refused by name, and the refusal says to load it.",
"#' Without this package, a formula that writes `car()`, `spde()` or
#' `sar()` is refused by name, and the refusal says to load it.")
rep1("R/frmtmb.spatial-package.R", "#' @importFrom frmtmb compat_rule_builder frm_message frm_stop
#'   frmtmb_register_compat frmtmb_register_covstruct",
"#' @importFrom frmtmb compat_rule_builder frm_message frm_stop
#'   frmtmb_register_compat frmtmb_register_covstruct
#'   frmtmb_register_autocor")
