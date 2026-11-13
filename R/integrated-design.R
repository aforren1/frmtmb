# Which rows of new data a prediction can compute without the values a
# Laplace approximation integrates out: read from the design, for the
# draws of frmtmb.sample's frm_sample(laplace = TRUE).

#' Rows of new data whose design reads no integrated value
#'
#' For each row of `newdata`, whether every entry of the design that
#' multiplies a random effect, a smooth or Gaussian-process coefficient
#' is exactly 0, in every linear predictor of the response: then a
#' prediction at that row reads none of them, whatever their values.
#' frmtmb.sample's draws from `frm_sample(laplace = TRUE)` hold none of
#' those values, and a prediction at such a row is the only one they can
#' answer that reads the group-level design. `(0 + x | g)` at `x = 0` is
#' one; `(1 + x | g)` at `x = -1` is not, though its two entries cancel
#' when the effects are equal.
#'
#' The answer is `FALSE` for every row wherever the design does not say:
#' a nonlinear body, a `me()` or `mi()` term, a latent residual or an
#' autocorrelation term on the response (the conditional mean reads the
#' fitted rows' residuals), a design that fails to build, or a design
#' whose rows are not `newdata`'s.
#'
#' @param fit A `frmtmb_fit` (or a fit at one posterior draw).
#' @param newdata A data frame of new rows.
#' @param resp The responses to read; `NULL` is every response.
#' @param allow_new_levels Passed to the design builder: a row at a level
#'   the fit did not see has a design too.
#' @return A logical vector, one element per row of `newdata`.
#' @keywords internal
#' @export
integrated_design_zero <- function(fit, newdata, resp = NULL,
                                   allow_new_levels = TRUE) {
  n <- NROW(newdata)
  none <- rep(FALSE, n)
  if (!is.data.frame(newdata) || !n) return(none)
  fr <- fit$frame
  resps <- resp %||% names(fit$spec$responses)
  if (length(fr[["mi_map"]] %||% list())) return(none)
  # a latent residual is a value per row that the design does not show
  keys <- paste0(resps, ".")
  for (bk in fr[["re_blocks"]] %||% list()) {
    if (!identical(bk[["covstruct"]], "latres")) next
    lk <- vapply(bk[["components"]], `[[`, "", "lp_key")
    if (any(outer(lk, keys, startsWith))) return(none)
  }
  for (r in resps) {
    if (!is.null((fr[["autocor"]] %||% list())[[r]])) return(none)
  }
  out <- rep(TRUE, n)
  for (lp in fr[["linpreds"]]) {
    if (!lp[["resp"]] %in% resps) next
    if (!is.null(lp[["nl_body"]]) || length(lp[["me"]]) ||
          length(lp[["mi"]])) {
      return(none)
    }
    pd <- tryCatch(suppressWarnings(
      pred_design(fit, lp, newdata, allow_new_levels, use_re = TRUE)),
      error = function(e) NULL)
    if (is.null(pd) || NROW(pd$X) != n) return(none)
    # a missing entry is not a 0: its product with an effect is unknown
    nonzero <- function(M) {
      M <- as.matrix(M)
      if (!ncol(M)) return(rep(FALSE, n))
      rowSums(is.na(M) | M != 0) > 0
    }
    for (rp in pd$re_parts) out[nonzero(rp$mm)] <- FALSE
    for (sp in pd$sm_parts) {
      out[nonzero(sp$Xr)] <- FALSE
      if (!is.null(sp$krig)) out[sp$krig$rows] <- FALSE
    }
  }
  out
}
