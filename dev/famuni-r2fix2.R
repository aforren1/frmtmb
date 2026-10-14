p <- "R/families-extreme.R"
x <- paste(readLines(p), collapse = "\n")
i0 <- regexpr("#' Warn when a fitted `xi` has run to -1, its edge.", x, fixed = TRUE)
i1 <- regexpr("  invisible(NULL)\n}", substring(x, i0), fixed = TRUE) + i0 + nchar("  invisible(NULL)\n}") - 1L
stopifnot(i0 > 0, i1 > i0)
new <- r"--(#' Warn when a fitted `xi` has run to -1, its edge.
#'
#' For `xi <= -1` the GEV likelihood has no maximum, and both links stop
#' the optimizer above -1: `log1p` by construction, `identity` by the
#' `NaN` it gives there. A fit that stops on the edge is the supremum of
#' a likelihood that rises toward it, not an interior maximum, and the
#' optimizer's own warning there, "false convergence (8)", does not say
#' which parameter or why. Measured on the review's cases
#' (dev/reviews/2026-10-13-famuni.md, fix 3 and fix 4): seeds 621 and 622
#' under the identity link end at xi = -1.0000, and seed 1402 under
#' `log1p` at -0.9999998.
#'
#' A fitted `xi` within 1e-3 of -1 is not by itself an edge: the review's
#' 720 correct-model fits held one true maximum at -0.999092 (identity
#' link, n = 1000, seed 39), and the first version of this check called
#' it a supremum. So a scalar `xi` is judged by the sign of
#' `d logLik / d xi` at the stop, from the objective's gradient carried
#' to the `xi` scale: negative, the likelihood is still rising toward -1
#' and the warning says so; otherwise the fit is at a maximum inside the
#' band and gets a message instead. A predicted `xi` has no one
#' direction to read, and gets a warning that names the distance only.
#'
#' @noRd
gev_fit_check <- function(fit, resp) {
  multi <- length(fit$spec$responses) > 1L
  xi <- tryCatch(
    frm_linpred(fit, dpar = "xi", type = "response",
                resp = if (multi) resp),
    error = function(e) NULL)
  if (!is.numeric(xi) || !length(xi) || !any(xi < -1 + 1e-3, na.rm = TRUE)) {
    return(invisible(NULL))
  }
  of <- if (multi) paste0("of response ", resp, " ")
  at <- format(min(xi, na.rm = TRUE), digits = 7)
  lp <- fit$frame[["linpreds"]][[linpred_key(resp, "xi")]]
  if (!is.null(lp[["constant"]])) return(invisible(NULL))
  dll <- if (!dpar_has_predictor(fit$spec, lp)) gev_xi_slope(fit, lp)
  if (is.null(dll) || !is.finite(dll)) {
    frm_warning("gen_extreme_value: xi ", of, "reaches ", at, ", within ",
                "1e-3 of the edge at -1. Below -1 the likelihood has no ",
                "maximum, and near it the standard errors do not have ",
                "their usual meaning. Check whether the fit stopped on ",
                "that edge", call. = FALSE, class = "frmtmb_gev_xi_edge")
  } else if (dll < 0) {
    frm_warning("gen_extreme_value: xi ", of, "ran to ", at, ", the edge ",
                "at -1, with the likelihood still rising toward it (d ",
                "logLik / d xi = ", format(dll, digits = 3), "). Below -1 ",
                "the likelihood has no maximum, so the fit stopped on ",
                "that edge, where the log-likelihood is a supremum and ",
                "the standard errors do not have their usual meaning. A ",
                "shorter upper tail than the GEV allows is the usual ",
                "cause; compare the fit with a bounded family",
                call. = FALSE, class = "frmtmb_gev_xi_edge")
  } else {
    frm_message("gen_extreme_value: xi ", of, "is ", at, ", within 1e-3 ",
                "of -1, at a maximum (d logLik / d xi = ",
                format(dll, digits = 3), "). Below -0.5 the maximum is ",
                "not regular, so the standard errors do not have their ",
                "usual meaning", class = "frmtmb_gev_xi_near_edge")
  }
  invisible(NULL)
}

#' `d logLik / d xi` at the fit's stop for a scalar `xi`: the objective's
#' gradient in the `xi` coefficient, divided by the link's `dxi / deta`.
#' `NULL` where it cannot be read.
#'
#' @noRd
gev_xi_slope <- function(fit, lp) {
  p <- fit$opt$par
  g <- tryCatch(drop(fit$obj$gr(p)), error = function(e) NULL)
  if (is.null(g) || length(g) != length(p)) return(NULL)
  fx <- if (identical(lp[["par"]], "betad")) fit$frame[["betad_fixed_idx"]]
  free <- setdiff(seq_along(fit$frame[["par_template"]][[lp[["par"]]]]), fx)
  j <- which(names(p) == lp[["par"]])[match(lp[["idx"]], free)]
  if (length(j) != 1L || is.na(j)) return(NULL)
  eta <- fit$estimates[[lp[["par"]]]][lp[["idx"]]]
  # the objective is -logLik; a positive unit rescales, not the sign
  -g[j] / lp[["link"]]$mu_eta(eta)
})--"
x <- paste0(substring(x, 1, i0 - 1L), new, substring(x, i1 + 1L))
writeLines(strsplit(x, "\n")[[1]], p)
