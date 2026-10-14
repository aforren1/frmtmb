p <- "R/families-extreme.R"
x <- paste(readLines(p), collapse = "\n")
r <- function(old, new) {
  stopifnot(grepl(old, x, fixed = TRUE))
  x <<- sub(old, new, x, fixed = TRUE)
}
r("    type = \"continuous\",
    post = list(
      mean_fn = function(dpars, aterms) {
        gev_mean(dpars[[\"mu\"]], dpars[[\"sigma\"]], dpars[[\"xi\"]])",
"    type = \"continuous\",
    post = list(
      fit_check = gev_fit_check,
      mean_fn = function(dpars, aterms) {
        gev_mean(dpars[[\"mu\"]], dpars[[\"sigma\"]], dpars[[\"xi\"]])")
r("#' Generalized extreme value family, dpars `mu` (location), `sigma`",
"#' Warn when a fitted `xi` has run to -1, its edge.
#'
#' For `xi <= -1` the GEV likelihood has no maximum, and both links stop
#' the optimizer above -1: `log1p` by construction, `identity` by the
#' `NaN` it gives there. A fit that stops on the edge is the supremum of
#' a likelihood that rises toward it, not an interior maximum, and the
#' optimizer's own warning there, \"false convergence (8)\", does not say
#' which parameter or why. Measured on the review's cases
#' (dev/reviews/2026-10-13-famuni.md, fix 3 and fix 4): seeds 621 and 622
#' under the identity link end at xi = -1.0000, and seed 1402 under
#' `log1p` at -0.9999998. 1e-3 above -1 is the edge.
#'
#' @noRd
gev_fit_check <- function(fit, resp) {
  multi <- length(fit$spec$responses) > 1L
  xi <- tryCatch(
    frm_linpred(fit, dpar = \"xi\", type = \"response\",
                resp = if (multi) resp),
    error = function(e) NULL)
  if (!is.numeric(xi) || !length(xi) || !any(xi < -1 + 1e-3, na.rm = TRUE)) {
    return(invisible(NULL))
  }
  frm_warning(\"gen_extreme_value: xi \",
              if (multi) paste0(\"of response \", resp, \" \"),
              \"ran to \", format(min(xi, na.rm = TRUE), digits = 7),
              \", the edge at -1. Below -1 the likelihood has no maximum, \",
              \"so the fit stopped on that edge, where the log-likelihood \",
              \"is a supremum and the standard errors do not have their \",
              \"usual meaning. A shorter upper tail than the GEV allows is \",
              \"the usual cause; compare the fit with a bounded family\",
              call. = FALSE, class = \"frmtmb_gev_xi_edge\")
  invisible(NULL)
}

#' Generalized extreme value family, dpars `mu` (location), `sigma`")
writeLines(strsplit(x, "\n")[[1]], p)
