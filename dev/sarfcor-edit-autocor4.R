p <- "R/autocor.R"
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1("#' The responses whose autocorrelation term is brms's default",
"#' A response's in-sample dpars with `mu` moved to the mean of the
#' response under its residual term: the `cov = FALSE` shift, or a
#' registry term's mean, `unwhiten()` at zero innovations, which for
#' brms's lag SAR is its `posterior_epred()`, `(I - rho W)^-1 mu`. The
#' draws keep `mu` itself (`autocor_cond_dpars()`), because a registry
#' term's draw applies its own map.
#'
#' @noRd
autocor_mean_dpars <- function(fit, resp, dp) {
  ac <- fit$frame[[\"autocor\"]][[resp]]
  if (autocor_is_seam(ac)) {
    dp[[\"mu\"]] <- autocor_seam_mean(ac, autocor_seam_theta(fit, ac),
                                    dp[[\"mu\"]], dp[[\"sigma\"]])
    return(dp)
  }
  autocor_cond_dpars(fit, resp, dp)
}

#' The responses whose autocorrelation term is brms's default")
rep1("autocor_cond_strip <- function(fit) {
  acs <- fit$frame[[\"autocor\"]]
  if (!length(acs)) return(fit)
  cond <- vapply(acs, autocor_is_cond, TRUE)",
"autocor_cond_strip <- function(fit) {
  acs <- fit$frame[[\"autocor\"]]
  if (!length(acs)) return(fit)
  # brms's incl_autocor = FALSE drops every autocorrelation term, so a
  # registry term goes too: its mean on a grid is mu, with no matrix
  # over the grid's rows
  cond <- vapply(acs, function(ac) {
    autocor_is_cond(ac) || autocor_is_seam(ac)
  }, TRUE)")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
