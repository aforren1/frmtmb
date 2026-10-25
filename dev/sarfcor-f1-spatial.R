# Fix round 1, frmtmb.spatial's sar(): one sigma, a matrix of
# innovations, the bounds slot, the complete-graph warning, and the
# docs for newdata and the rho interval. Rscript <this>
p <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/extensions/frmtmb.spatial/R/sar.R"
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1("  # sigma is constant here (sigma_varies = FALSE), but it may come as one
  # value per row
  lsig <- if (length(sigma) == 1L) length(y) * log(sigma) else sum(log(sigma))
  list(u = r / sigma, logjac = sum(log1p(-rho * aux[[\"ev\"]])) - lsig)",
"  # one sigma: sigma_varies = FALSE
  list(u = r / sigma,
       logjac = sum(log1p(-rho * aux[[\"ev\"]])) - length(y) * log(sigma))")
rep1("#' The response at innovations `u`: `A^-1 (mu + sigma u)` (lag) or
#' `mu + A^-1 sigma u` (error). On numbers.
#'
#' @noRd
sar_unwhiten <- function(u, mu, sigma, theta, aux) {
  rho <- sar_rho(theta, aux)
  A <- Matrix::Diagonal(nrow(aux[[\"W\"]])) - rho * aux[[\"W\"]]
  if (identical(aux[[\"type\"]], \"lag\")) {
    as.vector(Matrix::solve(A, mu + sigma * u))
  } else {
    mu + as.vector(Matrix::solve(A, sigma * u))
  }
}",
"#' The response at innovations `u`: `A^-1 (mu + sigma u)` (lag) or
#' `mu + A^-1 sigma u` (error). On numbers. `u` may be a matrix, one
#' column per set of innovations: one sparse factorization of `A` then
#' serves every column.
#'
#' @noRd
sar_unwhiten <- function(u, mu, sigma, theta, aux) {
  rho <- sar_rho(theta, aux)
  A <- Matrix::Diagonal(nrow(aux[[\"W\"]])) - rho * aux[[\"W\"]]
  out <- if (identical(aux[[\"type\"]], \"lag\")) {
    as.matrix(Matrix::solve(A, mu + sigma * u))
  } else {
    mu + as.matrix(Matrix::solve(A, sigma * u))
  }
  if (is.matrix(u)) out else as.vector(out)
}

#' The interval of rho, for brms's class prior on it and frmtmb.sample's
#' flat default (`frmtmb::frmtmb_autocor_pars()`).
#'
#' @noRd
sar_bounds <- function(aux) {
  cbind(lower = aux[[\"lo\"]], upper = aux[[\"hi\"]])
}")
rep1("    usage = \"sar(M, type = \\\"lag\\\")\", compat_note = sar_compat_note)",
"    bounds = sar_bounds, usage = \"sar(M, type = \\\"lag\\\")\",
    compat_note = sar_compat_note)")
# the complete graph: the extreme eigenvalue repeated in all but one
# direction
rep1("  list(W = M, ev = ev, lo = 1 / lo, hi = 1 / hi, type = spec[[\"type\"]])",
"  sar_warn_repeated(ev, tol)
  list(W = M, ev = ev, lo = 1 / lo, hi = 1 / hi, type = spec[[\"type\"]])")
rep1("#' rho from its unbounded parameter:",
"#' Warn when an extreme eigenvalue of `W` is repeated in all but one
#' direction, as on a complete graph. At that end of rho's interval
#' `I - rho W` keeps one direction, the constant one for a
#' row-standardized `W`, which an intercept spans: the residuals vanish
#' in every other direction and the likelihood of both types grows
#' without bound, so rho runs to the end. A warning and not a refusal:
#' a model without an intercept keeps its maximum, and `build` cannot
#' see the formula.
#'
#' @noRd
sar_warn_repeated <- function(ev, tol) {
  n <- length(ev)
  if (n < 3L) return(invisible(NULL))
  for (end in c(min(ev), max(ev))) {
    k <- sum(abs(ev - end) <= tol)
    if (k >= n - 1L) {
      frm_warning(\"sar(): the eigenvalue \", format(end, digits = 4),
                  \" of 'M' is repeated \", k, \" times in \", n, \" rows \",
                  \"(a complete graph, for one). At rho = \",
                  format(1 / end, digits = 4), \", the end of its \",
                  \"interval, I - rho M keeps one direction, and with an \",
                  \"intercept the likelihood has no maximum there: rho \",
                  \"runs to that end. brms bounds rho by the eigenvalues \",
                  \"themselves, which hides it. A weights matrix whose \",
                  \"rows differ avoids it\", call. = FALSE)
    }
  }
  invisible(NULL)
}

#' rho from its unbounded parameter:")
# docs
rep1("#' row-standardized, brms's interval reaches past
#' \\eqn{1 / \\lambda_{max}}, where the determinant changes sign and
#' brms's log-determinant is the log of a negative number.",
"#' row-standardized, brms's interval reaches past
#' \\eqn{1 / \\lambda_{max}}, where the determinant changes sign; there
#' brms's log-determinant is the log of a negative number and Stan
#' rejects the proposal, so its support is frmtmb's interval there.
#'
#' The consequence on a row-standardized `W`: when the data favor rho
#' below the smallest eigenvalue, frmtmb and brms fit different models.
#' On a row-standardized 8 x 8 queen lattice (eigenvalues -0.5006 to 1)
#' with data drawn at rho = -0.8, frmtmb's estimate is -1.103, with a
#' log-likelihood 3.74 above the best inside brms's bounds, and brms's
#' posterior piles at its bound, -0.5006 (median -0.443; a fifth of the
#' draws within 0.02 of it).
#'
#' A `W` whose extreme eigenvalue is repeated in all but one direction,
#' as on a complete graph, gives a likelihood with no maximum at that
#' end when the model has an intercept; frmtmb warns.")
rep1("#' model. [frmtmb::autocor_matrix()] returns \\eqn{(A'A)^{-1}}. A
#' prediction or simulation on `newdata` is refused, because `M` covers
#' the fitted rows only; brms refuses the same calls unless a matrix over
#' the new rows comes in its `newdata2`. [frmtmb::frm_linpred()] gives
#' the linear predictor \\eqn{\\eta} on new data.",
"#' model. [frmtmb::autocor_matrix()] returns \\eqn{(A'A)^{-1}}.
#' [fitted()], [predict()] and [simulate()] take `newdata` with as many
#' rows as the fit, read as the fitted locations in their order, as
#' brms takes it: a counterfactual such as `x + 1` everywhere gives the
#' lag model's spillover, \\eqn{A^{-1} \\eta_{new}}. Other `newdata` is
#' refused, because `M` covers the fitted rows only; brms refuses it
#' too unless a matrix over the new rows comes in its `newdata2`.
#' [frmtmb::frm_linpred()] gives the linear predictor \\eqn{\\eta} on any
#' new data.
#'
#' brms's class priors `lagsar` and `errorsar` work:
#' `set_prior(\"normal(0, 0.5)\", class = \"lagsar\")` is a density on rho,
#' carried to its internal parameter with the log-Jacobian of the map.
#' frmtmb.sample's default is brms's, flat on rho over its interval.")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
cat("ok\n")
