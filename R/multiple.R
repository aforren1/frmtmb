# Multiply-imputed data: fit per imputation, pool by Rubin's rules.

#' Fit a model across multiply imputed datasets
#'
#' The frequentist counterpart of brms's `brm_multiple()`: fits the
#' model on every imputed dataset and pools the fixed-effect estimates
#' by Rubin's rules with Barnard-Rubin adjusted degrees of freedom.
#' Accepts a plain list of data frames or a `mice::mids` object.
#'
#' Only the fixed effects are pooled (covariance parameters are
#' reported per fit through `$fits`). For missing-data mechanisms
#' beyond imputation, see the roadmap note on latent-variable `mi()`
#' in dev/feature-gaps.md.
#'
#' @param formula,... As in [frm()].
#' @param data A list of completed data frames, or a `mice::mids`.
#' @return A `frmtmb_multiple` object: `pooled` (the Rubin table) and
#'   `fits` (the per-imputation fits).
#' @examples
#' set.seed(8)
#' n <- 80
#' x <- rnorm(n)
#' y <- rnorm(n, 1 + 0.5 * x, 1)
#' x[sample(n, 15)] <- NA
#' imps <- lapply(1:3, function(i) {
#'   xi <- x
#'   xi[is.na(xi)] <- sample(x[!is.na(x)], sum(is.na(xi)), TRUE)
#'   data.frame(y = y, x = xi)
#' })
#' frm_multiple(bf(y ~ x) + gaussian(), data = imps)
#' @export
frm_multiple <- function(formula, data, ...) {
  if (inherits(data, "mids")) {
    if (!requireNamespace("mice", quietly = TRUE)) {
      stop("A mids object needs the 'mice' package", call. = FALSE)
    }
    data <- mice::complete(data, action = "all")
  }
  if (!is.list(data) || length(data) < 2L ||
      !all(vapply(data, is.data.frame, TRUE))) {
    stop("data must be a list of at least two data frames (or a ",
         "mice::mids object)", call. = FALSE)
  }
  fits <- lapply(data, function(d) frm(formula, data = d, ...))
  m <- length(fits)

  cf <- vapply(fits, function(f) {
    bd <- f$estimates$betad
    if (length(fx <- f$frame$betad_fixed_idx)) bd <- bd[-fx]
    c(f$estimates$beta, bd)
  }, numeric(length(estimated_coef_names(fits[[1]]))))
  rownames(cf) <- estimated_coef_names(fits[[1]])
  Vs <- lapply(fits, vcov)

  qbar <- rowMeans(cf)
  ubar <- diag(Reduce(`+`, Vs) / m)
  bvar <- apply(cf, 1, stats::var)
  tvar <- ubar + (1 + 1 / m) * bvar
  se <- sqrt(tvar)
  # Barnard-Rubin small-sample degrees of freedom
  lam <- (1 + 1 / m) * bvar / tvar
  df_old <- (m - 1) / pmax(lam, 1e-12)^2
  dfcom <- df.residual(fits[[1]])
  df_obs <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lam)
  df_br <- df_old * df_obs / (df_old + df_obs)
  tstat <- qbar / se
  tab <- data.frame(
    estimate = qbar, se = se, df = df_br, t = tstat,
    p = 2 * stats::pt(-abs(tstat), df_br), fmi = lam
  )
  structure(list(pooled = tab, fits = fits, m = m),
            class = "frmtmb_multiple")
}

#' @export
print.frmtmb_multiple <- function(x, digits = 4, ...) {
  cat("Pooled over", x$m, "imputations (Rubin's rules):\n\n")
  print(signif(x$pooled, digits))
  invisible(x)
}
