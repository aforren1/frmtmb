# Multiply-imputed data: fit per imputation, pool by Rubin's rules.

# Rubin's rules with Barnard-Rubin small-sample df. Q and U are
# quantity-by-imputation matrices of estimates and squared standard
# errors; dfcom is the complete-data residual df.
rubin_pool <- function(Q, U, dfcom) {
  m <- ncol(Q)
  qbar <- rowMeans(Q)
  ubar <- rowMeans(U)
  bvar <- apply(Q, 1, stats::var)
  tvar <- ubar + (1 + 1 / m) * bvar
  lam <- (1 + 1 / m) * bvar / tvar
  df_old <- (m - 1) / pmax(lam, 1e-12)^2
  df_obs <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lam)
  data.frame(estimate = qbar, se = sqrt(tvar),
             df = df_old * df_obs / (df_old + df_obs), fmi = lam)
}

#' Fit a model across multiply imputed datasets
#'
#' The frequentist counterpart of brms's `brm_multiple()`: fits the
#' model on every imputed dataset and pools by Rubin's rules with
#' Barnard-Rubin adjusted degrees of freedom. Accepts a plain list of
#' data frames or a `mice::mids` object.
#'
#' Fixed effects are pooled on the link scale, so distributional
#' coefficients like `sigma` appear in `$pooled` on their link (log)
#' scale, exactly as in `vcov()`. Random-effect SDs and correlations
#' are pooled on the scales where a Wald argument is defensible (log
#' for SDs and GP ranges, Fisher z for correlations, the
#' [confint_varcorr()] convention) and back-transformed, giving
#' `$pooled_varcorr`. For missing-data mechanisms beyond imputation,
#' see in-model `mi()`.
#'
#' @param formula,... As in [frm()].
#' @param data A list of completed data frames, or a `mice::mids`.
#' @param level Confidence level for the `$pooled_varcorr` interval.
#' @return A `frmtmb_multiple` object: `pooled` (the Rubin table for
#'   the fixed effects), `pooled_varcorr` (grp/term/type/estimate/
#'   lwr/upr/df/fmi for the random-effect SDs and correlations; `NULL`
#'   without random effects), and `fits` (the per-imputation fits).
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
frm_multiple <- function(formula, data, level = 0.95, ...) {
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
  dfcom <- df.residual(fits[[1]])

  nm <- estimated_coef_names(fits[[1]])
  cf <- vapply(fits, function(f) {
    bd <- f$estimates$betad
    if (length(fx <- f$frame$betad_fixed_idx)) bd <- bd[-fx]
    c(f$estimates$beta, bd)
  }, numeric(length(nm)))
  cf <- matrix(cf, nrow = length(nm), dimnames = list(nm, NULL))
  us <- vapply(fits, function(f) diag(vcov(f)), numeric(length(nm)))
  pl <- rubin_pool(cf, matrix(us, nrow = length(nm)), dfcom)
  tstat <- pl$estimate / pl$se
  tab <- data.frame(
    estimate = pl$estimate, se = pl$se, df = pl$df, t = tstat,
    p = 2 * stats::pt(-abs(tstat), pl$df), fmi = pl$fmi,
    row.names = nm
  )

  # variance components, pooled on the log/atanh scales and
  # back-transformed (interval via the Barnard-Rubin t quantile)
  vc <- lapply(fits, varcorr_trans_rows)
  pooled_varcorr <- NULL
  if (!is.null(vc[[1]])) {
    nq <- nrow(vc[[1]])
    Q <- matrix(vapply(vc, function(v) v$est_t, numeric(nq)), nrow = nq)
    U <- matrix(vapply(vc, function(v) v$se_t^2, numeric(nq)), nrow = nq)
    pv <- rubin_pool(Q, U, dfcom)
    tq <- stats::qt(1 - (1 - level) / 2, pv$df)
    pooled_varcorr <- data.frame(
      grp = vc[[1]]$block, term = vc[[1]]$term, type = vc[[1]]$type,
      estimate = varcorr_untrans(vc[[1]]$type, pv$estimate),
      lwr = varcorr_untrans(vc[[1]]$type, pv$estimate - tq * pv$se),
      upr = varcorr_untrans(vc[[1]]$type, pv$estimate + tq * pv$se),
      df = pv$df, fmi = pv$fmi
    )
  }
  structure(list(pooled = tab, pooled_varcorr = pooled_varcorr,
                 fits = fits, m = m, level = level),
            class = "frmtmb_multiple")
}

#' @export
print.frmtmb_multiple <- function(x, digits = 4, ...) {
  cat("Pooled over", x$m, "imputations (Rubin's rules):\n\n")
  print(signif(x$pooled, digits))
  if (!is.null(x$pooled_varcorr)) {
    cat("\nPooled variance components (log/Fisher-z scale, ",
        "back-transformed):\n\n", sep = "")
    vc <- x$pooled_varcorr
    vc[-(1:3)] <- lapply(vc[-(1:3)], signif, digits)
    print(vc, row.names = FALSE)
  }
  invisible(x)
}

#' @rdname hypothesis
#' @export
hypothesis.frmtmb_multiple <- function(x, hypothesis, alpha = 0.05, ...) {
  if (...length()) {
    warning("ignoring arguments unused by pooled hypothesis tests: ",
            paste(...names(), collapse = ", "), call. = FALSE)
  }
  exs <- lapply(hypothesis, hyp_parse)
  m <- x$m
  Q <- matrix(NA_real_, length(exs), m)
  U <- matrix(NA_real_, length(exs), m)
  for (j in seq_len(m)) {
    fit <- x$fits[[j]]
    pc <- hyp_par_cov(fit)
    for (i in seq_along(exs)) {
      ex <- exs[[i]]
      val <- hyp_eval(fit, ex, pc$vals, pc$comp)
      if (!is.numeric(val) || length(val) != 1L) {
        stop("Hypothesis '", hypothesis[i], "' must evaluate to a ",
             "single number", call. = FALSE)
      }
      g <- hyp_fd_grad(function(v) hyp_eval(fit, ex, v, pc$comp),
                       pc$vals)
      Q[i, j] <- val
      U[i, j] <- max(0, drop(t(g) %*% pc$V %*% g))
    }
  }
  pl <- rubin_pool(Q, U, df.residual(x$fits[[1]]))
  tq <- stats::qt(1 - alpha / 2, pl$df)
  tstat <- pl$estimate / pl$se
  out <- data.frame(
    hypothesis = hypothesis, estimate = pl$estimate, se = pl$se,
    lwr = pl$estimate - tq * pl$se, upr = pl$estimate + tq * pl$se,
    t = tstat, df = pl$df, p = 2 * stats::pt(-abs(tstat), pl$df)
  )
  rownames(out) <- NULL
  attr(out, "method") <- "wald"
  attr(out, "alpha") <- alpha
  class(out) <- c("frmtmb_hypothesis", "data.frame")
  out
}
