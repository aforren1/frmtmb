# Confidence intervals, convergence diagnostics, model comparison.

# Names of the outer (optimized) parameters, in obj$par order: template
# component order, minus `random` components, minus mapped entries.
outer_par_names <- function(fit) {
  tpl <- fit$frame$par_template
  random <- if (fit$REML) c("b", "beta") else "b"
  nm <- character(0)
  for (cp in names(tpl)) {
    if (cp %in% random) next
    v <- names(tpl[[cp]])
    if (is.null(v)) v <- paste0(cp, "_", seq_along(tpl[[cp]]))
    if (cp == "betad" && length(fit$frame$betad_fixed_idx)) {
      v <- v[-fit$frame$betad_fixed_idx]
    }
    nm <- c(nm, v)
  }
  nm
}

#' Confidence intervals for frmtmb fits
#'
#' Covariance parameters (`theta_*`) are reported on their internal
#' (unconstrained) scale.
#'
#' @param object A `frmtmb_fit`.
#' @param parm Parameter names (see `rownames` of the Wald result) or
#'   indices. Required for `"profile"` and `"uniroot"`; defaults to all
#'   parameters for `"wald"`.
#' @param level Confidence level.
#' @param method `"wald"` (fast, from the sdreport covariance),
#'   `"profile"` (likelihood profile via [TMB::tmbprofile()]), or
#'   `"uniroot"` (likelihood-root search via [TMB::tmbroot()]).
#' @param ... Passed to the TMB profiling functions.
#' @return A matrix with columns `lwr`, `upr`, `est`.
#' @export
confint.frmtmb_fit <- function(object, parm = NULL, level = 0.95,
                               method = c("wald", "profile", "uniroot"),
                               ...) {
  method <- match.arg(method)
  nm <- outer_par_names(object)
  est <- object$opt$par
  a <- (1 - level) / 2

  resolve_parm <- function(parm) {
    if (is.null(parm)) return(seq_along(nm))
    if (is.numeric(parm)) return(as.integer(parm))
    idx <- match(parm, nm)
    if (anyNA(idx)) {
      stop("Unknown parameter(s): ",
           paste(parm[is.na(idx)], collapse = ", "),
           ". Available: ", paste(nm, collapse = ", "), call. = FALSE)
    }
    idx
  }
  idx <- resolve_parm(parm)

  if (method == "wald") {
    se <- sqrt(diag(sdr_of(object)$cov.fixed))
    ci <- cbind(lwr = est + stats::qnorm(a) * se,
                upr = est + stats::qnorm(1 - a) * se,
                est = est)
    rownames(ci) <- nm
    return(ci[idx, , drop = FALSE])
  }

  if (is.null(parm)) {
    stop("`parm` is required for method = '", method, "'", call. = FALSE)
  }
  ci <- matrix(NA_real_, length(idx), 3,
               dimnames = list(nm[idx], c("lwr", "upr", "est")))
  for (k in seq_along(idx)) {
    i <- idx[k]
    if (method == "profile") {
      pr <- TMB::tmbprofile(object$obj, name = i, trace = FALSE, ...)
      ci[k, 1:2] <- unname(stats::confint(pr, level = level))
    } else {
      r <- TMB::tmbroot(object$obj, name = i,
                        target = 0.5 * stats::qchisq(level, df = 1), ...)
      ci[k, 1:2] <- unname(r)
    }
    ci[k, 3] <- est[i]
  }
  ci
}

#' Convergence diagnostics for a frmtmb fit
#'
#' @param fit A `frmtmb_fit`.
#' @param quiet If `TRUE`, return the diagnostics without printing.
#' @return Invisibly, a list of diagnostics.
#' @export
diagnose <- function(fit, quiet = FALSE) {
  stopifnot(inherits(fit, "frmtmb_fit"))
  nm <- outer_par_names(fit)
  gr <- drop(fit$obj$gr(fit$opt$par))
  se <- sqrt(diag(sdr_of(fit)$cov.fixed))
  ev <- tryCatch(eigen(sdr_of(fit)$cov.fixed, symmetric = TRUE,
                       only.values = TRUE)$values,
                 error = function(e) NULL)
  th <- fit$estimates$theta
  out <- list(
    convergence = fit$opt$convergence,
    message = fit$opt$message,
    max_grad = max(abs(gr)),
    worst_grad = nm[which.max(abs(gr))],
    pdHess = isTRUE(sdr_of(fit)$pdHess),
    bad_se = nm[!is.finite(se)],
    min_cov_eigenvalue = if (!is.null(ev)) min(ev),
    extreme_theta = which(abs(th) > 8)
  )
  if (!quiet) {
    cat("Optimizer convergence code:", out$convergence,
        if (!is.null(out$message)) paste0("(", out$message, ")"), "\n")
    cat("Max |gradient|:", format(out$max_grad, digits = 4),
        "at", out$worst_grad, "\n")
    cat("Hessian positive definite:", out$pdHess, "\n")
    if (length(out$bad_se)) {
      cat("Non-finite standard errors:",
          paste(out$bad_se, collapse = ", "), "\n")
    }
    if (length(out$extreme_theta)) {
      cat("Extreme covariance parameters (|theta| > 8) at index ",
          paste(out$extreme_theta, collapse = ", "),
          ": the fit is near-singular; consider simplifying the ",
          "random effects (e.g. diag() instead of a correlated term)\n",
          sep = "")
    }
    if (out$convergence == 0 && out$pdHess && !length(out$bad_se) &&
        out$max_grad < 1e-3) {
      cat("No convergence problems detected\n")
    }
  }
  invisible(out)
}

#' Likelihood-ratio tests between nested frmtmb fits
#'
#' @param object A `frmtmb_fit`.
#' @param ... Further `frmtmb_fit` objects, nested with `object`.
#' @return An `anova` table.
#' @export
anova.frmtmb_fit <- function(object, ...) {
  fits <- c(list(object), Filter(function(x) inherits(x, "frmtmb_fit"),
                                 list(...)))
  if (length(fits) < 2) {
    stop("anova() needs at least two frmtmb fits to compare", call. = FALSE)
  }
  if (any(vapply(fits, `[[`, TRUE, "REML"))) {
    stop("Likelihood-ratio tests require ML fits (REML = FALSE)",
         call. = FALSE)
  }
  ll <- vapply(fits, function(f) as.numeric(logLik(f)), 0)
  df <- vapply(fits, function(f) attr(logLik(f), "df"), 0L)
  ord <- order(df)
  fits <- fits[ord]; ll <- ll[ord]; df <- df[ord]
  chisq <- c(NA, 2 * diff(ll))
  ddf <- c(NA, diff(df))
  p <- stats::pchisq(chisq, ddf, lower.tail = FALSE)
  tab <- data.frame(
    Df = df, logLik = ll, AIC = -2 * ll + 2 * df,
    Chisq = chisq, `Chi Df` = ddf, `Pr(>Chisq)` = p,
    check.names = FALSE
  )
  rownames(tab) <- vapply(fits, function(f) deparse1(formula(f)), "")
  structure(tab, class = c("anova", "data.frame"),
            heading = "Likelihood-ratio tests\n")
}

#' @export
update.frmtmb_fit <- function(object, ..., evaluate = TRUE) {
  cl <- object$call
  extras <- match.call(expand.dots = FALSE)$...
  for (nm in names(extras)) cl[[nm]] <- extras[[nm]]
  if (evaluate) eval(cl, parent.frame()) else cl
}
