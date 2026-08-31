# S3 accessor and print methods for frmtmb_fit.

# Template-shaped Estimate / Std. Error lists. as.list(sdreport) fills the
# full parameter template regardless of map or REML, which keeps indexing
# by lp$idx valid; mapped-out entries get SE = NA.
par_est_se <- function(fit) {
  est <- fit$estimates
  se <- as.list(sdr_of(fit), what = "Std. Error")
  if (length(fit$frame$betad_fixed_idx)) {
    se$betad[fit$frame$betad_fixed_idx] <- NA_real_
  }
  list(est = est, se = se)
}

#' @export
print.frmtmb_fit <- function(x, ...) {
  if (inherits(x$bform, "frmtmb_mvformula")) {
    for (f in x$bform$forms) cat("frmtmb fit:", deparse1(f$formula), "\n")
  } else {
    cat("frmtmb fit:", deparse1(formula(x)), "\n")
  }
  fam_str <- paste(vapply(x$spec$responses,
                          function(r) r$family$family, ""),
                   collapse = ", ")
  cat("Family:", fam_str, "  Method:",
      paste0(if (x$REML) "REML" else "ML",
             if (!is.null(x$priors)) " (MAP)"), "\n")
  ll <- logLik(x)
  cat("logLik:", format(as.numeric(ll), digits = 6),
      " AIC:", format(stats::AIC(x), digits = 6),
      " nobs:", stats::nobs(x), "\n")
  cat("\nFixed effects:\n")
  for (nm in names(fixef(x))) {
    cat(" ", nm, ":\n", sep = "")
    print(fixef(x)[[nm]], digits = 4)
  }
  if (length(x$frame$re_blocks)) {
    cat("\nRandom effects:\n")
    print(VarCorr(x))
  }
  invisible(x)
}

coef_block_key <- function(fit, lp) {
  if (length(fit$spec$responses) > 1) {
    paste(lp$resp, lp$dpar, sep = "_")
  } else {
    lp$dpar
  }
}

#' @export
summary.frmtmb_fit <- function(object, ...) {
  ps <- par_est_se(object)
  coefs <- list()
  for (lp in object$frame$linpreds) {
    est <- ps$est[[lp$par]][lp$idx]
    se <- ps$se[[lp$par]][lp$idx]
    z <- est / se
    cm <- cbind(Estimate = est, `Std. Error` = se, `z value` = z,
                `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
    rownames(cm) <- colnames(lp$X)
    coefs[[coef_block_key(object, lp)]] <- cm
  }
  structure(
    list(call = object$call, family = family(object),
         formula = formula(object), nobs = stats::nobs(object),
         loglik = logLik(object), AIC = stats::AIC(object),
         BIC = stats::BIC(object), REML = object$REML,
         coefficients = coefs, varcor = VarCorr(object),
         rescor = rescor_matrix(object),
         smooth_edf = smooth_edf(object),
         extras = local({
           ex <- list()
           for (nm in object$frame$extra_names %||% character(0)) {
             cm <- cbind(Estimate = ps$est[[nm]],
                         `Std. Error` = ps$se[[nm]])
             rownames(cm) <- paste0(nm, "_", seq_len(nrow(cm)))
             ex[[nm]] <- cm
           }
           ex
         }),
         fixed_dpars = local({
           fx <- Filter(function(lp) !is.null(lp$constant),
                        object$frame$linpreds)
           stats::setNames(vapply(fx, `[[`, numeric(1), "constant"),
                           vapply(fx, function(lp) {
                             coef_block_key(object, lp)
                           }, character(1)))
         })),
    class = "summary.frmtmb_fit"
  )
}

#' @export
print.summary.frmtmb_fit <- function(x, ...) {
  cat("Family:", x$family$family, "\n")
  cat("Formula:", deparse1(x$formula), "\n")
  cat("Method:", if (x$REML) "REML" else "ML",
      "  nobs:", x$nobs, "\n")
  cat("logLik:", format(as.numeric(x$loglik), digits = 6),
      " AIC:", format(x$AIC, digits = 6),
      " BIC:", format(x$BIC, digits = 6), "\n")
  if (length(x$varcor)) {
    cat("\nRandom effects:\n")
    print(x$varcor)
  }
  if (!is.null(x$rescor)) {
    cat("\nResidual correlation:\n")
    print(signif(x$rescor, 4))
  }
  if (!is.null(x$smooth_edf)) {
    cat("\nSmooth terms (edf of the penalized part):\n")
    print(round(x$smooth_edf, 2))
  }
  for (nm in names(x$coefficients)) {
    if (nm %in% names(x$fixed_dpars)) next
    cat("\nCoefficients (", nm, "):\n", sep = "")
    stats::printCoefmat(x$coefficients[[nm]], signif.stars = FALSE,
                        na.print = "-")
  }
  for (i in seq_along(x$fixed_dpars)) {
    cat("\nFixed dpar: ", names(x$fixed_dpars)[i], " = ",
        x$fixed_dpars[i], "\n", sep = "")
  }
  for (nm in names(x$extras)) {
    cat("\nFamily parameters (", nm, ", internal scale):\n", sep = "")
    stats::printCoefmat(x$extras[[nm]], signif.stars = FALSE,
                        na.print = "-")
  }
  invisible(x)
}

#' @export
logLik.frmtmb_fit <- function(object, ...) {
  structure(-object$opt$objective,
            df = length(object$opt$par),
            nobs = object$frame$n_obs,
            REML = object$REML,
            class = "logLik")
}

#' @export
nobs.frmtmb_fit <- function(object, ...) object$frame$n_obs

#' @export
df.residual.frmtmb_fit <- function(object, ...) {
  object$frame$n_obs - length(object$opt$par)
}

#' @export
family.frmtmb_fit <- function(object, ...) {
  fams <- lapply(object$spec$responses, `[[`, "family")
  if (length(fams) == 1) fams[[1]] else fams
}

#' Estimated residual correlation matrix (rescor fits), else NULL
#' @param fit A `frmtmb_fit`.
#' @return A correlation matrix or `NULL`.
#' @export
rescor_matrix <- function(fit) {
  if (!isTRUE(fit$spec$rescor)) return(NULL)
  K <- length(fit$spec$responses)
  C <- us_chol_cor(fit$estimates$thetar, K)
  dimnames(C) <- list(names(fit$spec$responses),
                      names(fit$spec$responses))
  C
}

#' @export
formula.frmtmb_fit <- function(x, ...) {
  if (inherits(x$bform, "frmtmb_mvformula")) {
    x$bform$forms[[1]]$formula
  } else {
    x$bform$formula
  }
}

# Names of the estimated (non-mapped) fixed coefficients, in template order.
estimated_coef_names <- function(fit) {
  tpl <- fit$frame$par_template
  nm_betad <- names(tpl$betad)
  if (length(fit$frame$betad_fixed_idx)) {
    nm_betad <- nm_betad[-fit$frame$betad_fixed_idx]
  }
  c(names(tpl$beta), nm_betad)
}

#' Covariance matrix of the fixed-effect estimates
#'
#' Covers the estimated coefficients of every linear predictor; dpars
#' fixed to constants are excluded.
#'
#' @param object A `frmtmb_fit`.
#' @param full If `TRUE`, include covariance parameters (`theta`).
#' @param ... Unused.
#' @return A covariance matrix.
#' @export
vcov.frmtmb_fit <- function(object, full = FALSE, ...) {
  nm <- estimated_coef_names(object)
  if (!object$REML) {
    V <- sdr_of(object)$cov.fixed
    if (full) {
      return(V)
    }
    ord <- c(which(rownames(V) == "beta"), which(rownames(V) == "betad"))
    V <- V[ord, ord, drop = FALSE]
  } else {
    Q <- sdr_of(object)$jointPrecision
    Vall <- solve(Q)
    rn <- rownames(Q)
    ord <- c(which(rn == "beta"), which(rn == "betad"))
    V <- as.matrix(Vall[ord, ord, drop = FALSE])
    if (full) {
      warning("full = TRUE is not supported under REML; returning the ",
              "fixed-effect block", call. = FALSE)
    }
  }
  dimnames(V) <- list(nm, nm)
  V
}

#' @export
coef.frmtmb_fit <- function(object, ...) {
  unlist(fixef(object))
}

#' Extract fixed effects
#' @param object A `frmtmb_fit`.
#' @param ... Unused.
#' @return A named list of coefficient vectors, one per dpar.
#' @export
fixef <- function(object, ...) UseMethod("fixef")

#' @rdname fixef
#' @export
fixef.frmtmb_fit <- function(object, ...) {
  est <- object$estimates
  out <- list()
  for (lp in object$frame$linpreds) {
    v <- est[[lp$par]][lp$idx]
    names(v) <- colnames(lp$X)
    out[[coef_block_key(object, lp)]] <- v
  }
  out
}

#' Extract random-effect modes
#' @param object A `frmtmb_fit`.
#' @param ... Unused.
#' @return A named list of levels-by-coefficients matrices, one per
#'   random-effect term.
#' @export
ranef <- function(object, ...) UseMethod("ranef")

#' @rdname ranef
#' @export
ranef.frmtmb_fit <- function(object, ...) {
  b <- object$estimates[["b"]]
  out <- list()
  for (bk in object$frame$re_blocks) {
    M <- t(matrix(b[bk$b_idx], nrow = bk$dim))
    dimnames(M) <- list(bk$levels, bk$cnms)
    out[[bk$term_label]] <- M
  }
  out
}

#' Extract random-effect covariance matrices
#' @param x A `frmtmb_fit`.
#' @param ... Unused.
#' @return A named list of covariance matrices, one per random-effect term.
#' @export
VarCorr <- function(x, ...) UseMethod("VarCorr")

#' @rdname VarCorr
#' @export
VarCorr.frmtmb_fit <- function(x, ...) {
  th <- x$estimates$theta
  out <- lapply(x$frame$re_blocks, function(bk) {
    if (bk$covstruct == "smooth") {
      # one smoothing variance; the k x k identity blowup is noise
      matrix(exp(th[bk$theta_idx])^2, 1, 1,
             dimnames = list("sd(wiggle)", "sd(wiggle)"))
    } else {
      covstruct_registry[[bk$covstruct]]$vcov(th[bk$theta_idx], bk)
    }
  })
  names(out) <- vapply(x$frame$re_blocks, `[[`, "", "term_label")
  structure(out, class = "VarCorr_frmtmb")
}

#' @export
print.VarCorr_frmtmb <- function(x, ...) {
  for (nm in names(x)) {
    V <- x[[nm]]
    sdv <- sqrt(diag(V))
    cat(" ", nm, "\n")
    tab <- data.frame(Name = colnames(V), `Std.Dev.` = signif(sdv, 5),
                      check.names = FALSE)
    if (ncol(V) > 1) {
      C <- stats::cov2cor(V)
      corr <- format(signif(C, 3))
      corr[upper.tri(corr, diag = TRUE)] <- ""
      tab <- cbind(tab, Corr = corr[, -ncol(corr), drop = FALSE])
    }
    print(tab, row.names = FALSE)
  }
  invisible(x)
}
