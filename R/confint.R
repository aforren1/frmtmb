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

#' Natural-scale confidence intervals for covariance parameters
#'
#' Wald intervals for random-effect standard deviations (on the log
#' scale, back-transformed) and correlations (on the Fisher-z scale,
#' back-transformed), delta-method-propagated from the internal `theta`
#' covariance. One row per SD and per correlation of every block.
#'
#' @param fit A `frmtmb_fit`.
#' @param level Confidence level.
#' @return A data frame with columns `block`, `term`, `type`,
#'   `estimate`, `lwr`, `upr`.
#' @export
confint_varcorr <- function(fit, level = 0.95) {
  sdr <- sdr_of(fit)
  Vfull <- sdr$cov.fixed
  th_pos <- which(rownames(Vfull) == "theta")
  th <- fit$estimates$theta
  z <- stats::qnorm(1 - (1 - level) / 2)
  rows <- list()
  for (bk in fit$frame$re_blocks) {
    Vth <- Vfull[th_pos[bk$theta_idx], th_pos[bk$theta_idx],
                 drop = FALSE]
    t0 <- th[bk$theta_idx]
    if (bk$covstruct == "smooth") {
      se1 <- sqrt(Vth[1, 1])
      rows[[length(rows) + 1L]] <- data.frame(
        block = bk$term_label, term = "sd(wiggle)", type = "sd",
        estimate = exp(t0[1]),
        lwr = exp(t0[1] - z * se1), upr = exp(t0[1] + z * se1)
      )
      next
    }
    # g(theta): log-sds then atanh-correlations, via the block's vcov
    gfun <- function(tt) {
      V <- covstruct_registry[[bk$covstruct]]$vcov(tt, bk)
      sds <- sqrt(diag(V))
      out <- log(sds)
      if (nrow(V) > 1) {
        C <- stats::cov2cor(V)
        out <- c(out, atanh(pmin(pmax(C[lower.tri(C)], -0.9999), 0.9999)))
      }
      out
    }
    g0 <- gfun(t0)
    # numeric jacobian, central differences
    J <- vapply(seq_along(t0), function(i) {
      h <- 1e-5 * max(abs(t0[i]), 1)
      tp <- t0; tp[i] <- tp[i] + h
      tm <- t0; tm[i] <- tm[i] - h
      (gfun(tp) - gfun(tm)) / (2 * h)
    }, numeric(length(g0)))
    J <- matrix(J, nrow = length(g0))
    se_g <- sqrt(pmax(diag(J %*% Vth %*% t(J)), 0))

    d <- bk$dim
    sd_names <- if (bk$covstruct == "smooth") "sd(wiggle)" else bk$cnms
    n_sd <- length(g0) - if (d > 1) d * (d - 1) / 2 else 0
    for (i in seq_len(n_sd)) {
      rows[[length(rows) + 1L]] <- data.frame(
        block = bk$term_label,
        term = if (bk$covstruct == "smooth") "sd(wiggle)" else
          sd_names[min(i, length(sd_names))],
        type = "sd",
        estimate = exp(g0[i]),
        lwr = exp(g0[i] - z * se_g[i]),
        upr = exp(g0[i] + z * se_g[i])
      )
    }
    if (d > 1 && bk$covstruct != "smooth") {
      pairs <- which(lower.tri(diag(d)), arr.ind = TRUE)
      for (k in seq_len(nrow(pairs))) {
        i <- n_sd + k
        rows[[length(rows) + 1L]] <- data.frame(
          block = bk$term_label,
          term = paste0("cor(", bk$cnms[pairs[k, 2]], ",",
                        bk$cnms[pairs[k, 1]], ")"),
          type = "cor",
          estimate = tanh(g0[i]),
          lwr = tanh(g0[i] - z * se_g[i]),
          upr = tanh(g0[i] + z * se_g[i])
        )
      }
    }
  }
  do.call(rbind, rows)
}

# Effective degrees of freedom of the smooth blocks: for an iid wiggly
# block, edf = k - tr(posterior cov)/prior variance (the ridge identity).
smooth_edf <- function(fit) {
  blocks <- Filter(function(bk) bk$covstruct == "smooth",
                   fit$frame$re_blocks)
  if (!length(blocks)) return(NULL)
  sdr <- sdr_of(fit)
  dcr <- sdr$diag.cov.random
  if (is.null(dcr)) return(NULL)
  # par.random holds the `random` components in template order; b entries
  # are the ones named "b"
  b_pos <- which(names(sdr$par.random) == "b")
  th <- fit$estimates$theta
  out <- vapply(blocks, function(bk) {
    prior_var <- exp(th[bk$theta_idx])^2
    k <- bk$dim
    # +1 null-space columns live in beta; conventionally reported as the
    # penalized-part edf
    k - sum(dcr[b_pos[bk$b_idx]]) / prior_var
  }, numeric(1))
  stats::setNames(out, vapply(blocks, `[[`, "", "term_label"))
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
