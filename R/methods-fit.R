# S3 accessor and print methods for frmtmb_fit.

#' Template-shaped Estimate / Std. Error lists. as.list(sdreport) fills the
#' full parameter template regardless of map or REML, which keeps indexing
#' by lp$idx valid; mapped-out entries get SE = NA.
#'
#' @noRd
par_est_se <- function(fit, vcov = NULL) {
  est <- fit$estimates
  if (!is.null(vcov)) {
    # a supplied covariance (vcov_cluster(), say) replaces the SEs of
    # the outer parameters only; the inner ones have none to replace
    rv <- resolve_vcov_arg(fit, vcov, "summary")
    om <- outer_par_map(fit)
    se <- lapply(est, function(v) rep(NA_real_, length(v)))
    sd_all <- sqrt(pmax(0, diag(rv$V)))
    for (cp in unique(om$comp)) {
      i <- om$comp == cp
      pos <- seq_along(est[[cp]])
      if (cp == "betad" && length(fit$frame$betad_fixed_idx)) {
        pos <- setdiff(pos, fit$frame$betad_fixed_idx)
      }
      se[[cp]][pos] <- sd_all[i]
    }
    return(list(est = est, se = se))
  }
  # A degenerate fit (no free outer parameters, no random effects) gives
  # sdreport an empty summary, and as.list() has no rows to reshape
  degenerate <- !length(fit$opt$par) && !length(fit$obj$env$random)
  se <- if (degenerate) {
    lapply(est, function(v) rep(NA_real_, length(v)))
  } else {
    as.list(sdr_of(fit), what = "Std. Error")
  }
  if (length(fit$frame$betad_fixed_idx)) {
    se$betad[fit$frame$betad_fixed_idx] <- NA_real_
  }
  list(est = est, se = se)
}

#' @export
print.frmtmb_fit <- function(x, ...) {
  require_fitted(x, "print()")
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
             if (!is.null(x$prior)) " (MAP)"), "\n")
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

#' The name a linear predictor's coefficient block gets in output. A
#' multivariate fit needs the response in the key to stay unambiguous;
#' a univariate fit shows the dpar alone.
#'
#' @noRd
coef_block_key <- function(fit, lp) {
  if (length(fit$spec$responses) > 1) {
    paste(lp$resp, lp$dpar, sep = "_")
  } else {
    lp$dpar
  }
}

# `vcov` takes a covariance over the whole outer parameter vector -
# vcov_cluster(full = TRUE), or a function of the fit returning one -
# and reports its standard errors in the coefficient table, with a t
# reference when the matrix carries degrees of freedom. Documented on
# the vcov_cluster() page; the variance components keep the
# model-based standard errors either way.
#' @export
summary.frmtmb_fit <- function(object, vcov = NULL, ...) {
  rdf <- NULL
  if (!is.null(vcov)) {
    # resolve once: `vcov` may be a function of the fit
    rv <- resolve_vcov_arg(object, vcov, "summary")
    vcov <- rv$V
    rdf <- rv$df
  }
  ps <- par_est_se(object, vcov)
  coefs <- list()
  for (lp in object$frame$linpreds) {
    est <- ps$est[[lp$par]][lp$idx]
    se <- ps$se[[lp$par]][lp$idx]
    z <- est / se
    cm <- if (is.null(rdf)) {
      cbind(Estimate = est, `Std. Error` = se, `z value` = z,
            `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
    } else {
      cbind(Estimate = est, `Std. Error` = se, `t value` = z,
            `Pr(>|t|)` = 2 * stats::pt(-abs(z), rdf))
    }
    rownames(cm) <- colnames(lp$X)
    coefs[[coef_block_key(object, lp)]] <- cm
  }
  structure(
    list(call = object$call, family = family(object),
         formula = formula(object), nobs = stats::nobs(object),
         ngrps = ngrps(object),
         loglik = logLik(object), AIC = stats::AIC(object),
         BIC = stats::BIC(object), REML = object$REML,
         coefficients = coefs, varcor = VarCorr(object),
         rescor = rescor_matrix(object),
         # R-side residual correlation, on the natural scale with the
         # same delta-method interval confint_varcorr() reports
         autocor = local({
           tr <- autocor_trans_rows(object)
           if (is.null(tr)) NULL else {
             z <- stats::qnorm(0.975)
             m <- cbind(
               Estimate = varcorr_untrans(tr$type, tr$est_t),
               `2.5 %` = varcorr_untrans(tr$type, tr$est_t - z * tr$se_t),
               `97.5 %` = varcorr_untrans(tr$type, tr$est_t + z * tr$se_t))
             rownames(m) <- if (length(unique(tr$block)) > 1L) {
               paste(tr$block, tr$term)
             } else tr$term
             attr(m, "label") <- tr$block[1L]
             m
           }
         }),
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
  if (length(x$ngrps %||% integer(0))) {
    cat("Groups:", paste(names(x$ngrps), x$ngrps, sep = ", ",
                         collapse = "; "), "\n")
  }
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
  if (!is.null(x$autocor)) {
    cat("\nWithin-group residual correlation: ",
        attr(x$autocor, "label"), "\n", sep = "")
    m_ac <- x$autocor
    attr(m_ac, "label") <- NULL
    print(signif(m_ac, 4))
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

#' Estimated outer parameters: profiled betas (control profile = TRUE)
#' leave opt$par but are still estimated and must count toward df.
#'
#' @noRd
n_outer_est <- function(object) {
  n <- length(object$opt$par)
  if (isTRUE(object$control$profile)) {
    n <- n + length(object$frame$par_template$beta)
  }
  n
}

#' @export
logLik.frmtmb_fit <- function(object, ...) {
  require_fitted(object, "logLik() (and AIC(), BIC(), anova())")
  structure(-object$opt$objective,
            df = n_outer_est(object),
            nobs = object$frame$n_obs,
            REML = object$REML,
            class = "logLik")
}

#' @export
nobs.frmtmb_fit <- function(object, ...) object$frame$n_obs

#' @export
df.residual.frmtmb_fit <- function(object, ...) {
  object$frame$n_obs - n_outer_est(object)
}

#' @export
family.frmtmb_fit <- function(object, ...) {
  fams <- lapply(object$spec$responses, `[[`, "family")
  if (length(fams) == 1) fams[[1]] else fams
}

#' Estimated residual correlation matrix (rescor fits), else NULL
#' @param fit A `frmtmb_fit`.
#' @return A correlation matrix or `NULL`.
#' @examples
#' set.seed(2)
#' n <- 80
#' dd <- data.frame(x = rnorm(n))
#' # two responses that share a residual disturbance
#' e <- rnorm(n)
#' dd$y1 <- 1 + 0.5 * dd$x + e + rnorm(n, 0, 0.5)
#' dd$y2 <- 2 - 0.3 * dd$x + e + rnorm(n, 0, 0.5)
#' fit <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x), rescor = TRUE) + gaussian(),
#'            data = dd)
#' rescor_matrix(fit)
#'
#' # a fit without rescor has no residual correlation to report
#' fit0 <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x)) + gaussian(), data = dd)
#' rescor_matrix(fit0)
#' @export
rescor_matrix <- function(fit) {
  if (inherits(fit, "frmtmb_draws")) {
    stop("rescor_matrix() reads the fitted point estimate, so it ",
         "takes the frmtmb_fit, not draws: rescor_matrix(ds$fit). For ",
         "the posterior of the correlation, subset_draws() on the ",
         "rescor columns of as_draws(ds)", call. = FALSE)
  }
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

#' Names of the estimated (non-mapped) fixed coefficients, in template order.
#'
#' @noRd
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
#' `full = TRUE` is the joint covariance of the whole outer parameter
#' vector on its internal scale: the fixed-effect coefficients, the
#' covariance parameters `theta` (log standard deviations, Fisher-z
#' correlations, and whatever else a structure keeps there), and any
#' extra parameters such as the ordinal thresholds. It is the matrix a
#' delta-method calculation on a variance component needs, and it is
#' what [hypothesis()] uses for `method = "wald"` - so an ICC or a
#' heritability is usually easier to ask for through `hypothesis()`,
#' which names the components for you, than to assemble by hand from
#' this matrix.
#'
#' Under `REML = TRUE` (or `frmtmb_control(profile = TRUE)`) the fixed
#' effects are integrated out of the outer problem, so they are not
#' part of `full = TRUE`; the block comes from the joint precision and
#' carries exactly the parameters [confint.frmtmb_fit()] reports.
#' `vcov(object)`
#' is still the fixed-effect covariance there.
#'
#' Passing `cluster` forwards to [vcov_cluster()] for the
#' cluster-robust (sandwich) covariance, in the `sandwich::vcovCL()`
#' spelling: `vcov(fit, cluster = ~ g, type = "CR1")`.
#'
#' @param object A `frmtmb_fit`.
#' @param full If `TRUE`, include covariance parameters (`theta`),
#'   named as in `confint()` (the glmmTMB `vcov(full = TRUE)`
#'   convention).
#' @param cluster Optional clustering factor. When given, the result is
#'   [vcov_cluster()]'s cluster-robust covariance instead of the
#'   model-based one.
#' @param type Small-sample correction for `cluster`, see
#'   [vcov_cluster()].
#' @param ... Unused.
#' @return A covariance matrix.
#' @seealso [confint_varcorr()] for natural-scale intervals on the same
#'   covariance parameters, and [hypothesis()] for delta-method tests of
#'   expressions in them.
#'
#' @srrstats {RE4.6} The variance-covariance matrix of the model
#'   parameters is returned by `vcov()`: the fixed-effect block by
#'   default, and the covariance parameters as well under `full = TRUE`,
#'   named as in `confint()`. It comes from the inverse observed
#'   information produced by `RTMB::sdreport()`, or from the joint
#'   precision for a REML or profiled fit. A covariance that could not be
#'   recovered from the Hessian warns rather than returning silent `NaN`.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#' # standard errors of the fixed effects
#' sqrt(diag(vcov(fit)))
#' # the covariance parameters join the block on their internal scale
#' rownames(vcov(fit, full = TRUE))
#'
#' # the matrix is what a delta-method calculation needs
#' V <- vcov(fit)
#' a <- c(1, 2)                       # prediction at x = 2, no group
#' sqrt(drop(t(a) %*% V[1:2, 1:2] %*% a))
#' @export
vcov.frmtmb_fit <- function(object, full = FALSE, cluster = NULL,
                            type = "CR0", ...) {
  if (!is.null(cluster)) {
    return(vcov_cluster(object, cluster, type = type, full = full))
  }
  nm <- estimated_coef_names(object)
  if (!object$REML && !isTRUE(object$control$profile)) {
    V <- sdr_of(object)$cov.fixed
    # This branch reads an already-inverted covariance and used to
    # return whatever sdreport put there. When sdreport could not
    # invert the outer Hessian that is a matrix of NaN, on a fit whose
    # own convergence checks all passed - so say so here, or nothing
    # does. (The REML/profile branch below warns through
    # solve_joint_precision().)
    if (any(!is.finite(V))) warn_nonfinite_cov(object$cache)
    if (full) {
      # cov.fixed rows repeat the component names; the per-parameter
      # names (confint rows) are the useful labels
      dimnames(V) <- list(outer_par_names(object),
                          outer_par_names(object))
      return(V)
    }
    ord <- c(which(rownames(V) == "beta"), which(rownames(V) == "betad"))
    V <- V[ord, ord, drop = FALSE]
  } else {
    Q <- sdr_of(object)$jointPrecision
    Vall <- solve_joint_precision(Q, object$cache)
    rn <- rownames(Q)
    if (full) {
      # The outer parameter vector under REML (or control profile =
      # TRUE) does not contain beta: it is integrated out. So
      # full = TRUE returns exactly the parameters confint() reports,
      # which is what the naming invariant asks for, and the block
      # comes out of the joint precision rather than cov.fixed. Use
      # vcov(object) for the fixed-effect covariance.
      comps <- setdiff(names(object$frame$par_template),
                       c("b", "miss", "beta"))
      keep <- unlist(lapply(comps, function(cp) which(rn == cp)))
      onm <- outer_par_names(object)
      if (length(keep) == length(onm)) {
        Vf <- as.matrix(Vall[keep, keep, drop = FALSE])
        dimnames(Vf) <- list(onm, onm)
        return(Vf)
      }
      warning("full = TRUE could not align the joint-precision blocks ",
              "with the outer parameter names; returning the ",
              "fixed-effect block", call. = FALSE)
    }
    ord <- c(which(rn == "beta"), which(rn == "betad"))
    V <- as.matrix(Vall[ord, ord, drop = FALSE])
  }
  dimnames(V) <- list(nm, nm)
  V
}

#' Per-group coefficients (fixed effects plus conditional modes)
#'
#' Follows the lme4/glmmTMB/brms convention: for each random-effect
#' grouping factor, the fixed effects of its linear predictor broadcast
#' over the group levels, with the conditional modes added to the
#' matching columns. Use [fixef()] for the fixed effects alone.
#'
#' The result is a list of data frames keyed by grouping factor. When
#' random effects appear in more than one dpar (or response), an outer
#' layer keyed like [fixef()] is added. Smooth terms are excluded. A fit
#' without random effects returns [fixef()] (the single coefficient
#' vector when there is one linear predictor).
#'
#' @param object A `frmtmb_fit`.
#' @param ... Unused.
#' @return A named list of data frames, one per grouping factor, each
#'   with one row per group level and one column per coefficient. When
#'   random effects appear in more than one linear predictor, the list is
#'   nested one level deeper, keyed as in [fixef()]. A fit without random
#'   effects returns the [fixef()] value instead.
#'
#' @srrstats {RE4.2} Model coefficients are returned by `coef()`, in the
#'   lme4 and glmmTMB sense of per-group coefficients (fixed effects
#'   broadcast over the group levels with the conditional modes added),
#'   and by [fixef()] for the fixed effects alone. [ranef()] returns the
#'   conditional modes and [VarCorr()] the variance components.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#' # one row per group: the fixed effects with the modes added in
#' head(coef(fit)$g)
#' # which is fixef() plus ranef(), the lme4 identity
#' all.equal(coef(fit)$g[["(Intercept)"]],
#'           fixef(fit)$mu[["(Intercept)"]] + ranef(fit)$g[, 1],
#'           check.attributes = FALSE)
#'
#' # without random effects there are no groups, so coef() is fixef()
#' coef(frm(bf(y ~ x) + gaussian(), data = dd))
#' @export
coef.frmtmb_fit <- function(object, ...) {
  fe <- fixef(object)
  cvec <- coef_b(object)
  out <- list()
  for (bk in object$frame$re_blocks) {
    if (bk$covstruct == "smooth") next
    bmat <- t(matrix(cvec[bk$c_idx], nrow = bk$dim))
    for (cp in bk$components) {
      lp <- object$frame$linpreds[[cp$lp_key]]
      key <- coef_block_key(object, lp)
      # a second term on the same factor adds its modes to the same
      # frame; the fixed effects are broadcast only once
      gname <- bk$group_name
      df <- out[[key]][[gname]]
      if (!is.null(df) && !identical(rownames(df), bk$levels)) {
        gname <- bk$term_label
        df <- out[[key]][[gname]]
      }
      if (is.null(df)) {
        fev <- fe[[key]]
        df <- as.data.frame(
          matrix(fev, nrow = bk$n_levels, ncol = length(fev),
                 byrow = TRUE, dimnames = list(bk$levels, names(fev))),
          optional = TRUE
        )
      }
      for (j in seq_len(cp$dim)) {
        cn <- cp$cnms[j]
        bv <- bmat[, cp$offset + j]
        if (cn %in% colnames(df)) {
          df[[cn]] <- df[[cn]] + bv
        } else {
          df[[cn]] <- bv
        }
      }
      out[[key]][[gname]] <- df
    }
  }
  if (!length(out)) {
    # GLM-style fits: the mu vector alone, unless another dpar is
    # actually modeled with covariates
    if (length(object$spec$responses) == 1L && "mu" %in% names(fe)) {
      aux <- Filter(function(lp) lp$dpar != "mu",
                    object$frame$linpreds)
      simple <- all(vapply(aux, function(lp) {
        !is.null(lp$constant) ||
          (ncol(lp$X) == 1L && identical(colnames(lp$X), "(Intercept)"))
      }, TRUE))
      if (simple) return(fe[["mu"]])
    }
    return(fe)
  }
  if (length(out) == 1L) out[[1L]] else out
}

#' Extract fixed effects
#' @param object A `frmtmb_fit`.
#' @param ... Unused.
#' @return A named list of coefficient vectors, one per dpar.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#'
#' # one entry per distributional parameter, each on its link scale
#' fit <- frm(bf(y ~ x + (1 | g), sigma ~ x) + gaussian(), data = dd)
#' fixef(fit)
#' exp(fixef(fit)$sigma[["(Intercept)"]])   # sigma is modeled on the log
#'
#' # flatten to the vector confint() and hypothesis() name their rows by
#' unlist(fixef(fit))
#' @export
fixef <- function(object, ...) UseMethod("fixef")

#' @rdname fixef
#' @exportS3Method nlme::fixef
#' @export
fixef.frmtmb_fit <- function(object, ...) {
  require_fitted(object, "fixef()")
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
#' @param condVar If `TRUE`, attach the conditional SDs of the modes
#'   (from the Laplace posterior) as a `"condSD"` attribute on each
#'   matrix, in matching layout.
#' @param ... Unused.
#' @return A named list of levels-by-coefficients matrices, one per
#'   random-effect term. `as.data.frame()` gives the long form (with a
#'   `condsd` column when `condVar = TRUE` was used).
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#' # one matrix per random-effect term, levels by coefficients
#' ranef(fit)
#'
#' # condVar adds the conditional SDs a caterpillar plot needs
#' re <- as.data.frame(ranef(fit, condVar = TRUE))
#' head(re)
#' with(re[order(re$condval), ],
#'      plot(condval, seq_along(condval), pch = 16,
#'           xlim = range(condval - 2 * condsd, condval + 2 * condsd),
#'           xlab = "conditional mode", ylab = "group"))
#' @export
ranef <- function(object, ...) UseMethod("ranef")

#' @rdname ranef
#' @exportS3Method nlme::ranef
#' @export
ranef.frmtmb_fit <- function(object, condVar = FALSE, ...) {
  require_fitted(object, "ranef()")
  check_flag(condVar, "condVar")
  cvec <- coef_b(object)
  csd <- NULL
  if (condVar) {
    sdr <- sdr_of(object)
    dcr <- sdr$diag.cov.random
    if (!is.null(dcr)) {
      csd <- sqrt(pmax(dcr[names(sdr$par.random) == "b"], 0))
    }
  }
  out <- list()
  for (bk in object$frame$re_blocks) {
    M <- t(matrix(cvec[bk$c_idx], nrow = bk$dim))
    dimnames(M) <- list(bk$levels, bk$cnms)
    if (!is.null(csd)) {
      # rr factors live in a different space than the displayed
      # coefficients; no conditional SDs for those blocks
      S <- if (bk$covstruct == "rr") {
        matrix(NA_real_, nrow(M), ncol(M))
      } else {
        t(matrix(csd[bk$b_idx], nrow = bk$dim))
      }
      dimnames(S) <- dimnames(M)
      attr(M, "condSD") <- S
    }
    # appended, then named: `out[[label]] <- M` would DROP a block whose
    # label repeats (an animal model's (1 | gr(id, cov = A)) and its
    # permanent-environment (1 | id) both deparse to "1 | id")
    out[[length(out) + 1L]] <- M
  }
  names(out) <- vapply(object$frame$re_blocks, `[[`, "", "term_label")
  structure(out, class = "ranef_frmtmb")
}

#' @export
print.ranef_frmtmb <- function(x, ...) {
  for (i in seq_along(x)) {           # by position: labels can repeat
    cat("$", names(x)[i], "\n", sep = "")
    print(`attr<-`(x[[i]], "condSD", NULL))
    cat("\n")
  }
  invisible(x)
}

#' @export
as.data.frame.ranef_frmtmb <- function(x, ...) {
  rows <- lapply(seq_along(x), function(i) {   # by position: see print()
    nm <- names(x)[i]
    M <- x[[i]]
    S <- attr(M, "condSD")
    lv <- rownames(M) %||% as.character(seq_len(nrow(M)))
    df <- data.frame(
      grp = nm,
      term = rep(colnames(M), each = nrow(M)),
      level = rep(lv, ncol(M)),
      condval = as.vector(M)
    )
    if (!is.null(S)) df$condsd <- as.vector(S)
    df
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @export
as.data.frame.VarCorr_frmtmb <- function(x, ...) {
  rows <- list()
  # by position, not by name: two blocks can share a term label (an
  # animal model's (1 | gr(id, cov = A)) and its permanent-environment
  # (1 | id) both deparse to "1 | id"), and x[[nm]] would then return
  # the first block once per duplicate name
  for (i in seq_along(x)) {
    nm <- names(x)[i]
    V <- x[[i]]
    sds <- sqrt(diag(V))
    cn <- colnames(V)
    for (i in seq_along(sds)) {
      rows[[length(rows) + 1L]] <- data.frame(
        grp = nm, var1 = cn[i], var2 = NA_character_,
        vcov = V[i, i], sdcor = sds[i]
      )
    }
    if (ncol(V) > 1L) {
      C <- stats::cov2cor(V)
      for (i in seq_len(ncol(V) - 1L)) {
        for (j in seq(i + 1L, ncol(V))) {
          rows[[length(rows) + 1L]] <- data.frame(
            grp = nm, var1 = cn[i], var2 = cn[j],
            vcov = V[i, j], sdcor = C[i, j]
          )
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Extract random-effect covariance matrices
#' @param x A `frmtmb_fit`.
#' @param ... Unused.
#' @return A named list of covariance matrices, one per random-effect
#'   term. The names are the term labels, which can repeat when two
#'   blocks deparse the same way (`(1 | gr(id, cov = A)) + (1 | id)`, the
#'   animal model's genetic and permanent-environment terms). Index by
#'   position, not by name, when that is possible in your model.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(200), g = factor(rep(1:20, 10)))
#' u <- cbind(rnorm(20, 0, 0.8), rnorm(20, 0, 0.4))
#' dd$y <- rnorm(200, 1 + 0.5 * dd$x + u[dd$g, 1] + u[dd$g, 2] * dd$x, 1)
#' fit <- frm(bf(y ~ x + (x | g)) + gaussian(), data = dd)
#'
#' # the printed form shows SDs and correlations, as lme4 does
#' VarCorr(fit)
#' # the stored value is the covariance matrix itself
#' VarCorr(fit)[["x | g"]]
#' # tidy shape for broom.mixed-style code
#' as.data.frame(VarCorr(fit))
#' @export
VarCorr <- function(x, ...) UseMethod("VarCorr")

#' @rdname VarCorr
#' @exportS3Method nlme::VarCorr
#' @export
VarCorr.frmtmb_fit <- function(x, ...) {
  require_fitted(x, "VarCorr()")
  th <- x$estimates$theta
  out <- lapply(x$frame$re_blocks, function(bk) {
    if (bk$covstruct == "smooth") {
      # one smoothing variance; the k x k identity blowup is noise
      matrix(exp(th[bk$theta_idx])^2, 1, 1,
             dimnames = list("sd(wiggle)", "sd(wiggle)"))
    } else if (bk$covstruct %in% c("gp", "hsgp")) {
      # marginal GP sd; the range lives in confint_varcorr
      matrix(exp(th[bk$theta_idx[1]])^2, 1, 1,
             dimnames = list("sd(gp)", "sd(gp)"))
    } else {
      V <- covstruct_registry[[bk$covstruct]]$vcov(th[bk$theta_idx], bk)
      if (is_student_block(bk)) {
        # On a gr(dist = "student") block this is the SCALE matrix, not
        # the covariance: brms names the same quantity `sd_<group>__...`
        # and frmtmb keeps that name, so the matrix is tagged instead of
        # silently converted. The variance is `nu/(nu-2)` times it.
        attr(V, "dist_nu") <- bk$dist_nu
      }
      V
    }
  })
  names(out) <- vapply(x$frame$re_blocks, `[[`, "", "term_label")
  structure(out, class = "VarCorr_frmtmb")
}

#' @export
print.VarCorr_frmtmb <- function(x, ...) {
  # by position: duplicate term labels are legal (see
  # as.data.frame.VarCorr_frmtmb), and name lookup would print the
  # first block once per duplicate and never print the others
  for (i in seq_along(x)) {
    nm <- names(x)[i]
    V <- x[[i]]
    sdv <- sqrt(diag(V))
    nu <- attr(V, "dist_nu")
    cat(" ", nm, "\n")
    # a t block's diagonal is the SCALE, so the column is not headed
    # Std.Dev.: sd = scale * sqrt(nu/(nu-2)), printed alongside
    tab <- if (is.null(nu)) {
      data.frame(Name = colnames(V), `Std.Dev.` = signif(sdv, 5),
                 check.names = FALSE)
    } else {
      data.frame(Name = colnames(V), Scale = signif(sdv, 5),
                 `Std.Dev.` = signif(sdv * sqrt(nu / (nu - 2)), 5),
                 check.names = FALSE)
    }
    if (ncol(V) > 1) {
      C <- stats::cov2cor(V)
      corr <- format(signif(C, 3))
      corr[upper.tri(corr, diag = TRUE)] <- ""
      tab <- cbind(tab, Corr = corr[, -ncol(corr), drop = FALSE])
    }
    print(tab, row.names = FALSE)
    if (!is.null(nu)) {
      cat("   Student-t latent, nu = ", format(nu),
          " (fixed); the stored matrix is the scale\n", sep = "")
    }
  }
  invisible(x)
}
