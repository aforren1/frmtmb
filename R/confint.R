# Confidence intervals, convergence diagnostics, model comparison.

#' One-line label for a fit including dpar formulas, so two models that
#' share a primary formula (e.g. plain vs distributional) stay
#' distinguishable in anova tables.
#'
#' @noRd
model_label <- function(fit) {
  one <- function(bform) {
    parts <- deparse1(bform$formula)
    for (nm in names(bform$pforms)) {
      parts <- c(parts, deparse1(bform$pforms[[nm]]))
    }
    for (nm in names(bform$pfix)) {
      parts <- c(parts, paste(nm, "=", bform$pfix[[nm]]))
    }
    paste(parts, collapse = ", ")
  }
  bf0 <- fit$bform
  if (inherits(bf0, "frmtmb_mvformula")) {
    paste(vapply(bf0$forms, one, ""), collapse = " + ")
  } else {
    one(bf0)
  }
}

#' Names of the outer (optimized) parameters, in obj$par order: template
#' component order, minus `random` components, minus mapped entries.
#'
#' @noRd
outer_par_names <- function(fit) {
  tpl <- fit$frame$par_template
  # mirror the MakeADFun random= construction in fit_assembled: b and
  # the mi() latent component are always inner, beta under REML or
  # control profile = TRUE
  random <- c("b", "miss")
  if (fit$REML || isTRUE(fit$control$profile)) {
    random <- c(random, "beta")
  }
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
#'   parameters for `"wald"` and `"boot"`.
#' @param level Confidence level.
#' @param method `"wald"` (fast, from the sdreport covariance),
#'   `"profile"` (likelihood profile via [TMB::tmbprofile()]),
#'   `"uniroot"` (likelihood-root search via [TMB::tmbroot()]), or
#'   `"boot"` (parametric-bootstrap percentile intervals through
#'   [frm_bootstrap()], the `lme4::confint(method = "boot")` analog;
#'   like the other methods it works on the internal parameter scale).
#'   `"Wald"` is accepted as an alias for `"wald"`.
#' @param nsim,seed Bootstrap draws and seed for `method = "boot"`.
#' @param ... Passed to the TMB profiling functions, or to
#'   [frm_bootstrap()] for `method = "boot"` (e.g. `re.form`).
#' @return A matrix with columns `lwr`, `upr`, `est`.
#'
#' @srrstats {RE4.3} Confidence intervals on the model coefficients are
#'   returned by `confint()`, by four methods: Wald from the sdreport
#'   covariance, likelihood profile, likelihood-root search, and
#'   parametric bootstrap. `confint_varcorr()` gives natural-scale
#'   intervals for standard deviations and correlations. Row names match
#'   those of `vcov(full = TRUE)`, which the test suite asserts.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#' # Wald intervals for every parameter, covariance ones included
#' confint(fit)
#'
#' # the likelihood profile does not assume a quadratic log-likelihood,
#' # so it is the one to trust for a variance component
#' confint(fit, parm = "theta_1", method = "profile")
#'
#' # confint_varcorr() puts the same information on the SD scale
#' confint_varcorr(fit)
#' \donttest{
#' # a parametric bootstrap, the lme4 confint(method = "boot") analog
#' confint(fit, parm = "x", method = "boot", nsim = 50, seed = 1)
#' }
#' @export
confint.frmtmb_fit <- function(object, parm = NULL, level = 0.95,
                               method = c("wald", "Wald", "profile",
                                          "uniroot", "boot"),
                               nsim = 500, seed = NULL, ...) {
  method <- match.arg(method)
  if (method == "Wald") method <- "wald"
  nm <- outer_par_names(object)
  est <- object$opt$par
  a <- (1 - level) / 2

  resolve_parm <- function(parm) {
    if (is.null(parm)) return(seq_along(nm))
    if (is.numeric(parm)) return(as.integer(parm))
    idx <- match(parm, nm)
    if (anyNA(idx)) {
      stop("Unknown parameter(s) in confint(parm =): ",
           paste(parm[is.na(idx)], collapse = ", "),
           ". Available: ", paste(nm, collapse = ", "), call. = FALSE)
    }
    idx
  }
  idx <- resolve_parm(parm)

  if (method == "wald") {
    sdr <- sdr_of(object)
    se <- sqrt(diag(sdr$cov.fixed))
    # profiled betas are absent from BOTH opt$par and par.fixed in
    # current RTMB; the fallback stays as a defensive alignment only
    if (length(est) != length(se)) est <- sdr$par.fixed
    ci <- cbind(lwr = est + stats::qnorm(a) * se,
                upr = est + stats::qnorm(1 - a) * se,
                est = est)
    rownames(ci) <- nm
    return(ci[idx, , drop = FALSE])
  }

  if (method == "boot") {
    # refits share the fit's control, so opt$par lines up with nm
    # (profile = TRUE excludes beta from both)
    bs <- frm_bootstrap(object, FUN = function(f) f$opt$par,
                        nsim = nsim, seed = seed, ...)
    ci <- cbind(lwr = apply(bs$t, 2, stats::quantile, a, na.rm = TRUE),
                upr = apply(bs$t, 2, stats::quantile, 1 - a,
                            na.rm = TRUE),
                est = est)
    rownames(ci) <- nm
    return(ci[idx, , drop = FALSE])
  }

  if (isTRUE(object$control$profile)) {
    stop("confint(method = '", method, "') needs a fit without ",
         "frmtmb_control(profile = TRUE)", call. = FALSE)
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

#' Transformed-scale Wald rows for the covariance parameters of one fit:
#' one row per SD/range (log scale) and per correlation (Fisher-z
#' scale), with the delta-method se on that scale. These scales are
#' where a normal approximation is defensible, which makes the rows the
#' right currency both for confint_varcorr's intervals and for Rubin
#' pooling across imputations in frm_multiple.
#'
#' @noRd
varcorr_trans_rows <- function(fit) {
  sdr <- sdr_of(fit)
  Vfull <- sdr$cov.fixed
  th_pos <- which(rownames(Vfull) == "theta")
  th <- fit$estimates$theta
  rows <- list()
  add <- function(term, type, est_t, se_t, bk) {
    rows[[length(rows) + 1L]] <<- data.frame(
      block = bk$term_label, term = term, type = type,
      est_t = est_t, se_t = se_t
    )
  }
  for (bk in fit$frame$re_blocks) {
    Vth <- Vfull[th_pos[bk$theta_idx], th_pos[bk$theta_idx],
                 drop = FALSE]
    t0 <- th[bk$theta_idx]
    if (bk$covstruct == "smooth") {
      add("sd(wiggle)", "sd", t0[1], sqrt(Vth[1, 1]), bk)
      next
    }
    if (bk$covstruct %in% c("gp", "hsgp")) {
      se_t <- sqrt(diag(Vth))
      add("sd(gp)", "sd", t0[1], se_t[1], bk)
      # hsgp estimates the lengthscale on brms's rescaled inputs, but the
      # reported range belongs in data units. The scale factor is a data
      # constant, so the shift on the log scale is exact and the se rides
      # through unchanged. The exact gp keeps the raw scale (dmax NULL).
      log_dmax <- log(bk$gp_dmax %||% 1)
      # iso: one shared range; otherwise one per dimension
      nr <- length(t0) - 1L
      for (j in seq_len(nr)) {
        term_j <- if (nr == 1L) "range(gp)" else {
          paste0("range(gp, ", bk$gp_vars[j], ")")
        }
        add(term_j, "range", t0[1 + j] + log_dmax, se_t[1 + j], bk)
      }
      next
    }
    if (bk$covstruct == "car") {
      se_t <- sqrt(diag(Vth))
      add("sd(car)", "sd", t0[1], se_t[1], bk)
      if (length(t0) > 1L) {
        # brms's names for the two mixing parameters, both on (0, 1)
        nm <- if (identical(bk$car_type, "bym2")) "rhocar" else "car"
        add(nm, "prop", t0[2], se_t[2], bk)
      }
      next
    }
    if (bk$covstruct == "spde") {
      # sigma and range are analytic functions of (log tau, log kappa):
      # log sigma = -log tau - log kappa - log(4 pi) / 2 and
      # log range = log(8) / 2 - log kappa, so the delta method is one
      # exact linear map with no differencing
      g_sd <- c(-1, -1)
      add("sd(spde)", "sd", log(spde_sd(t0)),
          sqrt(max(drop(g_sd %*% Vth %*% g_sd), 0)), bk)
      add("range(spde)", "range", log(spde_range(t0)),
          sqrt(Vth[2, 2]), bk)
      next
    }
    if (bk$covstruct == "equalto") next   # nothing estimated
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
    n_sd <- length(g0) - if (d > 1) d * (d - 1) / 2 else 0
    for (i in seq_len(n_sd)) {
      add(bk$cnms[min(i, length(bk$cnms))], "sd", g0[i], se_g[i], bk)
    }
    if (d > 1) {
      pairs <- which(lower.tri(diag(d)), arr.ind = TRUE)
      for (k in seq_len(nrow(pairs))) {
        add(paste0("cor(", bk$cnms[pairs[k, 2]], ",",
                   bk$cnms[pairs[k, 1]], ")"),
            "cor", g0[n_sd + k], se_g[n_sd + k], bk)
      }
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Back-transform to the natural scale (elementwise over types): log for
#' scales, Fisher-z for correlations, logit for the CAR mixing
#' proportions (which live on (0, 1), as they do in brms).
#'
#' @noRd
varcorr_untrans <- function(type, v) {
  ifelse(type == "cor", tanh(v),
         ifelse(type == "prop", 1 / (1 + exp(-v)), exp(v)))
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
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(200), g = factor(rep(1:20, 10)))
#' u <- cbind(rnorm(20, 0, 0.8), rnorm(20, 0, 0.4))
#' dd$y <- rnorm(200, 1 + 0.5 * dd$x + u[dd$g, 1] + u[dd$g, 2] * dd$x, 1)
#' fit <- frm(bf(y ~ x + (x | g)) + gaussian(), data = dd)
#'
#' # one row per SD and per correlation, on the scale they are read on
#' confint_varcorr(fit)
#'
#' # confint() reports the same parameters on their internal scale, so
#' # the bounds there are log-SDs and Fisher-z correlations
#' confint(fit)[grep("^theta", rownames(confint(fit))), ]
#'
#' # a fit with no random effects has no covariance parameters
#' confint_varcorr(frm(bf(y ~ x) + gaussian(), data = dd))
#' @export
confint_varcorr <- function(fit, level = 0.95) {
  tr <- varcorr_trans_rows(fit)
  if (is.null(tr)) return(NULL)
  z <- stats::qnorm(1 - (1 - level) / 2)
  data.frame(
    block = tr$block, term = tr$term, type = tr$type,
    estimate = varcorr_untrans(tr$type, tr$est_t),
    lwr = varcorr_untrans(tr$type, tr$est_t - z * tr$se_t),
    upr = varcorr_untrans(tr$type, tr$est_t + z * tr$se_t)
  )
}

#' Effective degrees of freedom of the smooth blocks: for an iid wiggly
#' block, `edf = k - tr(posterior cov)/prior variance` (the ridge
#' identity).
#'
#' @noRd
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

# Families whose linear predictor lives on a bounded probability scale,
# so an unbounded coefficient is evidence of separated data rather than
# of a large effect.
separation_families <- c("binomial", "bernoulli", "beta_binomial",
                         "zero_inflated_binomial",
                         "zero_inflated_beta_binomial")

#' Complete (or quasi-complete) separation: the maximum likelihood sits
#' at infinity, so the optimizer stops wherever its tolerances bite and
#' reports a huge coefficient with a standard error to match. lme4 and
#' glmmTMB both flag the pair rather than either half, because a
#' genuinely large effect on a well-populated cell keeps a small se.
#' `[glmmTMB diagnose()]`
#'
#' @noRd
diagnose_separation <- function(fit, ps) {
  rows <- list()
  for (lp in fit$frame$linpreds) {
    fam <- fit$spec$responses[[lp$resp]]$family
    if (!fam$family %in% separation_families) next
    if (!lp$dpar %in% (fam$primary_dpars %||% "mu")) next
    if (is.null(lp$X) || !ncol(lp$X)) next
    est <- ps$est[[lp$par]][lp$idx]
    se <- ps$se[[lp$par]][lp$idx]
    hit <- which(abs(est) > 10 & (!is.finite(se) | se > 10))
    for (i in hit) {
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = paste0(coef_block_key(fit, lp), ": ",
                           colnames(lp$X)[i]),
        estimate = est[i], std.error = se[i]
      )
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Predictor columns whose spread is many orders of magnitude away from
#' one: the objective's curvature then spans the same range, and the
#' optimizer's convergence tolerances are absolute. `autoscale = TRUE`
#' fixes it without touching the model. `[glmmTMB diagnose()]`
#'
#' @noRd
diagnose_predictor_scale <- function(fit, tol = 3) {
  rows <- list()
  for (lp in fit$frame$linpreds) {
    if (is.null(lp$X) || !ncol(lp$X)) next
    X <- as.matrix(lp$X)
    for (j in seq_len(ncol(X))) {
      if (identical(colnames(X)[j], "(Intercept)")) next
      s <- stats::sd(X[, j])
      if (!is.finite(s) || s <= 0) next
      if (abs(log10(s)) <= tol) next
      rows[[length(rows) + 1L]] <- data.frame(
        column = paste0(coef_block_key(fit, lp), ": ", colnames(X)[j]),
        sd = s
      )
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' lme4's isSingular: a variance component sitting on the boundary of its
#' parameter space - a standard deviation at zero, or a correlation at
#' +/-1. The verdict is read off the estimates alone, so it stands even
#' when the Hessian is positive definite and every gradient is tiny (a
#' collapsed component is a well-behaved optimum of a model the data
#' cannot support). `[lme4 test-isSingular.R, #660]`
#'
#' @noRd
diagnose_singular <- function(fit, tol = 1e-4) {
  if (!length(fit$frame$re_blocks)) return(NULL)
  vc <- tryCatch(as.data.frame(VarCorr(fit)), error = function(e) NULL)
  if (is.null(vc) || !nrow(vc)) return(NULL)
  is_cor <- !is.na(vc$var2)
  bad <- ifelse(is_cor, abs(vc$sdcor) > 1 - tol, vc$sdcor < tol)
  bad[is.na(bad)] <- FALSE
  if (!any(bad)) return(NULL)
  out <- data.frame(
    block = vc$grp[bad],
    term = ifelse(is_cor[bad],
                  paste0("cor(", vc$var1[bad], ",", vc$var2[bad], ")"),
                  paste0("sd(", vc$var1[bad], ")")),
    value = vc$sdcor[bad]
  )
  rownames(out) <- NULL
  out
}

#' The theta components that are LOG STANDARD DEVIATIONS, named the way
#' confint_varcorr() names their rows.
#'
#' The `|theta| > 8` near-singularity heuristic only reads as a boundary
#' fit on a log sd: `e^-8` is a variance no data supports. The other
#' components live on their own scales, where the same magnitude is
#' ordinary - a CAR/BYM2 mixing proportion and an AR(1) phi are logit-
#' and arctan-like, so a rho legitimately at the boundary sits at
#' `|theta| >> 8`, and the SPDE's (log tau, log kappa) are a precision
#' and an inverse range, neither of which is a standard deviation.
#' Reading those as singular fits is a false alarm on a converged model.
#' Structures whose registry declares no sd (rr, equalto, spde)
#' contribute nothing.
#'
#' @noRd
log_sd_theta_index <- function(fit) {
  idx <- integer(0)
  nms <- character(0)
  for (bk in fit$frame$re_blocks %||% list()) {
    # blocks whose confint() rows are hand-written carry their own name
    shared <- switch(bk$covstruct,
                     smooth = "sd(wiggle)", gp = "sd(gp)",
                     hsgp = "sd(gp)", car = "sd(car)",
                     spde = NA_character_, equalto = NA_character_,
                     NULL)
    if (!is.null(shared)) {
      if (is.na(shared)) next
      idx <- c(idx, bk$theta_idx[1L])
      nms <- c(nms, paste0(bk$term_label, " ", shared))
      next
    }
    reg <- covstruct_registry[[bk$covstruct]]
    si <- if (is.null(reg)) integer(0) else {
      tryCatch(as.integer(reg$sd_idx(bk$dim)),
               error = function(e) integer(0))
    }
    for (i in seq_along(si)) {
      idx <- c(idx, bk$theta_idx[si[i]])
      # the generic confint() path labels sd rows by column name, and
      # falls back to the first when one sd is shared across columns
      nms <- c(nms, paste0(bk$term_label, " ",
                           bk$cnms[min(i, length(bk$cnms))]))
    }
  }
  th_n <- length(fit$estimates$theta %||% numeric(0))
  keep <- !is.na(idx) & idx >= 1L & idx <= th_n
  stats::setNames(idx[keep], nms[keep])
}

#' Convergence diagnostics for a frmtmb fit
#'
#' Reports the optimizer's own verdict plus four checks that a converged
#' fit can still fail: non-finite standard errors, complete separation
#' in a binomial-type fit, predictor columns scaled far from one, and
#' variance components on the boundary of their parameter space
#' (lme4's `isSingular()`, read off the estimates rather than the
#' Hessian).
#'
#' @param fit A `frmtmb_fit`.
#' @param quiet If `TRUE`, return the diagnostics without printing.
#' @return Invisibly, a list of diagnostics.
#'
#' @srrstats {RE2.4b} Perfect collinearity between the predictors and the
#'   response is reported as complete separation: a binomial-type fit
#'   whose coefficients diverge because a predictor perfectly predicts
#'   the response is flagged by name, with the offending estimate. The
#'   test suite checks both that a separating design is flagged and that
#'   a well-behaved binomial fit is not.
#' @srrstats {RE4.7} Convergence statistics are available from the model
#'   object. `fit$opt$convergence` and `fit$opt$message` carry the
#'   optimizer's verdict, and `diagnose()` returns the maximum absolute
#'   gradient, the worst-offending parameter, the positive-definiteness
#'   of the Hessian, non-finite standard errors, the smallest eigenvalue
#'   of the covariance, boundary (singular) variance components,
#'   separation, and predictor scaling. `frm_allfit()` refits across
#'   optimizers as a further convergence check, and `check_laplace()`
#'   audits the approximation itself.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#' diagnose(fit)
#'
#' # a random effect the data cannot support collapses to the boundary,
#' # which is a valid fit but a warning about the model
#' dd$h <- factor(rep(1:10, each = 10))
#' fit_s <- frm(bf(y ~ x + (1 | g) + (1 | h)) + gaussian(), data = dd)
#' d <- diagnose(fit_s, quiet = TRUE)
#' d$singular
#'
#' # a predictor scaled far from one slows the optimizer down; the
#' # remedy is frmtmb_control(autoscale = TRUE)
#' dd$xbig <- dd$x * 1e5
#' diagnose(frm(bf(xbig ~ 1) + gaussian(), data = dd), quiet = TRUE)$scale
#' @export
diagnose <- function(fit, quiet = FALSE) {
  stopifnot(inherits(fit, "frmtmb_fit"))
  nm <- outer_par_names(fit)
  # a degenerate fit (no free outer parameters) has no gradient, no
  # covariance and no theta to report on
  degenerate <- !length(fit$opt$par)
  gr <- if (degenerate) numeric(0) else drop(fit$obj$gr(fit$opt$par))
  V <- sdr_of(fit)$cov.fixed
  se <- if (!length(V)) numeric(0) else sqrt(diag(V))
  ev <- if (!length(V)) NULL else {
    tryCatch(eigen(V, symmetric = TRUE, only.values = TRUE)$values,
             error = function(e) NULL)
  }
  # theta is absent from fits with no random effects; abs(NULL) is an
  # error, not an empty result
  th <- fit$estimates$theta %||% numeric(0)
  # only the log-sd components; see log_sd_theta_index()
  sd_i <- log_sd_theta_index(fit)
  ps <- tryCatch(par_est_se(fit), error = function(e) NULL)
  out <- list(
    convergence = fit$opt$convergence,
    message = fit$opt$message,
    max_grad = if (length(gr)) max(abs(gr)) else NA_real_,
    worst_grad = if (length(gr)) nm[which.max(abs(gr))] else NA_character_,
    pdHess = isTRUE(sdr_of(fit)$pdHess),
    bad_se = nm[!is.finite(se)],
    min_cov_eigenvalue = if (!is.null(ev) && length(ev)) min(ev),
    extreme_theta = sd_i[abs(th[sd_i]) > 8],
    separation = if (!is.null(ps)) diagnose_separation(fit, ps),
    predictor_scale = diagnose_predictor_scale(fit),
    singular = diagnose_singular(fit)
  )
  if (!quiet) {
    cat("Optimizer convergence code:", out$convergence,
        if (!is.null(out$message)) paste0("(", out$message, ")"), "\n")
    if (degenerate) {
      cat("No free parameters: the model is degenerate and the ",
          "likelihood was evaluated once\n", sep = "")
    } else {
      cat("Max |gradient|:", format(out$max_grad, digits = 4),
          "at", out$worst_grad, "\n")
    }
    cat("Hessian positive definite:", out$pdHess, "\n")
    if (length(out$bad_se)) {
      cat("Non-finite standard errors:",
          paste(out$bad_se, collapse = ", "), "\n")
    }
    if (length(out$extreme_theta)) {
      cat("Extreme covariance parameters (|log sd| > 8): ",
          paste(paste0(names(out$extreme_theta), " (log sd ",
                       format(th[out$extreme_theta], digits = 3), ")"),
                collapse = "; "),
          "\n  The fit is near-singular; consider simplifying the ",
          "random effects (e.g. diag() instead of a correlated term)\n",
          sep = "")
    }
    if (!is.null(out$singular)) {
      cat("Singular fit: ",
          paste(paste0(out$singular$block, " ", out$singular$term,
                       " = ", format(out$singular$value, digits = 3)),
                collapse = "; "),
          "\n  A variance component is on the boundary of its ",
          "parameter space. The fit is valid but the random-effect ",
          "structure is more complex than the data support; drop the ",
          "collapsed term or use diag() instead of a correlated ",
          "block.\n", sep = "")
    }
    if (!is.null(out$separation)) {
      cat("Likely complete separation: ",
          paste(paste0(out$separation$parameter, " = ",
                       format(out$separation$estimate, digits = 3),
                       " (se ", format(out$separation$std.error,
                                       digits = 3), ")"),
                collapse = "; "),
          "\n  Coefficients this large on the link scale with standard ",
          "errors to match mean the maximum likelihood is at infinity: ",
          "some combination of the predictors separates the outcome ",
          "perfectly. Drop or pool the offending predictor, or add a ",
          "prior (see set_prior()).\n", sep = "")
    }
    if (!is.null(out$predictor_scale)) {
      cat("Badly scaled predictors: ",
          paste(paste0(out$predictor_scale$column, " (sd ",
                       format(out$predictor_scale$sd, digits = 3), ")"),
                collapse = "; "),
          "\n  Rescale the column, or refit with ",
          "frmtmb_control(autoscale = TRUE).\n", sep = "")
    }
    clean <- out$convergence == 0 && out$pdHess && !length(out$bad_se) &&
      is.null(out$singular) && is.null(out$separation) &&
      is.null(out$predictor_scale) &&
      (degenerate || out$max_grad < 1e-3)
    if (clean) cat("No convergence problems detected\n")
  }
  invisible(out)
}

#' Largest absolute entry of a matrix, and 0 for an empty one.
#'
#' @noRd
maxabs <- function(M) if (!length(M)) 0 else max(abs(M))

#' Residual of B after projecting onto the column space of A.
#'
#' @noRd
proj_resid <- function(A, B) if (!ncol(A)) B else qr.resid(qr(A), B)

#' Test whether two design matrices span the same column space, up to a
#' relative tolerance. Used to decide when two REML likelihoods are
#' comparable.
#'
#' @noRd
same_column_space <- function(A, B, tol = 1e-8) {
  if (nrow(A) != nrow(B)) return(FALSE)
  if (!ncol(A) && !ncol(B)) return(TRUE)
  s <- max(1, maxabs(A), maxabs(B))
  if (!ncol(A) || !ncol(B)) return(FALSE)
  maxabs(proj_resid(A, B)) <= tol * s &&
    maxabs(proj_resid(B, A)) <= tol * s
}

#' Designs REML integrates out: the primary-dpar linear predictors, whose
#' coefficients live in the `beta` template component (dpar formulas keep
#' their coefficients in `betad` and stay outer).
#'
#' @noRd
reml_designs <- function(fit) {
  parts <- list()
  for (lp in fit$frame$linpreds) {
    if (!identical(lp$par, "beta")) next
    X <- if (is.null(lp$X)) {
      matrix(numeric(0), fit$frame$n_obs, 0L)
    } else {
      as.matrix(lp$X)
    }
    parts[[linpred_key(lp$resp, lp$dpar)]] <- X
  }
  parts[order(names(parts))]
}

#' A REML likelihood carries a `-1/2 log|X' V^-1 X|` term, so it is a
#' likelihood for a DIFFERENT quantity - the error contrasts - once X
#' changes, and differencing two of them is meaningless. It is perfectly
#' meaningful when the error contrasts are the same, which is exactly
#' when the fixed-effect designs span the same column space; that is the
#' usual REML comparison of variance-component structures. Refusing every
#' REML fit (the old behavior) refused that case too. `[glmmTMB#776]`
#'
#' @noRd
reml_comparable <- function(fits) {
  d1 <- reml_designs(fits[[1]])
  for (f in fits[-1]) {
    d2 <- reml_designs(f)
    if (!identical(names(d1), names(d2))) return(FALSE)
    for (k in names(d1)) {
      if (!same_column_space(d1[[k]], d2[[k]])) return(FALSE)
    }
  }
  TRUE
}

#' The cheapest correct REML -> ML conversion: reuse the assembled
#' design (no formula parsing, no frame assembly) and warm-start the
#' optimizer at the REML estimates. The parameter template is the same
#' list either way; REML only decides whether `beta` is integrated out,
#' so the REML estimates are a valid ML start.
#'
#' @noRd
anova_refit_ml <- function(fit) {
  ctl <- fit$control %||% frmtmb_control()
  # one anova() call can trigger several refits; a verbose original fit
  # must not make each of them narrate itself
  ctl$verbose <- FALSE
  fit_assembled(fit$spec, fit$frame, fit$bform, fit$call,
                REML = FALSE, start = NULL, control = ctl, se = FALSE,
                lower = fit$lower, upper = fit$upper,
                priors = fit$priors,
                quadrature = isTRUE(fit$quadrature),
                template = fit$estimates,
                data2 = fit$data2 %||% list())
}

#' Likelihood-ratio tests between nested frmtmb fits
#'
#' ML fits compare freely. REML fits compare only with each other, and
#' only when their fixed-effect designs span the same column space: a
#' REML likelihood is a likelihood for the error contrasts of that
#' design, so two of them are on a common scale exactly when the design
#' is the same. That covers the usual REML use - testing
#' variance-component structures with the fixed effects held fixed - and
#' refuses the rest with the reason (glmmTMB#776).
#'
#' `refit = TRUE` is the lme4 convenience for the refused case: every
#' REML fit in the comparison is refit with `REML = FALSE` and the ML
#' fits are compared instead. lme4 does this silently by default; here
#' it is opt-in and the message names the models that were refit.
#'
#' When the smaller model removes a variance component, the null value
#' sits on the boundary of the parameter space and the usual chi-square
#' reference is wrong: the asymptotic null is a mixture (for one
#' component, half a point mass at zero and half a chi-square with one
#' df), so the reported p-value is conservative - up to a factor of two
#' for a single component. lme4 and glmmTMB report the same naive
#' p-value; halve it for the one-component case, or use
#' [frm_bootstrap()] for a simulation-based reference.
#'
#' @param object A `frmtmb_fit`.
#' @param ... Further `frmtmb_fit` objects, nested with `object`.
#' @param refit If `TRUE`, refit every REML fit in the comparison with
#'   ML and compare those, with a message naming what was refit. The
#'   refits reuse the assembled design and warm-start at the REML
#'   estimates. `FALSE` (the default) keeps the REML fits and refuses
#'   the comparisons a restricted likelihood cannot support.
#' @return An `anova` table.
#'
#' @srrstats {RE4.11} Goodness-of-fit statistics are available for the
#'   fitted model. `logLik()` reports the log-likelihood with its degrees
#'   of freedom and `nobs`, so `AIC()` and `BIC()` work through the
#'   `stats` defaults; `extractAIC()` and `deviance()` are implemented;
#'   `anova()` gives likelihood-ratio tests between nested fits and
#'   `drop1()` single-term deletions. Effect sizes with the coefficients
#'   come from `summary()` (estimate, standard error, z, p) and
#'   `confint()`. The boundary problem for variance-component tests is
#'   documented above rather than left implicit.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' dd <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:20, 10)))
#' u <- cbind(rnorm(20, 0, 0.8), rnorm(20, 0, 0.5))
#' dd$y <- rnorm(n, 1 + 0.5 * dd$x + u[dd$g, 1] + u[dd$g, 2] * dd$x, 1)
#'
#' m0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#' m1 <- frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd)
#' anova(m0, m1)
#'
#' # dropping a variance component puts the null on the boundary, so
#' # this p-value is conservative by up to a factor of two
#' m2 <- frm(bf(y ~ x) + gaussian(), data = dd)
#' anova(m2, m0)
#'
#' # REML fits compare only when the fixed-effect designs agree, which
#' # is the case for a variance-component test
#' r0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd, REML = TRUE)
#' r1 <- frm(bf(y ~ x + (x | g)) + gaussian(), data = dd, REML = TRUE)
#' anova(r0, r1)
#' # differing designs are refused; refit = TRUE compares ML fits instead
#' rz <- frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd, REML = TRUE)
#' try(anova(r0, rz))
#' anova(r0, rz, refit = TRUE)
#' @export
anova.frmtmb_fit <- function(object, ..., refit = FALSE) {
  fits <- c(list(object), Filter(function(x) inherits(x, "frmtmb_fit"),
                                 list(...)))
  if (length(fits) < 2) {
    stop("anova() needs at least two frmtmb fits to compare", call. = FALSE)
  }
  reml <- vapply(fits, `[[`, TRUE, "REML")
  if (any(reml) && isTRUE(refit)) {
    labs <- vapply(fits[reml], model_label, "")
    message("anova(): refitting ", length(labs), " REML model",
            if (length(labs) != 1L) "s" else "", " with ML: ",
            paste(labs, collapse = "; "))
    fits[reml] <- lapply(fits[reml], anova_refit_ml)
    reml[] <- FALSE
  }
  if (any(reml)) {
    if (!all(reml)) {
      stop("anova() cannot mix REML and ML fits: their likelihoods are ",
           "for different quantities. Refit them all with the same ",
           "REML setting, or pass refit = TRUE to compare them as ML ",
           "fits", call. = FALSE)
    }
    if (!reml_comparable(fits)) {
      stop("REML likelihoods are comparable only between fits whose ",
           "fixed-effect designs span the same column space; these do ",
           "not. Pass refit = TRUE (or refit with REML = FALSE) to ",
           "compare fixed effects, or hold the fixed effects fixed to ",
           "compare random-effect structures", call. = FALSE)
    }
  }
  # Likelihoods computed on different data are not on a common scale, so
  # the LRT would be meaningless (and can come out negative). lme4 keys
  # its equivalent check off the `data` argument in the call, which both
  # false-positives on identical frames and misses NA-dropped rows;
  # comparing the response actually used catches the real cases.
  # [lme4#622]
  nobs_all <- vapply(fits, function(f) as.integer(f$frame$n_obs), 0L)
  if (length(unique(nobs_all)) > 1L) {
    stop("anova() needs fits with the same number of observations (got ",
         paste(unique(nobs_all), collapse = ", "),
         "); models fit to different data or with different NA rows ",
         "dropped are not comparable", call. = FALSE)
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
  rownames(tab) <- make.unique(vapply(fits, model_label, ""))
  structure(tab, class = c("anova", "data.frame"),
            heading = "Likelihood-ratio tests\n")
}

#' Single-term deletions
#'
#' Drops each fixed-effect term of the primary (`mu`) formula in turn,
#' refits, and tabulates AIC (and likelihood-ratio tests with
#' `test = "Chisq"`), following [stats::drop1()] and lme4's
#' `drop1.merMod`. Random-effect, smooth, and `mo()`/`mi()` terms are
#' not part of the deletion scope.
#'
#' @param object A `frmtmb_fit` from an ML fit (`REML = FALSE`) of a
#'   univariate model.
#' @param scope Terms to drop: a character vector or a right-hand-side
#'   formula. Defaults to all fixed-effect terms that marginality
#'   allows ([stats::drop.scope()]).
#' @param test `"Chisq"` adds likelihood-ratio tests.
#' @param k AIC penalty per parameter.
#' @param ... Unused.
#' @return An `anova` table with one row per dropped term.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), z = rnorm(100),
#'                  g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#' fit <- frm(bf(y ~ x * z + (1 | g)) + gaussian(), data = dd)
#'
#' # marginality keeps the main effects out of scope while x:z is in it
#' drop1(fit)
#' drop1(fit, test = "Chisq")
#'
#' # name the terms to override the default scope
#' drop1(fit, scope = ~ x + z, test = "Chisq")
#' @export
drop1.frmtmb_fit <- function(object, scope, test = c("none", "Chisq"),
                             k = 2, ...) {
  test <- match.arg(test)
  if (object$REML) {
    stop("drop1() compares fixed effects; refit with REML = FALSE",
         call. = FALSE)
  }
  if (length(object$spec$responses) > 1) {
    stop("drop1() is not supported for multivariate fits", call. = FALSE)
  }
  tt <- terms(object)
  labs <- attr(tt, "term.labels")
  if (missing(scope)) {
    scope <- stats::drop.scope(tt)
  } else if (!is.character(scope)) {
    scope <- attr(stats::terms(stats::update.formula(
      stats::formula(tt), scope)), "term.labels")
  }
  bad <- setdiff(scope, labs)
  if (length(bad)) {
    stop("scope is not a subset of the term labels: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }

  ll0 <- logLik(object)
  df0 <- attr(ll0, "df")
  aic0 <- -2 * as.numeric(ll0) + k * df0
  n_sc <- length(scope)
  ddf <- aic <- lrt <- rep(NA_real_, n_sc)
  for (i in seq_len(n_sc)) {
    # rebuild the bf() object with the term removed; the stored model
    # frame carries every variable, so the refit does not depend on
    # the original data still being visible
    nb <- object$bform
    nf <- stats::update.formula(nb$formula,
                                paste(". ~ . -", scope[i]))
    environment(nf) <- environment(nb$formula)
    nb$formula <- nf
    cl <- object$call
    cl$formula <- nb
    cl$data <- object$frame$data_frame
    # same reason for data2: the stored structural objects go in by
    # value, so the refit does not need the names the user passed to
    # still resolve where the call is evaluated
    if (length(object$data2)) cl$data2 <- object$data2
    fit_i <- eval(cl, environment(nb$formula) %||% parent.frame())
    ll_i <- logLik(fit_i)
    ddf[i] <- df0 - attr(ll_i, "df")
    aic[i] <- -2 * as.numeric(ll_i) + k * attr(ll_i, "df")
    lrt[i] <- 2 * (as.numeric(ll0) - as.numeric(ll_i))
  }
  tab <- data.frame(Df = c(NA, ddf), AIC = c(aic0, aic),
                    row.names = c("<none>", scope), check.names = FALSE)
  if (test == "Chisq") {
    tab$LRT <- c(NA, lrt)
    tab$`Pr(>Chi)` <- c(NA, stats::pchisq(lrt, ddf, lower.tail = FALSE))
  }
  structure(tab, class = c("anova", "data.frame"),
            heading = c("Single term deletions\n",
                        paste("Model:", model_label(object)), ""))
}

#' @export
update.frmtmb_fit <- function(object, ..., evaluate = TRUE) {
  cl <- object$call
  extras <- match.call(expand.dots = FALSE)$...
  for (nm in names(extras)) cl[[nm]] <- extras[[nm]]
  if (!evaluate) return(cl)
  # the stored structural objects go into the call by value, so an
  # update in a session where the original data2 names are gone still
  # assembles; an explicit data2 = in the update wins
  if (length(object$data2) && !("data2" %in% names(extras))) {
    cl$data2 <- object$data2
  }
  eval(cl, parent.frame())
}

#' Likelihood profiles
#'
#' Wraps [TMB::tmbprofile()] per parameter. The returned objects have
#' `plot()` and `confint()` methods (from TMB).
#'
#' @param fitted A `frmtmb_fit`.
#' @param parm Parameter names (as in `confint()` rownames) or indices.
#'   Required; profiling is not free, so there is no all-parameters
#'   default.
#' @param ... Passed to [TMB::tmbprofile()].
#' @return A `tmbprofile` data frame, or a named list of them when
#'   `parm` has length above one.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#' # parameter names are the confint() row names
#' rownames(confint(fit))
#' pr <- profile(fit, "theta_1")
#' plot(pr)
#' # TMB's confint() reads the interval off the profile
#' confint(pr)
#'
#' # several parameters at once return a named list
#' prs <- profile(fit, c("x", "theta_1"))
#' names(prs)
#' @export
profile.frmtmb_fit <- function(fitted, parm, ...) {
  nm <- outer_par_names(fitted)
  idx <- if (is.numeric(parm)) as.integer(parm) else match(parm, nm)
  if (anyNA(idx)) {
    stop("Unknown parameter(s) in profile(parm =): ",
         paste(parm[is.na(match(parm, nm))], collapse = ", "),
         ". Available: ", paste(nm, collapse = ", "), call. = FALSE)
  }
  out <- lapply(idx, function(i) {
    TMB::tmbprofile(fitted$obj, name = i, trace = FALSE, ...)
  })
  names(out) <- nm[idx]
  if (length(out) == 1L) out[[1L]] else out
}

## hypothesis(): the expression environment and its parameter mapping.

#' Make a term or group name usable as a variable name in a hypothesis
#' expression: drop parentheses and every other character an R name
#' cannot carry.
#'
#' @noRd
hyp_san <- function(s) gsub("[^[:alnum:]_.]", "", gsub("[()]", "", s))

#' Parameter values without any covariance machinery (usable on refits
#' inside a bootstrap without triggering sdreport).
#'
#' @noRd
hyp_vals_only <- function(fit) {
  est <- fit$estimates
  bd <- est$betad
  if (length(fx <- fit$frame$betad_fixed_idx)) bd <- bd[-fx]
  list(
    vals = c(est$beta, bd, est$theta, est$thetar),
    comp = c(rep("beta", length(est$beta)), rep("betad", length(bd)),
             rep("theta", length(est$theta)),
             rep("thetar", length(est$thetar)))
  )
}

#' Values plus joint covariance of (beta, estimated betad, theta,
#' thetar). ML: straight from cov.fixed in opt$par order (outer_pos maps
#' back into the full outer vector for tmbroot lincombs). REML: beta is
#' integrated out, so the blocks come from the joint precision.
#'
#' @noRd
hyp_par_cov <- function(fit) {
  comps <- c("beta", "betad", "theta", "thetar")
  if (!fit$REML && !isTRUE(fit$control$profile)) {
    sdr <- sdr_of(fit)
    V <- sdr$cov.fixed
    rn <- rownames(V)
    keep <- which(rn %in% comps)
    # par.fixed equals opt$par at the optimum (this branch is never
    # taken under control profile = TRUE)
    list(vals = unname(sdr$par.fixed[keep]), comp = rn[keep],
         V = V[keep, keep, drop = FALSE], outer_pos = keep,
         n_outer = length(fit$opt$par))
  } else {
    Q <- sdr_of(fit)$jointPrecision
    Vall <- solve_joint_precision(Q, fit$cache)
    rn <- rownames(Q)
    keep <- which(rn %in% comps)
    vo <- hyp_vals_only(fit)
    vals <- numeric(length(keep))
    cnt <- stats::setNames(integer(length(comps)), comps)
    for (i in seq_along(keep)) {
      k <- rn[keep[i]]
      cnt[k] <- cnt[k] + 1L
      vals[i] <- vo$vals[vo$comp == k][cnt[k]]
    }
    list(vals = vals, comp = rn[keep],
         V = as.matrix(Vall[keep, keep, drop = FALSE]), outer_pos = NULL,
         n_outer = length(fit$opt$par))
  }
}

#' Named list the hypothesis expressions are evaluated in: fixed
#' coefficients under their vcov() names (parentheses stripped),
#' natural-scale random-effect summaries (`sd_<group>__<term>`,
#' `cor_<group>__<t1>__<t2>`), and `sigma` when it is a scalar.
#'
#' @noRd
hyp_env_vals <- function(fit, vals, comp) {
  env <- list()
  cf <- c(vals[comp == "beta"], vals[comp == "betad"])
  cn <- gsub("[()]", "", estimated_coef_names(fit))
  for (i in seq_along(cn)) env[[cn[i]]] <- cf[i]

  th <- vals[comp == "theta"]
  for (bk in fit$frame$re_blocks) {
    if (bk$covstruct %in% c("smooth", "gr_cov", "gr_prec",
                            "gp", "hsgp", "equalto", "car", "spde")) next
    V <- covstruct_registry[[bk$covstruct]]$vcov(th[bk$theta_idx], bk)
    tn <- hyp_san(bk$cnms)
    g <- hyp_san(bk$group_name)
    sds <- sqrt(diag(V))
    for (j in seq_along(sds)) {
      nm <- paste0("sd_", g, "__", tn[j])
      if (is.null(env[[nm]])) env[[nm]] <- sds[j]
    }
    if (nrow(V) > 1L) {
      C <- stats::cov2cor(V)
      for (j in seq_len(nrow(V) - 1L)) {
        for (k in seq(j + 1L, nrow(V))) {
          nm <- paste0("cor_", g, "__", tn[j], "__", tn[k])
          if (is.null(env[[nm]])) env[[nm]] <- C[j, k]
        }
      }
    }
  }

  if (length(fit$spec$responses) == 1L && is.null(env[["sigma"]])) {
    # the guard keeps a covariate literally named `sigma` visible
    for (lp in fit$frame$linpreds) {
      if (lp$dpar != "sigma") next
      if (!is.null(lp$constant)) {
        env[["sigma"]] <- lp$constant
      } else if (ncol(lp$X) == 1L &&
                 identical(colnames(lp$X), "(Intercept)") &&
                 is.null(lp$Z) && lp$par == "betad") {
        tpl_len <- length(fit$frame$par_template$betad)
        rk <- match(lp$idx, setdiff(seq_len(tpl_len),
                                    fit$frame$betad_fixed_idx))
        bd <- vals[comp == "betad"]
        if (!is.na(rk)) env[["sigma"]] <- lp$link$linkinv(bd[rk])
      }
    }
  }
  env
}

#' Turn one hypothesis string into a language object to evaluate. An
#' `"lhs = rhs"` hypothesis becomes the difference of the two sides, so
#' every hypothesis is then tested against zero.
#'
#' @noRd
hyp_parse <- function(h) {
  eq <- strsplit(h, "=", fixed = TRUE)[[1L]]
  txt <- if (length(eq) == 2L) {
    paste0("(", eq[1L], ") - (", eq[2L], ")")
  } else if (length(eq) == 1L) {
    h
  } else {
    stop("A hypothesis has at most one '=': '", h, "'", call. = FALSE)
  }
  str2lang(txt)
}

#' Evaluate a parsed hypothesis at one parameter vector. A failure lists
#' the available names, because an unknown name is the usual cause.
#'
#' @noRd
hyp_eval <- function(fit, ex, vals, comp) {
  ev <- hyp_env_vals(fit, vals, comp)
  tryCatch(eval(ex, ev), error = function(e) {
    stop(conditionMessage(e), "\nAvailable names: ",
         paste(names(ev), collapse = ", "), call. = FALSE)
  })
}

#' Central-difference gradient of a scalar function of the parameter
#' vector. The delta method needs a gradient, and a hypothesis is an
#' arbitrary R expression with no derivative available.
#'
#' @noRd
hyp_fd_grad <- function(f, v) {
  vapply(seq_along(v), function(i) {
    step <- max(1e-5, 1e-5 * abs(v[i]))
    vp <- v; vp[i] <- vp[i] + step
    vm <- v; vm[i] <- vm[i] - step
    (f(vp) - f(vm)) / (2 * step)
  }, numeric(1))
}

#' Hypothesis tests on parameter expressions
#'
#' The frequentist analog of brms's `hypothesis()`: evaluates
#' expressions of the model parameters at the estimates and tests them
#' against zero. A hypothesis is `"expr"` (tested against 0) or
#' `"expr = rhs"`, e.g. `"x1 - x2 = 0"` or `"exp(Intercept) = 1"`.
#'
#' Available names: the fixed-effect coefficients under their `vcov()`
#' row names with parentheses stripped (`Intercept`, `x`,
#' `sigma_Intercept`, ...), natural-scale random-effect summaries
#' `sd_<group>__<term>` and `cor_<group>__<t1>__<t2>` (brms naming),
#' and `sigma` when the residual SD is a scalar. So an ICC is
#' `"sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)"`.
#' [variables()] lists every usable name for a fit.
#'
#' Methods:
#' - `"wald"` (default): delta-method z-test, finite-difference
#'   gradient against the joint parameter covariance (under REML, from
#'   the joint precision).
#' - `"profile"`: profile-likelihood interval via [TMB::tmbroot()] with
#'   a `lincomb` direction. Only for hypotheses that are linear in the
#'   parameters, and only for ML fits; `se`, `z`, and `p` stay
#'   Wald-based - the method changes the interval.
#' - `"boot"`: parametric bootstrap through [frm_bootstrap()]
#'   (percentile interval; `p` is the two-sided percentile p-value,
#'   whose resolution is limited by `nsim`; `se` is the bootstrap SD).
#'   Handles any expression, including the variance-component names,
#'   whose sampling distributions Wald approximates poorly.
#'
#' For a [frm_multiple()] result the Wald estimate and delta-method
#' variance are computed per imputation and pooled by Rubin's rules
#' with Barnard-Rubin degrees of freedom; the returned table carries
#' `t` and `df` columns in place of `z` (reference t distribution, not
#' normal), and only Wald inference is available.
#'
#' @param x A `frmtmb_fit`, or a `frmtmb_multiple` for pooled tests.
#' @param hypothesis Character vector of hypotheses.
#' @param alpha Test level; the reported interval covers `1 - alpha`
#'   (brms spelling).
#' @param method `"wald"`, `"profile"`, or `"boot"`.
#' @param nsim Bootstrap draws for `method = "boot"`; all hypotheses
#'   share one bootstrap run.
#' @param seed Optional seed for `method = "boot"`.
#' @param ... Backend controls: passed to [TMB::tmbprofile()] for
#'   `method = "profile"` (e.g. `ytol`, `ystep`, `maxit`,
#'   `parm.range`) and to [frm_bootstrap()] for `method = "boot"`
#'   (e.g. `re.form = NULL` for a conditional bootstrap). Unused for
#'   `"wald"` (a warning).
#' @return A `frmtmb_hypothesis` object: a data frame with one row per
#'   hypothesis (`estimate`, `se`, `lwr`, `upr`, `z`, `p`) carrying the
#'   method payload in attributes - the bootstrap draws matrix
#'   (`attr(., "draws")`) or the profile curves (`attr(., "profiles")`).
#'   `plot()` shows the bootstrap distribution, the profile curve, or
#'   the implied Wald normal density, one panel per hypothesis.
#'   Subsetting with `[` drops the attributes; keep the full object for
#'   plotting.
#' @examples
#' set.seed(4)
#' dd <- data.frame(x1 = rnorm(120), x2 = rnorm(120),
#'                  g = factor(rep(1:10, 12)))
#' dd$y <- rnorm(120, 1 + 0.6 * dd$x1 + 0.4 * dd$x2 +
#'                 rnorm(10, 0, 0.5)[dd$g], 1)
#' fit <- frm(bf(y ~ x1 + x2 + (1 | g)) + gaussian(), data = dd)
#' hypothesis(fit, c("x1 - x2 = 0", "exp(Intercept)"))
#' # variance-component expressions: an ICC with bootstrap intervals
#' hypothesis(fit, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)",
#'            method = "boot", nsim = 20, seed = 1)
#' @export
hypothesis <- function(x, ...) UseMethod("hypothesis")

#' Usable parameter names
#'
#' The names that [hypothesis()] expressions (and `set_prior()`
#' targeting) accept: fixed-effect coefficients under their `vcov()`
#' names with parentheses stripped, natural-scale random-effect
#' summaries (`sd_<group>__<term>`, `cor_<group>__<t1>__<t2>`), and
#' `sigma` when the residual SD is a scalar. The brms spelling; for
#' sampled fits, `variables()` on the [frm_sample()] result lists the
#' draw columns instead.
#'
#' @param x A `frmtmb_fit` or `frmtmb_draws`.
#' @param ... Unused.
#' @return A character vector.
#' @examples
#' dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
#' dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#' variables(fit)
#' @export
variables <- function(x, ...) UseMethod("variables")

#' @rdname variables
#' @export
variables.frmtmb_fit <- function(x, ...) {
  vo <- hyp_vals_only(x)
  names(hyp_env_vals(x, vo$vals, vo$comp))
}

#' @rdname hypothesis
#' @export
hypothesis.frmtmb_fit <- function(x, hypothesis, alpha = 0.05,
                                  method = c("wald", "profile", "boot"),
                                  nsim = 500, seed = NULL, ...) {
  method <- match.arg(method)
  if (method == "wald" && ...length()) {
    warning("ignoring arguments unused by method = 'wald': ",
            paste(...names(), collapse = ", "), call. = FALSE)
  }
  z_q <- stats::qnorm(1 - alpha / 2)
  exs <- lapply(hypothesis, hyp_parse)
  vo <- hyp_vals_only(x)
  vals0 <- vapply(seq_along(exs), function(i) {
    val <- hyp_eval(x, exs[[i]], vo$vals, vo$comp)
    if (!is.numeric(val) || length(val) != 1L) {
      stop("Hypothesis '", hypothesis[i], "' must evaluate to a single ",
           "number at the fitted estimates", call. = FALSE)
    }
    val
  }, numeric(1))

  hyp_result <- function(out, extra = list()) {
    rownames(out) <- NULL
    attr(out, "method") <- method
    attr(out, "alpha") <- alpha
    for (nm in names(extra)) attr(out, nm) <- extra[[nm]]
    class(out) <- c("frmtmb_hypothesis", "data.frame")
    out
  }

  if (method == "boot") {
    FUN <- function(ft) {
      w <- hyp_vals_only(ft)
      vapply(exs, function(ex) hyp_eval(ft, ex, w$vals, w$comp),
             numeric(1))
    }
    bs <- frm_bootstrap(x, FUN, nsim = nsim, seed = seed, ...)
    colnames(bs$t) <- hypothesis
    rows <- lapply(seq_along(exs), function(i) {
      t_i <- bs$t[, i]
      t_i <- t_i[is.finite(t_i)]
      nb <- length(t_i)
      se <- stats::sd(t_i)
      p <- min(1, 2 * min((1 + sum(t_i <= 0)) / (1 + nb),
                          (1 + sum(t_i >= 0)) / (1 + nb)))
      data.frame(hypothesis = hypothesis[i], estimate = vals0[i],
                 se = se,
                 lwr = unname(stats::quantile(t_i, alpha / 2)),
                 upr = unname(stats::quantile(t_i, 1 - alpha / 2)),
                 z = vals0[i] / se, p = p)
    })
    return(hyp_result(do.call(rbind, rows),
                      list(draws = bs$t, nsim = nsim,
                           converged = bs$converged)))
  }

  pc <- hyp_par_cov(x)
  profiles <- vector("list", length(exs))
  rows <- vector("list", length(exs))
  for (i in seq_along(exs)) {
    ex <- exs[[i]]
    fn <- function(v) hyp_eval(x, ex, v, pc$comp)
    g <- hyp_fd_grad(fn, pc$vals)
    se <- sqrt(max(0, drop(t(g) %*% pc$V %*% g)))
    zv <- vals0[i] / se
    lwr <- vals0[i] - z_q * se
    upr <- vals0[i] + z_q * se
    if (method == "profile") {
      if (x$REML) {
        stop("method = 'profile' requires an ML fit (REML integrates ",
             "the fixed effects out of the outer problem)",
             call. = FALSE)
      }
      if (isTRUE(x$control$profile)) {
        stop("hypothesis(method = 'profile') needs a fit without ",
             "frmtmb_control(profile = TRUE)", call. = FALSE)
      }
      g2 <- hyp_fd_grad(fn, pc$vals + 0.1 * (1 + abs(pc$vals)))
      if (max(abs(g - g2)) > 1e-4 * max(1, max(abs(g)))) {
        stop("Hypothesis '", hypothesis[i], "' is not linear in the ",
             "parameters; use method = 'boot'", call. = FALSE)
      }
      v <- numeric(pc$n_outer)
      v[pc$outer_pos] <- g
      const <- vals0[i] - sum(g * pc$vals)
      pargs <- utils::modifyList(
        list(obj = x$obj, lincomb = v, trace = FALSE,
             ytol = 0.5 * stats::qchisq(1 - alpha, 1) + 1),
        list(...)
      )
      pr <- do.call(TMB::tmbprofile, pargs)
      ci <- stats::confint(pr, level = 1 - alpha)
      pr[[1L]] <- pr[[1L]] + const
      profiles[[i]] <- pr
      lwr <- unname(ci[1]) + const
      upr <- unname(ci[2]) + const
    }
    rows[[i]] <- data.frame(hypothesis = hypothesis[i],
                            estimate = vals0[i], se = se,
                            lwr = lwr, upr = upr, z = zv,
                            p = 2 * stats::pnorm(-abs(zv)))
  }
  hyp_result(do.call(rbind, rows),
             if (method == "profile") {
               list(profiles = stats::setNames(profiles, hypothesis))
             } else {
               list()
             })
}

#' @export
print.frmtmb_hypothesis <- function(x, digits = 4, ...) {
  method <- attr(x, "method")
  cat("Hypothesis tests (method = ", method, ")\n", sep = "")
  if (identical(method, "boot")) {
    cat("  bootstrap draws: ", attr(x, "nsim"), " (",
        sum(!attr(x, "converged")), " failed or not converged)\n",
        sep = "")
  }
  df <- x
  class(df) <- "data.frame"
  df[-1] <- lapply(df[-1], signif, digits)
  print(df, row.names = FALSE)
  invisible(x)
}

#' @export
plot.frmtmb_hypothesis <- function(x, ask = NULL, ...) {
  method <- attr(x, "method") %||% "wald"
  alpha <- attr(x, "alpha") %||% 0.05
  n <- nrow(x)
  ask <- ask %||% (n > 1L && grDevices::dev.interactive())
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }
  mark <- function(i) {
    graphics::abline(v = x$estimate[i], lwd = 2)
    graphics::abline(v = c(x$lwr[i], x$upr[i]), lty = 2)
    graphics::abline(v = 0, col = 2)
  }
  for (i in seq_len(n)) {
    h <- x$hypothesis[i]
    if (method %in% c("boot", "posterior")) {
      d <- attr(x, "draws")[, i]
      d <- d[is.finite(d)]
      graphics::hist(d, freq = FALSE, breaks = "FD", main = h,
                     xlab = if (method == "boot") "bootstrap value" else
                       "posterior value",
                     col = "gray90", border = "gray60")
      if (length(unique(d)) > 1L) {
        graphics::lines(stats::density(d), lwd = 2)
      }
      mark(i)
    } else if (method == "profile") {
      pr <- attr(x, "profiles")[[i]]
      dnll <- pr$value - min(pr$value, na.rm = TRUE)
      graphics::plot(pr[[1L]], dnll, type = "l", lwd = 2, main = h,
                     xlab = "value",
                     ylab = "profile neg. log-likelihood change")
      graphics::abline(h = 0.5 * stats::qchisq(1 - alpha, 1), lty = 3)
      mark(i)
    } else {
      xs <- seq(x$estimate[i] - 4 * x$se[i], x$estimate[i] + 4 * x$se[i],
                length.out = 200)
      graphics::plot(xs, stats::dnorm(xs, x$estimate[i], x$se[i]),
                     type = "l", lwd = 2, main = h, xlab = "value",
                     ylab = "Wald (normal) density")
      mark(i)
    }
  }
  invisible(x)
}
