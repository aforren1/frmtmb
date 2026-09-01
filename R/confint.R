# Confidence intervals, convergence diagnostics, model comparison.

# One-line label for a fit including dpar formulas, so two models that
# share a primary formula (e.g. plain vs distributional) stay
# distinguishable in anova tables.
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

# Names of the outer (optimized) parameters, in obj$par order: template
# component order, minus `random` components, minus mapped entries.
outer_par_names <- function(fit) {
  tpl <- fit$frame$par_template
  random <- if (fit$REML) c("b", "beta") else "b"
  # control profile = TRUE moves beta into the inner problem
  if (isTRUE(fit$control$profile)) random <- union(random, "beta")
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
    sdr <- sdr_of(object)
    se <- sqrt(diag(sdr$cov.fixed))
    # under control profile = TRUE the profiled betas are absent from
    # opt$par but present in par.fixed
    if (length(est) != length(se)) est <- sdr$par.fixed
    ci <- cbind(lwr = est + stats::qnorm(a) * se,
                upr = est + stats::qnorm(1 - a) * se,
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
    if (bk$covstruct %in% c("gp", "hsgp")) {
      se_t <- sqrt(diag(Vth))
      rows[[length(rows) + 1L]] <- data.frame(
        block = bk$term_label, term = "sd(gp)", type = "sd",
        estimate = exp(t0[1]),
        lwr = exp(t0[1] - z * se_t[1]), upr = exp(t0[1] + z * se_t[1])
      )
      # iso: one shared range; otherwise one per dimension
      nr <- length(t0) - 1L
      for (j in seq_len(nr)) {
        term_j <- if (nr == 1L) "range(gp)" else {
          paste0("range(gp, ", bk$gp_vars[j], ")")
        }
        rows[[length(rows) + 1L]] <- data.frame(
          block = bk$term_label, term = term_j, type = "range",
          estimate = exp(t0[1 + j]),
          lwr = exp(t0[1 + j] - z * se_t[1 + j]),
          upr = exp(t0[1 + j] + z * se_t[1 + j])
        )
      }
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
  rownames(tab) <- make.unique(vapply(fits, model_label, ""))
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
#' @export
profile.frmtmb_fit <- function(fitted, parm, ...) {
  nm <- outer_par_names(fitted)
  idx <- if (is.numeric(parm)) as.integer(parm) else match(parm, nm)
  if (anyNA(idx)) {
    stop("Unknown parameter(s): ",
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

hyp_san <- function(s) gsub("[^[:alnum:]_.]", "", gsub("[()]", "", s))

# Parameter values without any covariance machinery (usable on refits
# inside a bootstrap without triggering sdreport).
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

# Values plus joint covariance of (beta, estimated betad, theta,
# thetar). ML: straight from cov.fixed in opt$par order (outer_pos maps
# back into the full outer vector for tmbroot lincombs). REML: beta is
# integrated out, so the blocks come from the joint precision.
hyp_par_cov <- function(fit) {
  comps <- c("beta", "betad", "theta", "thetar")
  if (!fit$REML && !isTRUE(fit$control$profile)) {
    sdr <- sdr_of(fit)
    V <- sdr$cov.fixed
    rn <- rownames(V)
    keep <- which(rn %in% comps)
    # par.fixed equals opt$par at the optimum and, unlike opt$par,
    # includes betas profiled by control profile = TRUE
    list(vals = unname(sdr$par.fixed[keep]), comp = rn[keep],
         V = V[keep, keep, drop = FALSE], outer_pos = keep,
         n_outer = length(fit$opt$par))
  } else {
    Q <- sdr_of(fit)$jointPrecision
    Vall <- solve(Q)
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

# Named list the hypothesis expressions are evaluated in: fixed
# coefficients under their vcov() names (parentheses stripped),
# natural-scale random-effect summaries (sd_<group>__<term>,
# cor_<group>__<t1>__<t2>), and `sigma` when it is a scalar.
hyp_env_vals <- function(fit, vals, comp) {
  env <- list()
  cf <- c(vals[comp == "beta"], vals[comp == "betad"])
  cn <- gsub("[()]", "", estimated_coef_names(fit))
  for (i in seq_along(cn)) env[[cn[i]]] <- cf[i]

  th <- vals[comp == "theta"]
  for (bk in fit$frame$re_blocks) {
    if (bk$covstruct %in% c("smooth", "gr_cov", "gr_prec",
                            "gp", "hsgp", "equalto")) next
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

  if (length(fit$spec$responses) == 1L) {
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

hyp_eval <- function(fit, ex, vals, comp) {
  ev <- hyp_env_vals(fit, vals, comp)
  tryCatch(eval(ex, ev), error = function(e) {
    stop(conditionMessage(e), "\nAvailable names: ",
         paste(names(ev), collapse = ", "), call. = FALSE)
  })
}

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
#' @param x A `frmtmb_fit`.
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
           "number", call. = FALSE)
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
