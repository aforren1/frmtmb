# Prediction, fitted values, residuals, and simulation.
#
# Every linear predictor's Z spans the full b vector, so Z column indices
# equal b indices; |ID|-merged blocks need no special casing on the
# in-sample path. For newdata, each block contribution is rebuilt from the
# block's components (one component per contributing linear predictor).

# Memoized joint covariance of all estimated parameters (fixed + random),
# with row/col component names. Cached in fit$cache.
get_joint_cov <- function(fit) {
  cache <- fit$cache
  if (!is.null(cache$Vjoint)) return(cache$Vjoint)
  Q <- sdr_of(fit)$jointPrecision
  if (is.null(Q) && !is.null(sdr_of(fit)$par.random) &&
      length(sdr_of(fit)$par.random)) {
    Q <- autoscale_sdreport(fit, jp = TRUE)$jointPrecision
  }
  if (is.null(Q)) {
    V <- sdr_of(fit)$cov.fixed
    rn <- rownames(V)
  } else {
    V <- as.matrix(Matrix::solve(Q))
    rn <- rownames(Q)
  }
  cache$Vjoint <- list(V = V, names = rn)
  cache$Vjoint
}

# Numeric coefficient-space vector for a fitted model (rr factors
# expanded through the loadings; identity otherwise).
coef_b <- function(fit, b = fit$estimates[["b"]]) {
  if (is.null(b)) return(b)
  expand_b(fit$frame, b, fit$estimates$theta)
}

# Numeric mo() column values: D times the cumulative simplex at the
# category codes, evaluated at the current simplex estimates.
mo_col_values <- function(fit, mi, codes = mi$codes) {
  zeta <- exp(c(0, fit$estimates[[mi$zeta]]))
  zeta <- zeta / sum(zeta)
  cz0 <- c(0, cumsum(zeta))
  mi$D * cz0[codes + 1L]
}

# Category codes of a mo() variable in new data, validated against the
# fitted range.
mo_codes <- function(fit, lp, mi, newdata) {
  env <- fit$spec$responses[[lp$resp]]$formula_env
  v <- eval(mi$expr, newdata, env)
  if (!is.null(mi$levels)) {
    v <- factor(v, levels = mi$levels, ordered = TRUE)
    if (anyNA(v)) {
      stop("mo(): new data contain unknown categories", call. = FALSE)
    }
    return(as.integer(v) - 1L)
  }
  codes <- as.integer(round(v))   # grids may land between categories
  if (any(codes < 0 | codes > mi$D)) {
    stop("mo(): new data outside the fitted range 0..", mi$D,
         call. = FALSE)
  }
  codes
}

# Observed-or-latent values of a mi() response at the estimates.
mi_values <- function(fit, vn) {
  xv <- fit$frame$y[[vn]]
  mm_ <- fit$frame$mi_map[[vn]]
  if (!is.null(mm_)) xv[mm_$rows] <- fit$estimates$miss[mm_$idx]
  xv
}

# Fill the zero placeholder columns of a stored design matrix with the
# mo() and mi() values at the current estimates.
patch_mo_cols <- function(fit, lp, X) {
  for (mi in lp$mo %||% list()) {
    v <- mo_col_values(fit, mi)
    if (!is.null(mi$mult)) v <- v * mi$mult
    X[, mi$col] <- v
  }
  for (mt in lp$mi %||% list()) {
    v <- mi_values(fit, mt$var)
    if (!is.null(mt$mult)) v <- v * mt$mult
    X[, mt$col] <- v
  }
  X
}

# xlevels restricted to the variables a terms object actually uses;
# extra entries make model.frame warn.
xlev_for <- function(xlevels, tt) {
  xlevels[intersect(names(xlevels), all.vars(tt))]
}

# Rebuild the design pieces of one linear predictor for new data:
# dense X (parametric + smooth null-space columns), per-block RE
# component designs with level indices, and the offset.
pred_design <- function(fit, lp, newdata, allow_new_levels = FALSE,
                        use_re = TRUE) {
  env <- fit$spec$responses[[lp$resp]]$formula_env
  tt <- patch_predvars(lp$terms, fit$frame$predvar_map)
  mfp <- stats::model.frame(tt, newdata, na.action = stats::na.pass,
                            xlev = xlev_for(lp$xlevels, tt))
  # sparse_x fits keep newdata designs sparse too; NA rows must stay NA
  # in eta, and sparse.model.matrix silently zeroes NA factor rows, so
  # frames with NAs fall back to the dense builder
  X <- if (isTRUE(fit$frame$sparse_x) && !anyNA(mfp)) {
    sparse_mm(tt, mfp, contrasts.arg = lp$contrasts)
  } else {
    stats::model.matrix(tt, mfp, contrasts.arg = lp$contrasts)
  }
  # Estimability against a rank-deficient fit. A prediction is a linear
  # functional of beta, so it is identified only when the new design row
  # is orthogonal to every direction the fit could not resolve. Testing
  # the frozen null space (not merely "a dropped column is nonzero")
  # keeps rows that restate a kept column - x2 = 2 * x, or a cell whose
  # aliased indicator is implied by the kept ones - exact.
  nonest <- rep(FALSE, nrow(X))
  if (!is.null(lp$alias_null)) {
    Xa <- as.matrix(X[, rownames(lp$alias_null), drop = FALSE])
    scl <- max(abs(Xa[is.finite(Xa)]), 1)
    # NA rows already predict NA; keep them out of the estimability vote
    nonest <- rowSums(abs(Xa %*% lp$alias_null) > 1e-8 * scl,
                      na.rm = TRUE) > 0
  }
  # frozen intercept-drop / rank-deficiency column set from fit time
  X <- X[, lp$param_colnames, drop = FALSE]
  off <- extract_offset(tt, mfp, env)

  # Smooths: rebuild the (wiggly, fixed) split of the basis for newdata.
  # PredictMat %*% U, scaled by D, has the wiggly columns first (in rand
  # order) and the null-space columns last.
  sm_parts <- list()
  for (si in lp$smooths %||% list()) {
    M <- mgcv::PredictMat(si$sm, newdata) %*% si$U
    M <- sweep(M, 2, si$D, `*`)
    pos <- 0L
    for (r in seq_along(si$nr)) {
      Xr_new <- M[, pos + seq_len(si$nr[r]), drop = FALSE]
      pos <- pos + si$nr[r]
      sm_parts[[length(sm_parts) + 1L]] <- list(
        bk = fit$frame$re_blocks[[si$block_ids[r]]],
        Xr = Xr_new
      )
    }
    if (si$nf > 0L) {
      X <- cbind(X, M[, pos + seq_len(si$nf), drop = FALSE])
    }
  }

  # mo()/mi() columns come after the smooth null-space columns, matching
  # the fitted X layout; mo values use the current simplex, mi values
  # must be supplied complete in newdata
  # gp() contributions rebuild their basis at the new positions and
  # ride the smooth-parts machinery (curve kept at population level).
  # Exact gp at unseen positions kriges: conditional-mean weights
  # K* K^-1 at the fitted kernel slot into Xr, and the conditional
  # variance diag(K** - K* K^-1 K*') rides along for se.fit.
  for (gi in lp$gps %||% list()) {
    Xc <- do.call(cbind, lapply(gi$exprs, function(ex) {
      as.numeric(eval(ex, newdata, env))
    }))
    bk <- fit$frame$re_blocks[[gi$block_id]]
    extra_var <- NULL
    if (gi$type == "hsgp") {
      # dmax/center/L are the fitted scaling, so an in-sample newdata
      # row rebuilds its fitted basis row bit for bit
      Xr <- hsgp_basis(sweep(Xc / gi$dmax, 2, gi$center), gi$omega, gi$L)
    } else {
      pos <- gi$positions
      j <- match(pos_rowkey(Xc), pos_rowkey(pos))
      if (anyNA(j)) {
        th <- fit$estimates$theta[bk$theta_idx]
        K <- unname(covstruct_registry$gp$vcov(th, bk))
        Ks <- gp_cross_cov(th, bk, Xc, pos)
        Xr <- t(solve(K, t(Ks)))
        # observed rows reduce to indicators (K* row = K row), so their
        # conditional variance vanishes; the clamp absorbs roundoff
        kss <- exp(2 * th[1]) * (1 + 1e-6)
        extra_var <- pmax(kss - rowSums(Xr * Ks), 0)
      } else {
        # every row observed: exact indicator fast path
        Xr <- as.matrix(Matrix::sparseMatrix(i = seq_len(nrow(Xc)),
                                             j = j, x = 1,
                                             dims = c(nrow(Xc),
                                                      nrow(pos))))
      }
    }
    if (is.null(dim(Xr))) Xr <- matrix(Xr, nrow = nrow(Xc))
    sm_parts[[length(sm_parts) + 1L]] <- list(
      bk = bk,
      Xr = Xr,
      extra_var = extra_var
    )
  }

  nd_mult <- function(mult_expr) {
    if (is.null(mult_expr)) return(1)
    # same type gate as fit time, so a newdata column that changed type
    # reports the type rather than a downstream all-NA column [brms#1828]
    m <- check_special_mult(eval(mult_expr, newdata, env), mult_expr, "mo/mi")
    if (anyNA(m)) {
      stop("Interaction multiplier '", deparse1(mult_expr),
           "' has missing values in newdata", call. = FALSE)
    }
    m
  }
  for (mi in lp$mo %||% list()) {
    v <- mo_col_values(fit, mi, mo_codes(fit, lp, mi, newdata)) *
      nd_mult(mi$mult_expr)
    X <- cbind(X, matrix(v, ncol = 1, dimnames = list(NULL, mi$label)))
  }
  for (mt in lp$mi %||% list()) {
    v <- newdata[[mt$var]]
    if (is.null(v) || anyNA(v)) {
      stop("mi(", mt$var, "): newdata must supply complete values",
           call. = FALSE)
    }
    X <- cbind(X, matrix(as.numeric(v) * nd_mult(mt$mult_expr),
                         ncol = 1, dimnames = list(NULL, mt$label)))
  }

  if (!use_re) {
    # smooth wiggly parts are part of the curve, not group-level effects:
    # they stay in even for population-level predictions
    return(list(X = X, off = off, re_parts = list(), sm_parts = sm_parts,
                nonest = nonest))
  }

  re_parts <- list()
  for (bk in fit$frame$re_blocks) {
    if (bk$covstruct %in% c("smooth", "gp", "hsgp")) next
    for (comp in bk$components) {
      if (comp$lp_key != linpred_key(lp$resp, lp$dpar)) next
      tt2 <- stats::terms(stats::as.formula(call("~", comp$bar[[2]]),
                                            env = env))
      tt2 <- patch_predvars(tt2, fit$frame$predvar_map)
      mf2 <- stats::model.frame(tt2, newdata, na.action = stats::na.pass,
                                xlev = xlev_for(lp$xlevels, tt2))
      mm <- stats::model.matrix(tt2, mf2)
      if (!identical(colnames(mm), comp$cnms)) {
        stop("Random-effect design for `", comp$label, "` does not match ",
             "the fitted model (columns: ",
             paste(colnames(mm), collapse = ", "), " vs ",
             paste(comp$cnms, collapse = ", "), ")", call. = FALSE)
      }
      gv <- as.character(eval(comp$bar[[3]], newdata, env))
      j <- match(gv, bk$levels)
      if (anyNA(j) && !allow_new_levels) {
        stop("New levels in grouping factor `", deparse1(comp$bar[[3]]),
             "`: ", paste(unique(gv[is.na(j)]), collapse = ", "),
             ". Use allow_new_levels = TRUE to predict them at the ",
             "population level", call. = FALSE)
      }
      re_parts[[length(re_parts) + 1L]] <- list(bk = bk, comp = comp,
                                                mm = mm, j = j)
    }
  }
  list(X = X, off = off, re_parts = re_parts, sm_parts = sm_parts,
       nonest = nonest)
}

# RE contribution to eta for one linear predictor, given the full
# coefficient-space vector (see coef_b).
re_eta <- function(re_parts, cvec, n) {
  eta <- numeric(n)
  for (rp in re_parts) {
    bk <- rp$bk
    B <- t(matrix(cvec[bk$c_idx], nrow = bk$dim))  # levels x D
    cols <- rp$comp$offset + seq_len(rp$comp$dim)
    contrib <- rowSums(rp$mm * B[rp$j, cols, drop = FALSE])
    # only unmatched LEVELS predict at the population value; an NA in
    # the RE design data itself must propagate, not silently zero
    contrib[is.na(rp$j)] <- 0
    eta <- eta + contrib
  }
  eta
}

# Smooth wiggly contribution to eta for one linear predictor.
sm_eta <- function(sm_parts, cvec) {
  eta <- 0
  for (sp in sm_parts) {
    eta <- eta + drop(sp$Xr %*% cvec[sp$bk$c_idx])
  }
  eta
}

# Numeric dpar values at the estimates (optionally with a supplied b),
# for the training data. Nested: out[[resp]][[dpar]].
eval_dpars <- function(fit, b = fit$estimates[["b"]]) {
  est <- fit$estimates
  if (!is.null(b)) b <- expand_b(fit$frame, b, est$theta)
  out <- list()
  for (lp in fit$frame$linpreds) {
    if (!is.null(lp$nl_body)) {
      ev <- c(out[[lp$resp]], lp$data_list)
      eta <- eval(lp$nl_body, ev, lp$nl_env)
      out[[lp$resp]][[lp$dpar]] <- lp$link$linkinv(eta)
      next
    }
    eta <- if (ncol(lp$X)) {
      drop(as.matrix(patch_mo_cols(fit, lp, lp$X) %*%
                       est[[lp$par]][lp$idx]))
    } else {
      numeric(fit$frame$n_obs)
    }
    if (!is.null(lp$Z)) {
      eta <- eta + as.numeric(lp$Z %*% b)
    }
    if (!is.null(lp$offset)) eta <- eta + lp$offset
    out[[lp$resp]][[lp$dpar]] <- lp$link$linkinv(eta)
  }
  out
}

# Sparse RE design (n x n_c) for the newdata delta method, columns at
# global coefficient-space positions.
re_design_matrix <- function(re_parts, n, q) {
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  for (rp in re_parts) {
    bk <- rp$bk
    D <- bk$dim
    ok <- which(!is.na(rp$j))
    for (k in seq_len(rp$comp$dim)) {
      ii <- c(ii, ok)
      jj <- c(jj, bk$c_idx[(rp$j[ok] - 1L) * D + rp$comp$offset + k])
      xx <- c(xx, rp$mm[ok, k])
    }
  }
  Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, q))
}

# Jacobians of the coefficient-space expansion for rr fits: d cvec/d b
# (sparse; identity except the rr blocks' loadings) and d cvec/d theta
# for the rr loading parameters (finite differences on expand_b).
rr_jacobians <- function(fit) {
  frame <- fit$frame
  est <- fit$estimates
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  th_cols <- list()
  for (bk in frame$re_blocks) {
    if (bk$covstruct == "rr") {
      L <- rr_loadings(est$theta[bk$theta_idx], bk$dim, bk$rank)
      for (l in seq_len(bk$n_levels)) {
        rows <- bk$c_idx[(l - 1L) * bk$dim + seq_len(bk$dim)]
        cols <- bk$b_idx[(l - 1L) * bk$rank + seq_len(bk$rank)]
        ii <- c(ii, rep(rows, times = bk$rank))
        jj <- c(jj, rep(cols, each = bk$dim))
        xx <- c(xx, as.vector(L))
      }
      for (j in bk$theta_idx) {
        h <- 1e-6 * max(1, abs(est$theta[j]))
        tp <- est$theta; tp[j] <- tp[j] + h
        tn <- est$theta; tn[j] <- tn[j] - h
        dvec <- (expand_b(frame, est[["b"]], tp) -
                   expand_b(frame, est[["b"]], tn)) / (2 * h)
        th_cols[[length(th_cols) + 1L]] <- list(j = j, dvec = dvec)
      }
    } else {
      ii <- c(ii, bk$c_idx)
      jj <- c(jj, bk$b_idx)
      xx <- c(xx, rep(1, length(bk$b_idx)))
    }
  }
  list(Jb = Matrix::sparseMatrix(i = ii, j = jj, x = xx,
                                 dims = c(frame$n_c,
                                          length(est[["b"]]))),
       th_cols = th_cols)
}

# The single response of a univariate fit, or an informative error.
uni_resp <- function(fit, what) {
  if (length(fit$spec$responses) > 1) {
    stop(what, " is not supported yet for multivariate fits",
         call. = FALSE)
  }
  fit$spec$responses[[1]]
}

# Whether the family's response mean is just the mu dpar. All the
# built-in identity-mean families spell their mean_fn as `dpars$mu`, so
# the structural check is exact for them and conservative (general
# path) for custom families.
mean_is_mu <- function(fam) {
  is.null(fam$post$mean_fn) ||
    identical(body(fam$post$mean_fn), quote(dpars$mu))
}

# Addition-term values (trials, se, ...) re-evaluated on new data, for
# the expected-response prediction path. Terms the family mean cannot
# use (censoring, truncation, structural flags) are skipped; trials and
# se must evaluate because omitting them silently changes the mean.
aterms_for_newdata <- function(rspec, newdata) {
  skip <- c("cens", "cens_y2", "trunc_lb", "trunc_ub", "se_sigma",
            "mi", "mi_sd", "weights")
  av <- list()
  for (nm in setdiff(names(rspec$aterms), skip)) {
    v <- tryCatch(
      as.numeric(eval(rspec$aterms[[nm]], newdata, rspec$formula_env)),
      error = function(e) NULL
    )
    if (is.null(v) && nm %in% c("trials", "se")) {
      stop("Addition term ", nm, "(",
           deparse1(rspec$aterms[[nm]]), ") could not be evaluated on ",
           "newdata; supply the variable or use type = \"conditional\"",
           call. = FALSE)
    }
    if (!is.null(v)) av[[nm]] <- v
  }
  if (!is.null(rspec$aterms$se_sigma)) {
    av$se_sigma <- rspec$aterms$se_sigma
  }
  av
}

# Expected response over all dpars: the family mean at predicted dpar
# values (fitted()'s convention, extended to newdata and re.form).
predict_mean_response <- function(fit, rspec, newdata, re.form,
                                  allow_new_levels) {
  fam <- rspec$family
  rn <- rspec$resp_name
  if (is.null(newdata) && is.null(re.form)) {
    # exactly fitted(): dpars at the estimates, conditional on the modes
    dp <- eval_dpars(fit)[[rn]]
    out <- fam$post$mean_fn(dp, fit$frame$aterm_values[[rn]])
    return(napred(fit, out))
  }
  dp <- list()
  for (dnm in names(rspec$dpars)) {
    dp[[dnm]] <- as.vector(predict(fit, newdata = newdata, dpar = dnm,
                                   resp = rn, re.form = re.form,
                                   type = "response",
                                   allow_new_levels = allow_new_levels))
  }
  av <- if (is.null(newdata)) {
    # in-sample dpar predictions come back napredict-ed; pad the
    # per-observation aterm values the same way (a no-op under na.omit)
    lapply(fit$frame$aterm_values[[rn]], function(v) {
      if (is.numeric(v) && length(v) == fit$frame$n_obs) {
        napred(fit, v)
      } else {
        v
      }
    })
  } else {
    aterms_for_newdata(rspec, newdata)
  }
  fam$post$mean_fn(dp, av)
}

#' Predictions from a frmtmb fit
#'
#' @param object A `frmtmb_fit`.
#' @param newdata Optional data frame to predict on. Defaults to the
#'   training data.
#' @param type `"link"` for the linear predictor, `"response"` for the
#'   expected response (which equals [fitted()] on the training data;
#'   for zero-inflated, hurdle, and similar families this is the
#'   response mean, not the `mu` dpar). When `dpar` is given,
#'   `"response"` is that dpar on its natural scale. The glmmTMB
#'   spellings `"conditional"` (the `mu` dpar on its natural scale),
#'   `"zprob"`/`"zlink"` (the zero-inflation/hurdle probability on the
#'   response/link scale), and `"disp"` (the dispersion dpar) are
#'   accepted as aliases.
#' @param dpar Which distributional parameter to predict; defaults to the
#'   family's first location parameter (`"mu"` for most families).
#' @param resp For multivariate fits: which response to predict (defaults
#'   to the first).
#' @param re.form `NULL` (default) includes random effects; `NA` or `~0`
#'   gives population-level predictions.
#' @param se.fit If `TRUE`, return a list with elements `fit` and `se.fit`
#'   (delta-method standard errors accounting for fixed-effect and
#'   random-effect uncertainty). Exact `gp()` terms predict unseen
#'   positions by kriging: the conditional mean at the fitted kernel,
#'   with the GP conditional variance added to the standard errors.
#' @param allow_new_levels Predict unseen grouping-factor levels at the
#'   population level instead of erroring.
#' @param ... Unused.
#' @details
#' When the fixed-effect design was rank deficient, the aliased columns
#' were dropped at fit time and some coefficient combinations are not
#' estimable. Rows of `newdata` that load on a dropped direction get
#' `NA` (and `NA` standard errors), with one warning naming the dropped
#' columns; every other row is unaffected. The test is the one
#' [stats::predict.lm()] uses: a row is non-estimable when it is not
#' orthogonal to the null space of the fitted design, up to a relative
#' tolerance of `1e-8`. Two limits follow. It is a numerical test, so
#' near-aliased designs sit on a threshold rather than a clean
#' yes/no. And it covers the parametric fixed-effect block only:
#' smooth null-space, `gp()`, `mo()` and `mi()` columns are appended
#' after the rank check and are never dropped.
#' @return A numeric vector, or a list when `se.fit = TRUE`.
#' @export
predict.frmtmb_fit <- function(object, newdata = NULL,
                               type = c("link", "response",
                                        "conditional", "zprob", "zlink",
                                        "disp"),
                               dpar = NULL, resp = NULL, re.form = NULL,
                               se.fit = FALSE,
                               allow_new_levels = FALSE, ...) {
  dots <- list(...)
  # lme4/glmmTMB spell allow_new_levels with dots; accept it silently
  if ("allow.new.levels" %in% names(dots)) {
    allow_new_levels <- isTRUE(dots[["allow.new.levels"]])
    dots[["allow.new.levels"]] <- NULL
  }
  if (length(dots)) {
    warning("ignoring unknown arguments to predict(): ",
            paste(names(dots), collapse = ", "), call. = FALSE)
  }
  type <- match.arg(type)
  use_re <- is.null(re.form) ||
    (inherits(re.form, "formula") && !identical(deparse1(re.form[[2]]), "0"))
  if (!is.null(re.form) && !inherits(re.form, "formula")) use_re <- FALSE

  resp <- resp %||% names(object$spec$responses)[1]
  rspec <- object$spec$responses[[resp]]
  if (is.null(rspec)) {
    stop("Unknown response: '", resp, "'. Available: ",
         paste(names(object$spec$responses), collapse = ", "),
         call. = FALSE)
  }
  # glmmTMB type aliases resolve to a dpar plus scale
  if (type %in% c("zprob", "zlink", "disp")) {
    if (!is.null(dpar)) {
      stop("type = '", type, "' selects its own dpar; drop dpar =",
           call. = FALSE)
    }
    dpar <- switch(type,
      zprob = ,
      zlink = intersect(c("zi", "hu"), names(rspec$dpars))[1],
      disp = intersect(c("sigma", "shape", "phi"),
                       names(rspec$dpars))[1]
    )
    if (is.na(dpar)) {
      stop("type = '", type, "' needs a family with a ",
           if (type == "disp") "dispersion" else "zero-inflation/hurdle",
           " parameter; family '", rspec$family$family, "' has none",
           call. = FALSE)
    }
    type <- if (type == "zlink") "link" else "response"
  } else if (type == "conditional") {
    type <- "response"   # the conditional mean is the mu dpar
  } else if (type == "response" && is.null(dpar) &&
             !mean_is_mu(rspec$family)) {
    # the response mean is not the mu dpar (zi, hurdle, lognormal,
    # trials-binomial, ...): "response" means the expected response,
    # the fitted()/glmmTMB/brms-epred convention. Per-dpar values stay
    # available through dpar = or type = "conditional".
    if (se.fit) {
      stop("se.fit is not available for the expected response of ",
           "family '", rspec$family$family, "'; use type = ",
           "\"conditional\" (the mu dpar) or a parametric bootstrap",
           call. = FALSE)
    }
    return(predict_mean_response(object, rspec, newdata, re.form,
                                 allow_new_levels))
  }
  dpar <- dpar %||% if ("mu" %in% names(rspec$dpars)) "mu" else
    rspec$primary_dpars[1]
  key <- linpred_key(resp, dpar)
  lp <- object$frame$linpreds[[key]]
  if (is.null(lp)) {
    stop("Unknown dpar: '", dpar, "' for response '", resp,
         "'. Available: ", paste(names(rspec$dpars), collapse = ", "),
         call. = FALSE)
  }

  if (!is.null(lp$nl_body)) {
    if (se.fit) {
      stop("se.fit is not supported for the nonlinear predictor yet; ",
           "request the nonlinear parameters (dpar = '",
           rspec$nlpars[1], "', ...) instead", call. = FALSE)
    }
    vals <- list()
    for (np in rspec$nlpars) {
      vals[[np]] <- predict(object, newdata = newdata, dpar = np,
                            resp = resp, re.form = re.form,
                            allow_new_levels = allow_new_levels)
    }
    dl <- if (is.null(newdata)) {
      lp$data_list
    } else {
      lapply(stats::setNames(names(lp$data_list), names(lp$data_list)),
             function(v) {
               if (is.null(newdata[[v]])) {
                 stop("Variable '", v, "' missing from newdata",
                      call. = FALSE)
               }
               newdata[[v]]
             })
    }
    eta <- eval(lp$nl_body, c(vals, dl), lp$nl_env)
    out <- if (type == "response") lp$link$linkinv(eta) else eta
    return(if (is.null(newdata)) napred(object, out) else out)
  }

  est <- object$estimates
  re_parts <- list()
  sm_parts <- list()
  nonest <- FALSE
  sm_blocks <- Filter(function(bk) {
    bk$covstruct %in% c("smooth", "gp", "hsgp") &&
      any(vapply(bk$components, function(cp) cp$lp_key == key, TRUE))
  }, object$frame$re_blocks)

  if (is.null(newdata)) {
    X <- patch_mo_cols(object, lp, lp$X)
    off <- lp$offset
    n <- object$frame$n_obs
    eta <- drop(as.matrix(X %*% est[[lp$par]][lp$idx]))
    if (!is.null(lp$Z)) {
      cvec <- coef_b(object)
      if (use_re) {
        eta <- eta + as.numeric(lp$Z %*% cvec)
      } else if (length(sm_blocks)) {
        # population-level: drop group effects but keep the smooth curve
        for (bk in sm_blocks) {
          eta <- eta + as.numeric(lp$Z[, bk$c_idx, drop = FALSE] %*%
                                    cvec[bk$c_idx])
        }
      }
    }
  } else {
    pd <- pred_design(object, lp, newdata, allow_new_levels,
                      use_re = use_re)
    X <- pd$X
    off <- pd$off
    re_parts <- pd$re_parts
    sm_parts <- pd$sm_parts
    nonest <- pd$nonest
    n <- nrow(X)
    eta <- drop(as.matrix(X %*% est[[lp$par]][lp$idx]))
    cvec <- coef_b(object)
    if (use_re && length(re_parts)) {
      eta <- eta + re_eta(re_parts, cvec, n)
    }
    if (length(sm_parts)) {
      eta <- eta + sm_eta(sm_parts, cvec)
    }
  }
  if (!is.null(off)) eta <- eta + off
  # A non-estimable row has no defined fixed-effect part, so returning
  # the remaining (random-effect, smooth) contributions alone would look
  # like a valid prediction. One warning per call names the culprits.
  if (any(nonest)) {
    warning("Rank-deficient fit: ", sum(nonest), " row(s) of newdata are ",
            "not estimable because they load on the dropped column(s) ",
            paste(lp$dropped_colnames, collapse = ", "),
            "; predicting NA there", call. = FALSE)
    eta[nonest] <- NA_real_
  }

  if (!se.fit) {
    out <- if (type == "response") lp$link$linkinv(eta) else eta
    return(if (is.null(newdata)) napred(object, out) else out)
  }

  # Delta method: var(eta) = A V A' over the estimated coefficients (and b
  # when random effects are included). For rr fits the Z matrices span
  # the coefficient space, so the b columns go through the Jacobian of
  # the loadings expansion, and the rr loading parameters (theta)
  # contribute their own columns.
  has_rr <- isTRUE(object$frame$has_rr)
  rrj <- if (has_rr) rr_jacobians(object)
  jc <- get_joint_cov(object)
  rn <- jc$names
  add_b_cols <- function(A, coef_pos, Zc, b_pos, th_pos) {
    if (has_rr) {
      A <- Matrix::cbind2(A, Zc %*% rrj$Jb)
      coef_pos <- c(coef_pos, b_pos)
      for (tc in rrj$th_cols) {
        A <- Matrix::cbind2(A, Zc %*% tc$dvec)
        coef_pos <- c(coef_pos, th_pos[tc$j])
      }
    } else {
      A <- Matrix::cbind2(A, Zc)
      coef_pos <- c(coef_pos, b_pos)
    }
    list(A = A, coef_pos = coef_pos)
  }
  if (lp$par == "beta") {
    coef_pos <- which(rn == "beta")[lp$idx]
    A <- X
  } else {
    tpl_len <- length(object$frame$par_template$betad)
    est_rank <- match(lp$idx,
                      setdiff(seq_len(tpl_len),
                              object$frame$betad_fixed_idx))
    keep <- !is.na(est_rank)
    coef_pos <- which(rn == "betad")[est_rank[keep]]
    A <- X[, keep, drop = FALSE]
  }
  if (length(object$frame$re_blocks)) {
    b_pos <- which(rn == "b")
    th_pos <- which(rn == "theta")
    if (is.null(newdata)) {
      if (use_re && !is.null(lp$Z)) {
        upd <- add_b_cols(A, coef_pos, lp$Z, b_pos, th_pos)
        A <- upd$A
        coef_pos <- upd$coef_pos
      } else if (!use_re && length(sm_blocks) && !is.null(lp$Z)) {
        for (bk in sm_blocks) {
          A <- Matrix::cbind2(A, lp$Z[, bk$c_idx, drop = FALSE])
          coef_pos <- c(coef_pos, b_pos[bk$b_idx])
        }
      }
    } else {
      if (use_re && length(re_parts)) {
        Zn <- re_design_matrix(re_parts, n,
                               object$frame$n_c %||% length(est[["b"]]))
        upd <- add_b_cols(A, coef_pos, Zn, b_pos, th_pos)
        A <- upd$A
        coef_pos <- upd$coef_pos
      }
      for (sp in sm_parts) {
        A <- Matrix::cbind2(A, sp$Xr)
        coef_pos <- c(coef_pos, b_pos[sp$bk$b_idx])
      }
    }
  }
  V <- jc$V[coef_pos, coef_pos, drop = FALSE]
  A <- as.matrix(A)
  var_eta <- pmax(rowSums((A %*% V) * A), 0)
  # new grouping levels (allow_new_levels) contribute their block's
  # marginal variance: the population-prediction-interval convention.
  # For |ID|-merged blocks this is the JOINT block's slice for this
  # component.
  if (use_re && length(re_parts)) {
    th <- object$estimates$theta
    for (rp in re_parts) {
      nas <- which(is.na(rp$j))
      if (!length(nas)) next
      bk <- rp$bk
      if (bk$covstruct %in% c("gr_cov", "gr_prec")) next  # levels ARE the structure
      S <- covstruct_registry[[bk$covstruct]]$vcov(th[bk$theta_idx], bk)
      cols <- rp$comp$offset + seq_len(rp$comp$dim)
      Scc <- S[cols, cols, drop = FALSE]
      mmn <- rp$mm[nas, , drop = FALSE]
      var_eta[nas] <- var_eta[nas] + rowSums((mmn %*% Scc) * mmn)
    }
  }
  # exact-gp kriging at unseen positions: the GP's own conditional
  # variance adds to the delta-method variance (same convention as the
  # new-level marginal variance above); zero at observed positions
  for (sp in sm_parts) {
    if (!is.null(sp$extra_var)) var_eta <- var_eta + sp$extra_var
  }
  se_eta <- sqrt(var_eta)
  # the kept columns still have a finite variance, but it is not the
  # standard error of anything the fit estimates
  if (any(nonest)) se_eta[nonest] <- NA_real_

  out <- if (type == "response") {
    list(fit = lp$link$linkinv(eta),
         se.fit = abs(lp$link$mu_eta(eta)) * se_eta)
  } else {
    list(fit = eta, se.fit = se_eta)
  }
  if (is.null(newdata)) {
    out$fit <- napred(object, out$fit)
    out$se.fit <- napred(object, out$se.fit)
  }
  out
}

#' @export
model.frame.frmtmb_fit <- function(formula, ...) {
  # the stored combined frame survives even when the caller's data
  # environment is gone (lme4 test-formulaEval.R bug class)
  formula$frame$data_frame
}

# Reinsert NAs for na.exclude fits (napredict is a no-op for na.omit).
napred <- function(fit, x) {
  stats::napredict(fit$frame$na_action, x)
}

#' @export
fitted.frmtmb_fit <- function(object, ...) {
  rspec <- uni_resp(object, "fitted()")
  dp <- eval_dpars(object)[[rspec$resp_name]]
  if (!"mu" %in% names(dp) && is.null(rspec$family$post$mean_fn)) {
    stop("fitted() is not defined for family '", rspec$family$family, "'",
         call. = FALSE)
  }
  fam <- rspec$family
  out <- if (!is.null(fam$post$mean_fn)) {
    fam$post$mean_fn(dp, object$frame$aterm_values[[rspec$resp_name]])
  } else {
    dp$mu
  }
  napred(object, out)
}

#' Residuals from a frmtmb fit
#'
#' `"osa"` gives one-step-ahead (conditional quantile) residuals via
#' [TMB::oneStepPredict()]: standard-normal under a correctly specified
#' model, valid under correlated observations where pearson residuals
#' mislead.
#'
#' @param object A `frmtmb_fit`.
#' @param type `"response"`, `"pearson"`, or `"osa"`.
#' @param osa_method Method for [TMB::oneStepPredict()]; defaults to
#'   `"fullGaussian"` for gaussian models and `"oneStepGeneric"`
#'   otherwise.
#' @param ... For `type = "osa"`: passed to [TMB::oneStepPredict()].
#' @return A numeric vector.
#' @export
residuals.frmtmb_fit <- function(object, type = c("response", "pearson",
                                                  "osa"),
                                 osa_method = NULL, ...) {
  type <- match.arg(type)
  rspec <- uni_resp(object, "residuals()")
  fam <- rspec$family
  if (type == "osa") {
    method <- osa_method %||%
      if (identical(fam$family, "gaussian")) "fullGaussian" else
        "oneStepGeneric"
    args <- list(obj = object$obj, observation.name = ".frm_obs",
                 method = method, trace = FALSE, ...)
    if (method == "oneStepGeneric") {
      args$discrete <- identical(fam$type, "discrete")
      if (args$discrete && is.null(args$range)) args$range <- c(0, Inf)
    }
    osa <- do.call(RTMB::oneStepPredict, args)
    return(napred(object, osa$residual))
  }
  dp <- eval_dpars(object)[[rspec$resp_name]]
  mu <- if (!is.null(fam$post$mean_fn)) {
    fam$post$mean_fn(dp, object$frame$aterm_values[[rspec$resp_name]])
  } else {
    dp$mu
  }
  r <- object$frame$y[[rspec$resp_name]] - mu
  if (type == "pearson") {
    if (is.null(fam$post$var_fn)) {
      stop("Family '", fam$family, "' has no variance function; ",
           "pearson residuals are unavailable", call. = FALSE)
    }
    v <- fam$post$var_fn(dp, object$frame$aterm_values[[rspec$resp_name]])
    r <- r / sqrt(v)
  }
  napred(object, r)
}

# One draw of the full b vector from its estimated distribution N(0, Sigma).
draw_b <- function(fit) {
  th <- fit$estimates$theta
  b <- numeric(length(fit$estimates[["b"]] %||% numeric(0)))
  for (bk in fit$frame$re_blocks) {
    if (bk$covstruct == "rr") {
      # standard-normal factors; eval_dpars expands them via loadings
      b[bk$b_idx] <- stats::rnorm(length(bk$b_idx))
      next
    }
    if (bk$covstruct == "gr_cov") {
      # correlation is across levels, not within them
      S <- covstruct_registry$gr_cov$vcov(th[bk$theta_idx], bk)
      K <- kronecker(bk$aux_A, S)
      b[bk$b_idx] <- drop(crossprod(chol(K),
                                    stats::rnorm(nrow(K))))
      next
    }
    if (bk$covstruct == "gr_prec") {
      # x = sd * U^-1 z with U'U = Q has covariance sd^2 Q^-1
      U <- Matrix::chol(bk$aux_Q)
      z <- stats::rnorm(bk$n_levels)
      b[bk$b_idx] <- exp(th[bk$theta_idx]) *
        as.vector(Matrix::solve(U, z))
      next
    }
    V <- covstruct_registry[[bk$covstruct]]$vcov(th[bk$theta_idx], bk)
    L <- chol(V)
    U <- matrix(stats::rnorm(bk$n_levels * bk$dim), bk$n_levels) %*% L
    b[bk$b_idx] <- as.vector(t(U))   # level-major
  }
  b
}

#' Simulate responses from a frmtmb fit
#'
#' @param object A `frmtmb_fit`.
#' @param nsim Number of simulated response vectors.
#' @param seed Optional RNG seed. Follows the [stats::simulate()]
#'   contract: the global RNG state is restored afterwards, and the
#'   seed used is attached as the `"seed"` attribute.
#' @param re.form `NULL` (default) conditions on the estimated random
#'   effects; `NA` redraws them from their estimated distribution
#'   (marginal simulation).
#' @param ... Unused.
#' @return A data frame with `nsim` columns and a `"seed"` attribute.
#' @export
simulate.frmtmb_fit <- function(object, nsim = 1, seed = NULL,
                                re.form = NULL, ...) {
  # the stats::simulate seed contract (as in simulate.lm)
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1)
  }
  if (is.null(seed)) {
    rng_state <- get(".Random.seed", envir = globalenv())
  } else {
    saved_seed <- get(".Random.seed", envir = globalenv())
    set.seed(seed)
    rng_state <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", saved_seed, envir = globalenv()))
  }
  rspec <- uni_resp(object, "simulate()")
  fam <- rspec$family
  if (is.null(fam$sim)) {
    stop("Family '", fam$family, "' has no simulator yet", call. = FALSE)
  }
  mg <- object$frame$mix_g[[rspec$resp_name]]
  marginal <- !is.null(re.form) && !inherits(re.form, "formula") &&
    is.na(re.form)
  n <- stats::nobs(object)
  out <- vector("list", nsim)
  av <- object$frame$aterm_values[[rspec$resp_name]]
  for (s in seq_len(nsim)) {
    b_use <- if (marginal && length(object$frame$re_blocks)) {
      draw_b(object)
    } else {
      object$estimates[["b"]]
    }
    dp <- eval_dpars(object, b = b_use)[[rspec$resp_name]]
    out[[s]] <- if (is.null(mg)) {
      fam$sim(dp, av, n)
    } else {
      # latent-class mixture: one class draw per group, then each
      # observation simulates from its group's component
      lps <- fam$mix$log_pi(dp)
      Pg <- vapply(lps, function(l) {
        exp(rep(l, length.out = n)[mg$first])
      }, numeric(length(mg$first)))
      kg <- vapply(seq_len(nrow(Pg)), function(g_) {
        sample.int(fam$mix$K, 1L, prob = Pg[g_, ])
      }, integer(1))
      kk <- kg[mg$gindex]
      ys <- numeric(n)
      for (k in seq_len(fam$mix$K)) {
        idx <- which(kk == k)
        if (length(idx)) {
          dk <- lapply(fam$mix$comp_dpars(dp, k), function(v) {
            rep(v, length.out = n)[idx]
          })
          ys[idx] <- fam$mix$comp_sim(dk, av, length(idx), k)
        }
      }
      ys
    }
  }
  names(out) <- paste0("sim_", seq_len(nsim))
  out <- as.data.frame(out)
  attr(out, "seed") <- rng_state
  out
}
