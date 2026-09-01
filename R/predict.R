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
    # same degradation vcov() uses: a singular joint precision gives NaN
    # standard errors and one warning naming diagnose(), not a raw
    # LAPACK message from deep inside predict()
    V <- as.matrix(solve_joint_precision(Q, cache))
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
      gvr <- eval(comp$bar[[3]], newdata, env)
      # an spde block's levels are mesh ROW NUMBERS, so the node has to
      # be read as a number here too: as.character() on a double would
      # spell node 100000 as "1e+05" and lose the column
      gv <- if (bk$covstruct == "spde") spde_node_labels(gvr) else {
        as.character(gvr)
      }
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

# Whether a response carries trunc() bounds, from the spec (the stored
# expressions) rather than from evaluated values.
has_trunc <- function(rspec) {
  any(c("trunc_lb", "trunc_ub") %in% names(rspec$aterms))
}

# Addition-term values (trials, se, trunc bounds, ...) re-evaluated on
# new data, for the expected-response prediction path. Terms the family
# mean cannot use (censoring, structural flags) are skipped; trials, se
# and the truncation bounds must evaluate because omitting them silently
# changes the mean.
# An addition term as the user wrote it, for error messages.
aterm_label <- function(nm, ex) {
  switch(nm,
    trunc_lb = paste0("trunc(lb = ", deparse1(ex), ")"),
    trunc_ub = paste0("trunc(ub = ", deparse1(ex), ")"),
    # vint(a, b) is stored one argument per aterm as vint1, vint2, ...
    paste0(sub("[0-9]+$", "", nm), "(", deparse1(ex), ")")
  )
}

# A custom family's lpdf/mean_fn reads its vint()/vreal() payload out of
# aterms, so an omitted one is not a missing covariate but a missing
# argument: the family returns a zero-length prediction.
is_custom_data_aterm <- function(nm) grepl("^v(int|real)[0-9]+$", nm)

aterms_for_newdata <- function(rspec, newdata) {
  skip <- c("cens", "cens_y2", "se_sigma", "mi", "mi_sd", "weights")
  need <- c("trials", "se", "trunc_lb", "trunc_ub")
  nd_n <- nrow(newdata)
  av <- list()
  for (nm in setdiff(names(rspec$aterms), skip)) {
    ex <- rspec$aterms[[nm]]
    v <- tryCatch(
      as.numeric(eval(ex, newdata, rspec$formula_env)),
      error = function(e) NULL
    )
    # a bound the model frame supplied but newdata did not can still
    # resolve in the formula environment, to the FITTED rows; the length
    # check catches that rather than silently pairing the wrong bounds
    if (!is.null(v) && !is.null(nd_n) && !length(v) %in% c(1L, nd_n)) {
      v <- NULL
    }
    if (is.null(v)) {
      label <- aterm_label(nm, ex)
      missed <- setdiff(all.vars(ex), names(newdata))
      if (nm %in% need || is_custom_data_aterm(nm)) {
        stop("Addition term ", label, " could not be evaluated on ",
             "newdata",
             if (length(missed)) {
               paste0(": newdata has no column ",
                      paste(missed, collapse = ", "))
             } else "",
             "; supply the variable or use type = \"conditional\"",
             call. = FALSE)
      }
      # anything else is dropped, but never silently: an aterm the
      # family reads and this function omits is a wrong prediction
      warning("Addition term ", label, " could not be evaluated on ",
              "newdata and is omitted from the prediction",
              if (length(missed)) {
                paste0(" (newdata has no column ",
                       paste(missed, collapse = ", "), ")")
              } else "",
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
    out <- response_mean(fam, dp, fit$frame$aterm_values[[rn]])
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
  response_mean(fam, dp, av)
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
#' @section Truncated responses:
#' For a response with `trunc()` bounds, `type = "response"` (and
#' [fitted()]) report the truncated mean `E[Y | lb <= Y <= ub]`, matching
#' the likelihood the model was fitted with. Predictions of a
#' distributional parameter (`type = "link"`, `dpar = `, or
#' `type = "conditional"`) stay **untruncated**: they are statements
#' about the latent parameter, not about the observed, truncated
#' response. Bounds are re-evaluated on `newdata` the same way `trials()`
#' and `se()` are: a literal bound carries over unchanged, and a bound
#' given as a variable must be a column of `newdata` of the right length.
#' @param resp For multivariate fits: which response to predict (defaults
#'   to the first).
#' @param re.form `NULL` (default) includes random effects; `NA` or `~0`
#'   gives population-level predictions.
#' @param se.fit If `TRUE`, return a list with elements `fit` and `se.fit`
#'   (delta-method standard errors accounting for fixed-effect and
#'   random-effect uncertainty). Exact `gp()` terms predict unseen
#'   positions by kriging: the conditional mean at the fitted kernel,
#'   with the GP conditional variance added to the standard errors.
#' @section Standard errors of the expected response:
#' For a family whose mean is the `mu` dpar, `se.fit` on
#' `type = "response"` is the usual one-predictor delta method:
#' `|dmu/deta| * se(eta)`.
#'
#' When the mean is a function of several dpars (zero-inflated and
#' hurdle families, `lognormal`, a `trials()` binomial, or any
#' `trunc()`ed response), the delta method runs jointly over every
#' dpar's linear predictor: `se^2 = g' V g`, where row `i` of `g`
#' stacks `dm_i/deta_k` times the design row of predictor `k`, and `V`
#' is the joint covariance of all the coefficients (`vcov()`'s
#' `jointPrecision` block, so the cross-predictor covariances and the
#' shared random-effect block are included). The gradients
#' `dm/deta_k` are central differences of the family mean, taken one
#' predictor at a time with a relative step.
#'
#' Random effects enter conditional on their modes, the same convention
#' `se.fit` uses for the linear predictor. Unseen grouping levels
#' (`allow_new_levels = TRUE`) add their block's marginal variance,
#' propagated through the same gradients.
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
             (!mean_is_mu(rspec$family) || has_trunc(rspec))) {
    # the response mean is not the mu dpar (zi, hurdle, lognormal,
    # trials-binomial, ...), or the response is truncated so the
    # expected response is E[Y | lb <= Y <= ub]: "response" means the
    # expected response, the fitted()/glmmTMB/brms-epred convention.
    # Per-dpar values stay available through dpar = or
    # type = "conditional".
    if (se.fit) {
      return(predict_mean_se(object, rspec, newdata, use_re,
                             allow_new_levels))
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

  ed <- lp_eta_design(object, lp, newdata, use_re, allow_new_levels)
  eta <- ed$eta
  n <- ed$n

  if (!se.fit) {
    out <- if (type == "response") lp$link$linkinv(eta) else eta
    return(if (is.null(newdata)) napred(object, out) else out)
  }

  has_rr <- isTRUE(object$frame$has_rr)
  rrj <- if (has_rr) rr_jacobians(object)
  jc <- get_joint_cov(object)
  da <- lp_delta_A(object, lp, ed, newdata, use_re, jc, has_rr, rrj)
  V <- jc$V[da$coef_pos, da$coef_pos, drop = FALSE]
  A <- as.matrix(da$A)
  var_eta <- pmax(rowSums((A %*% V) * A), 0)
  ev <- lp_extra_var(object, ed, use_re)
  for (nl in ev$new_levels) {
    Scc <- nl$S[nl$cols, nl$cols, drop = FALSE]
    mmn <- nl$mm[nl$nas, , drop = FALSE]
    var_eta[nl$nas] <- var_eta[nl$nas] + rowSums((mmn %*% Scc) * mmn)
  }
  for (gv in ev$gp) var_eta <- var_eta + gv
  se_eta <- sqrt(var_eta)
  # the kept columns still have a finite variance, but it is not the
  # standard error of anything the fit estimates
  if (any(ed$nonest)) se_eta[ed$nonest] <- NA_real_

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

# eta and the design pieces of one linear predictor, in sample or on
# newdata. Shared by predict() and by the joint delta method for the
# expected response, so both see exactly the same eta.
lp_eta_design <- function(object, lp, newdata, use_re, allow_new_levels) {
  est <- object$estimates
  key <- linpred_key(lp$resp, lp$dpar)
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
  list(eta = eta, X = X, off = off, re_parts = re_parts,
       sm_parts = sm_parts, sm_blocks = sm_blocks, nonest = nonest, n = n)
}

# A standard error that leaves out a variance component is only honest
# if it says so, so this is a warning rather than a note in the docs.
warn_modes_conditional_se <- function() {
  warning("The fitted objective marginalizes the random effects ",
          "(quadrature = TRUE), so its covariance carries no ",
          "random-effect block: se.fit is conditional on the ",
          "conditional modes and omits random-effect uncertainty. ",
          "Refit with quadrature = FALSE for the full delta method",
          call. = FALSE)
  invisible(NULL)
}

# Delta method: var(eta) = A V A' over the estimated coefficients (and b
# when random effects are included). For rr fits the Z matrices span the
# coefficient space, so the b columns go through the Jacobian of the
# loadings expansion, and the rr loading parameters (theta) contribute
# their own columns. Returns A and the positions of its columns in the
# joint covariance, so several linear predictors can be combined.
lp_delta_A <- function(object, lp, ed, newdata, use_re, jc, has_rr, rrj) {
  est <- object$estimates
  rn <- jc$names
  X <- ed$X
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
    would_add <- if (is.null(newdata)) {
      !is.null(lp$Z) && (use_re || length(ed$sm_blocks) > 0L)
    } else {
      (use_re && length(ed$re_parts) > 0L) || length(ed$sm_parts) > 0L
    }
    if (!length(b_pos) && would_add) {
      # A marginalized (quadrature) objective has no b in its parameter
      # vector, so the joint covariance carries no random-effect rows to
      # pair the Z columns with. Adding them anyway made A wider than V
      # and the delta method died non-conformable. Report the standard
      # error conditional on the modes instead, and say so.
      warn_modes_conditional_se()
      return(list(A = A, coef_pos = coef_pos))
    }
    if (is.null(newdata)) {
      if (use_re && !is.null(lp$Z)) {
        upd <- add_b_cols(A, coef_pos, lp$Z, b_pos, th_pos)
        A <- upd$A
        coef_pos <- upd$coef_pos
      } else if (!use_re && length(ed$sm_blocks) && !is.null(lp$Z)) {
        for (bk in ed$sm_blocks) {
          A <- Matrix::cbind2(A, lp$Z[, bk$c_idx, drop = FALSE])
          coef_pos <- c(coef_pos, b_pos[bk$b_idx])
        }
      }
    } else {
      if (use_re && length(ed$re_parts)) {
        Zn <- re_design_matrix(ed$re_parts, ed$n,
                               object$frame$n_c %||% length(est[["b"]]))
        upd <- add_b_cols(A, coef_pos, Zn, b_pos, th_pos)
        A <- upd$A
        coef_pos <- upd$coef_pos
      }
      for (sp in ed$sm_parts) {
        A <- Matrix::cbind2(A, sp$Xr)
        coef_pos <- c(coef_pos, b_pos[sp$bk$b_idx])
      }
    }
  }
  list(A = A, coef_pos = coef_pos)
}

# Variance sources that are not coefficient uncertainty. A new grouping
# level (allow_new_levels) contributes its block's marginal variance -
# the population-prediction-interval convention; for |ID|-merged blocks
# that is the JOINT block's slice for this component. An exact gp() at
# an unseen position contributes the GP's own conditional variance
# (zero at observed positions).
lp_extra_var <- function(object, ed, use_re) {
  nl <- list()
  if (use_re && length(ed$re_parts)) {
    th <- object$estimates$theta
    for (rp in ed$re_parts) {
      nas <- which(is.na(rp$j))
      if (!length(nas)) next
      bk <- rp$bk
      # the levels ARE the structure there, so there is no marginal
      # variance to hand an unseen one
      if (bk$covstruct %in% c("gr_cov", "gr_prec", "car", "spde")) next
      nl[[length(nl) + 1L]] <- list(
        bk = bk,
        S = covstruct_registry[[bk$covstruct]]$vcov(th[bk$theta_idx], bk),
        cols = rp$comp$offset + seq_len(rp$comp$dim),
        mm = rp$mm, nas = nas
      )
    }
  }
  gp <- list()
  for (sp in ed$sm_parts) {
    if (!is.null(sp$extra_var)) gp[[length(gp) + 1L]] <- sp$extra_var
  }
  list(new_levels = nl, gp = gp)
}

# Central-difference gradient of the expected response with respect to
# one dpar's linear predictor. Analytic gradients exist for the simple
# mean_fn forms, but the family set (and the truncated means) is wide
# enough that one differencing rule beats a table of hand derivatives;
# every mean_fn is elementwise, so a whole column of the Jacobian costs
# two evaluations. The step is relative so it survives both tiny and
# large etas.
mean_eta_grad <- function(fam, dp, av, dnm, link, eta) {
  h <- 1e-5 * pmax(1, abs(eta))
  dp_hi <- dp
  dp_lo <- dp
  dp_hi[[dnm]] <- link$linkinv(eta + h)
  dp_lo[[dnm]] <- link$linkinv(eta - h)
  (response_mean(fam, dp_hi, av) - response_mean(fam, dp_lo, av)) / (2 * h)
}

# Delta-method SEs for the expected response of a family whose mean is
# not the mu dpar (zero-inflation, hurdles, lognormal, trials-binomial,
# or any truncated response). The mean runs through EVERY dpar's linear
# predictor, so the gradient row stacks dm/deta_k times each predictor's
# own A matrix and the quadratic form is taken over the JOINT
# coefficient covariance - the cross-dpar covariances (and the shared b
# block) are part of the answer, not an afterthought. Random effects
# enter conditional on their modes, the same convention predict() uses
# for eta.
predict_mean_se <- function(object, rspec, newdata, use_re,
                            allow_new_levels) {
  fam <- rspec$family
  rnm <- rspec$resp_name
  jc <- get_joint_cov(object)
  has_rr <- isTRUE(object$frame$has_rr)
  rrj <- if (has_rr) rr_jacobians(object)
  dnames <- names(rspec$dpars)
  eds <- list()
  das <- list()
  evs <- list()
  dp <- list()
  for (dnm in dnames) {
    lp <- object$frame$linpreds[[linpred_key(rnm, dnm)]]
    if (!is.null(lp$nl_body)) {
      stop("se.fit is not supported for the nonlinear predictor yet; ",
           "request the nonlinear parameters (dpar = '",
           rspec$nlpars[1], "', ...) instead", call. = FALSE)
    }
    ed <- lp_eta_design(object, lp, newdata, use_re, allow_new_levels)
    eds[[dnm]] <- ed
    das[[dnm]] <- lp_delta_A(object, lp, ed, newdata, use_re, jc,
                             has_rr, rrj)
    evs[[dnm]] <- lp_extra_var(object, ed, use_re)
    dp[[dnm]] <- lp$link$linkinv(ed$eta)
  }
  n <- eds[[1L]]$n
  av <- if (is.null(newdata)) {
    object$frame$aterm_values[[rnm]]
  } else {
    aterms_for_newdata(rspec, newdata)
  }
  m <- response_mean(fam, dp, av)
  if (is.null(m) || !is.numeric(m) || !is.null(dim(m)) || length(m) != n) {
    stop("se.fit is not available for the expected response of family '",
         fam$family, "': its mean is not one number per observation",
         call. = FALSE)
  }

  grad <- list()
  for (dnm in dnames) {
    lp <- object$frame$linpreds[[linpred_key(rnm, dnm)]]
    grad[[dnm]] <- mean_eta_grad(fam, dp, av, dnm, lp$link, eds[[dnm]]$eta)
  }

  # one gradient row per observation over the union of coefficient
  # positions; predictors that share a coefficient (the b block) add
  pos_all <- unique(unlist(lapply(das, `[[`, "coef_pos")))
  G <- matrix(0, n, length(pos_all))
  for (dnm in dnames) {
    cp <- das[[dnm]]$coef_pos
    if (!length(cp)) next
    cols <- match(cp, pos_all)
    Ak <- as.matrix(das[[dnm]]$A) * grad[[dnm]]
    if (anyDuplicated(cols)) {
      for (j in seq_along(cols)) {
        G[, cols[j]] <- G[, cols[j]] + Ak[, j]
      }
    } else {
      G[, cols] <- G[, cols] + Ak
    }
  }
  V <- jc$V[pos_all, pos_all, drop = FALSE]
  var_m <- pmax(rowSums((G %*% V) * G), 0)

  # New grouping levels: a block whose components sit in several linear
  # predictors enters once, through the summed gradient over its own
  # component space, so the within-block cross-dpar covariance is kept.
  blocks <- list()
  for (dnm in dnames) {
    for (nl in evs[[dnm]]$new_levels) {
      bkey <- as.character(nl$bk$c_idx[1])
      B <- blocks[[bkey]] %||% list(S = nl$S,
                                    M = matrix(0, n, nl$bk$dim),
                                    rows = integer(0))
      B$M[nl$nas, nl$cols] <- B$M[nl$nas, nl$cols] +
        grad[[dnm]][nl$nas] * nl$mm[nl$nas, , drop = FALSE]
      B$rows <- union(B$rows, nl$nas)
      blocks[[bkey]] <- B
    }
  }
  for (B in blocks) {
    Mr <- B$M[B$rows, , drop = FALSE]
    var_m[B$rows] <- var_m[B$rows] + rowSums((Mr %*% B$S) * Mr)
  }
  for (dnm in dnames) {
    for (gv in evs[[dnm]]$gp) {
      var_m <- var_m + grad[[dnm]]^2 * gv
    }
  }

  se_m <- sqrt(var_m)
  names(se_m) <- names(m)   # the one-predictor path labels both alike
  nonest <- Reduce(`|`, lapply(eds, function(e) {
    rep(e$nonest, length.out = n)
  }))
  if (any(nonest)) se_m[nonest] <- NA_real_
  out <- list(fit = m, se.fit = se_m)
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
  out <- response_mean(fam, dp,
                       object$frame$aterm_values[[rspec$resp_name]])
  napred(object, out)
}

# Number of ordinal categories, from the threshold vector rather than
# from the data: the top category may be unobserved.
ordinal_ncat <- function(fit) {
  raw <- fit$estimates[["tau_raw"]]
  if (is.null(raw)) {
    rspec <- uni_resp(fit, "residuals()")
    return(max(fit$frame$y[[rspec$resp_name]]))
  }
  length(raw) + 1L
}

# OSA integration window and row split for a censored response, or NULL
# when nothing is censored.
#
# A censored row contributes a probability MASS: on the tape its
# likelihood no longer depends on the observation, so oneStepPredict
# either inverts a singular system (fullGaussian) or integrates a flat
# slice to infinity (oneStepGeneric). Both are real: the "observation"
# on such a row is an event, not a value, and has no one-step CDF.
# What IS well defined is the CDF of the uncensored rows conditional on
# the censoring events, which is what subset/conditional buys. Under
# type-I censoring (one censoring point per side) an uncensored row is
# exactly a draw that landed inside the window, so its PIT renormalizes
# on that window just as a trunc() fit's does.
osa_cens_domain <- function(av, y) {
  cen <- av[["cens"]]
  if (is.null(cen) || !any(cen != 0)) return(NULL)
  if (any(cen == 2)) {
    stop("residuals(type = \"osa\") does not support interval censoring ",
         "(cens code 2): an interval-censored row observes an event, not ",
         "a value, and the uncensored rows' observation window is then ",
         "not a single interval", call. = FALSE)
  }
  i_obs <- which(cen == 0)
  if (!length(i_obs)) {
    stop("residuals(type = \"osa\") needs at least one uncensored ",
         "observation", call. = FALSE)
  }
  point <- function(idx, side) {
    p <- unique(y[idx])
    if (length(p) > 1L) {
      stop("residuals(type = \"osa\") on a cens() fit needs one ",
           side, "-censoring point shared by every censored row ",
           "(type-I censoring); got ", length(p), " distinct points. ",
           "With row-varying censoring times the distribution of an ",
           "uncensored response is not identified without a model for ",
           "the censoring process. dharma_residuals() is not an ",
           "alternative: simulate() draws the LATENT uncensored ",
           "response, so its draws are not comparable with the ",
           "observed censored values, and simulate(censored = TRUE) ",
           "needs the same single censoring point this message is ",
           "about", call. = FALSE)
    }
    p
  }
  i_r <- which(cen == 1)
  i_l <- which(cen == -1)
  hi <- if (length(i_r)) point(i_r, "right") else Inf
  lo <- if (length(i_l)) point(i_l, "left") else -Inf
  if (any(y[i_obs] < lo) || any(y[i_obs] > hi)) {
    stop("residuals(type = \"osa\") on a cens() fit found uncensored ",
         "responses outside the censoring window [", lo, ", ", hi,
         "]; the censoring is not type-I and the one-step CDF has no ",
         "well-defined domain", call. = FALSE)
  }
  list(lo = lo, hi = hi, subset = i_obs, conditional = c(i_l, i_r))
}

#' Residuals from a frmtmb fit
#'
#' `"osa"` gives one-step-ahead (conditional quantile) residuals via
#' [TMB::oneStepPredict()]: standard-normal under a correctly specified
#' model, valid under correlated observations where pearson residuals
#' mislead.
#'
#' On a `trunc()`ed response, `"response"` residuals are taken against
#' the truncated mean `E[Y | lb <= Y <= ub]`. `"pearson"` divides by the
#' untruncated family variance, so it is conservative there. `"osa"`
#' builds its conditional CDF on `[lb, ub]` (see `osa_method`).
#'
#' On a `cens()`ed response, `"osa"` returns `NA` for every censored
#' row: what is observed there is an event (`Y > c`), not a value, and
#' an event has no one-step CDF. The uncensored rows get residuals
#' conditional on the censoring events, which needs one censoring point
#' per side (type-I censoring); row-varying censoring times and interval
#' censoring are refused. `dharma_residuals()` is not a substitute on a
#' censored fit, because [simulate.frmtmb_fit()] draws the latent
#' uncensored response by default (as brms's `posterior_predict()`
#' does) and those draws are not comparable with the observed censored
#' values; `simulate(censored = TRUE)` makes them comparable, but the
#' resulting point mass at each censoring point is not a distribution
#' DHARMa's rank transform can use.
#'
#' Ordinal responses use `"oneStepGeneric"` over the discrete support
#' `1..K`, which makes the residuals randomized quantile residuals.
#'
#' @section Deviance residuals:
#' `"deviance"` returns `sign(y - E[Y]) * sqrt(w * d)`, where the unit
#' deviance `d = 2 * (loglik of the saturated fit - loglik at the fitted
#' value)` is taken with the dispersion parameter held at its estimate,
#' and `w` is the `weights()` addition term (1 by default). For the
#' exponential-dispersion families this is the glm unit deviance, so a
#' fixed-effect fit reproduces `residuals(glm(...), type = "deviance")`
#' exactly.
#'
#' Supported families: `gaussian`, `poisson`, `binomial`, `bernoulli`,
#' `Gamma`, `exponential`, `inverse.gaussian`, `negbinomial`
#' (`nbinom2`), `nbinom1`, `geometric`, `beta`, and `tweedie`. Every
#' other family is refused: ordinal, mixture, multinomial, hurdle,
#' zero-inflated and location-shift families have no standard unit
#' deviance. `nbinom1` follows glmmTMB and evaluates the
#' negative-binomial size at the fitted row's `mu / phi`; letting the
#' size follow the saturated mean is not a deviance (the difference goes
#' negative). `trunc()`ed and `cens()`ed responses are refused as well,
#' because the fitted likelihood there is not the family's own density.
#'
#' A gaussian response with `se()` has no common dispersion for a raw
#' squared residual to be measured against, so the known variance
#' enters as a glm prior weight `sigma^2 / s_i^2` on top of `w`, where
#' `s_i` is the row's residual sd (the quantity `"pearson"` divides
#' by). Without `se()` that weight is 1 and nothing changes; `se(x)`
#' alone maps `sigma` out at 1, leaving the familiar `1 / se_i^2` of a
#' known-variance weighted fit.
#'
#' In a mixed model the residuals are conditional on the random-effect
#' modes, the glmmTMB convention: `E[Y]` is [fitted()], not the
#' population-level mean.
#'
#' `deviance(fit)` is unrelated: it stays `-2 * logLik(fit)` (the lme4
#' convention), which for a mixed model is the Laplace-approximated
#' marginal deviance and does **not** equal `sum(residuals(fit, type =
#' "deviance")^2)`.
#'
#' @param object A `frmtmb_fit`.
#' @param type `"response"`, `"pearson"`, `"deviance"`, or `"osa"`.
#' @param osa_method Method for [TMB::oneStepPredict()]; defaults to
#'   `"fullGaussian"` for gaussian models and `"oneStepGeneric"`
#'   otherwise. A truncated, censored or ordinal response always uses
#'   `"oneStepGeneric"` (a truncated gaussian is not gaussian) with the
#'   integration domain and discrete support taken from the `trunc()`
#'   bounds or the censoring window, which must then be the same for
#'   every row.
#' @param ... For `type = "osa"`: passed to [TMB::oneStepPredict()].
#' @return A numeric vector, `NA` on censored rows.
#' @export
residuals.frmtmb_fit <- function(object, type = c("response", "pearson",
                                                  "deviance", "osa"),
                                 osa_method = NULL, ...) {
  type <- match.arg(type)
  rspec <- uni_resp(object, "residuals()")
  fam <- rspec$family
  if (type == "osa") {
    av0 <- object$frame$aterm_values[[rspec$resp_name]]
    tb <- trunc_bounds(av0, object$frame$n_obs)
    cb <- osa_cens_domain(av0, object$frame$y[[rspec$resp_name]])
    ordinal <- identical(fam$type, "ordinal")
    method <- osa_method %||%
      if (!is.null(tb) || !is.null(cb) || ordinal) "oneStepGeneric"
      else if (identical(fam$family, "gaussian")) "fullGaussian"
      else "oneStepGeneric"
    if (!is.null(cb) && !identical(method, "oneStepGeneric")) {
      # a censored row's contribution is a probability MASS, so its
      # observation drops out of the tape; every method that
      # differentiates the observation hits a singular system there
      stop("residuals(type = \"osa\") on a cens() fit needs ",
           "osa_method = \"oneStepGeneric\"", call. = FALSE)
    }
    args <- list(obj = object$obj, observation.name = ".frm_obs",
                 method = method, trace = FALSE, ...)
    if (method == "oneStepGeneric") {
      args$discrete <- identical(fam$type, "discrete") || ordinal
      if (ordinal) {
        # the taped lpdf is a proper pmf on 1..K once the category is
        # selected arithmetically (see ord_cat_sel)
        if (is.null(args[["discreteSupport"]])) {
          args$discreteSupport <- seq_len(ordinal_ncat(object))
        }
      }
      lo <- -Inf
      hi <- Inf
      if (!is.null(tb)) {
        # The taped density integrates to 1 only over [lb, ub], so the
        # conditional CDF must be built on that domain: over the whole
        # line it sums to 1/P(lb <= Y <= ub) and the residuals come out
        # shrunk. TMB needs one domain for every row.
        lo <- unique(tb$lb)
        hi <- unique(tb$ub)
        if (length(lo) > 1L || length(hi) > 1L) {
          stop("residuals(type = \"osa\") needs trunc() bounds that are ",
               "the same for every observation; got row-varying bounds",
               call. = FALSE)
        }
      }
      if (!is.null(cb)) {
        # An uncensored row is only observed because it fell inside the
        # censoring window, so its PIT is F(y) / P(window) - the same
        # renormalization trunc() needs, on the window's domain.
        lo <- max(lo, cb$lo)
        hi <- min(hi, cb$hi)
        args$subset <- cb$subset
        args$conditional <- cb$conditional
      }
      if (!is.null(tb) || !is.null(cb)) {
        if (is.null(args[["range"]])) args$range <- c(lo, hi)
        if (args$discrete) {
          if (is.null(args[["discreteSupport"]]) && is.finite(hi)) {
            args$discreteSupport <- seq(lo, hi)
          }
        } else if (is.null(args[["splineApprox"]])) {
          # the spline approximation of the transformed density loses
          # about six digits on a bounded domain; exact quadrature
          # reproduces the analytic PIT
          args$splineApprox <- FALSE
        }
      } else if (args$discrete && is.null(args[["range"]]) && !ordinal) {
        args$range <- c(0, Inf)
      }
    }
    osa <- do.call(RTMB::oneStepPredict, args)
    r <- osa$residual
    if (!is.null(cb)) {
      # censored rows carry no residual of their own; they enter only as
      # the conditioning event
      full <- rep(NA_real_, object$frame$n_obs)
      full[cb$subset] <- r
      r <- full
    }
    return(napred(object, r))
  }
  dp <- eval_dpars(object)[[rspec$resp_name]]
  av <- object$frame$aterm_values[[rspec$resp_name]]
  yv <- object$frame$y[[rspec$resp_name]]
  if (type == "deviance") {
    return(napred(object, deviance_residuals(fam, yv, dp, av,
                                             object$frame$n_obs)))
  }
  # on a truncated response the residual is against the truncated mean,
  # which is what the data were actually drawn from
  mu <- response_mean(fam, dp, av)
  r <- yv - mu
  if (type == "pearson") {
    if (is.null(fam$post$var_fn)) {
      stop("Family '", fam$family, "' has no variance function; ",
           "pearson residuals are unavailable", call. = FALSE)
    }
    # the scale stays the untruncated family variance; only the centering
    # is truncation-aware, so pearson residuals on a truncated model are
    # slightly conservative
    v <- fam$post$var_fn(dp, av)
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
      # x = U^-1 z with U'U = Q has covariance Q^-1; for correlated
      # slopes Q is the level-major Kronecker precision Q (x) Sigma^-1
      Qb <- if (bk$dim == 1L) {
        exp(-2 * th[bk$theta_idx[1]]) * bk$aux_Q
      } else {
        S <- covstruct_registry$gr_prec$vcov(th[bk$theta_idx], bk)
        Matrix::kronecker(bk$aux_Q, methods::as(solve(unname(S)),
                                                "generalMatrix"))
      }
      U <- Matrix::chol(Qb)
      b[bk$b_idx] <- as.vector(Matrix::solve(U,
                                             stats::rnorm(length(bk$b_idx))))
      next
    }
    if (bk$covstruct == "car") {
      # the whole field is one draw from its (dense) covariance
      K <- car_cov(th[bk$theta_idx], bk)
      b[bk$b_idx] <- drop(crossprod(chol(K), stats::rnorm(nrow(K))))
      next
    }
    if (bk$covstruct == "spde") {
      U <- Matrix::chol(spde_prec(th[bk$theta_idx], bk))
      b[bk$b_idx] <- as.vector(Matrix::solve(U,
                                             stats::rnorm(bk$n_levels)))
      next
    }
    V <- covstruct_registry[[bk$covstruct]]$vcov(th[bk$theta_idx], bk)
    L <- chol(V)
    U <- matrix(stats::rnorm(bk$n_levels * bk$dim), bk$n_levels) %*% L
    b[bk$b_idx] <- as.vector(t(U))   # level-major
  }
  b
}

# The observation window the censoring mechanism imposes: the censoring
# point is the response value on a censored row, and it must be the
# same for every censored row on a side (type-I censoring). With
# row-varying censoring times the point is unknown on the rows that
# were NOT censored - the data only say it was never reached - so the
# mechanism cannot be applied to their draws at all.
cens_window <- function(av, yobs) {
  cen <- av$cens
  if (any(cen == 2)) {
    stop("simulate(censored = TRUE) is defined for left- and ",
         "right-censored rows; an interval-censored observation is an ",
         "interval, not a value", call. = FALSE)
  }
  point <- function(idx, side) {
    p <- unique(yobs[idx])
    if (length(p) > 1L) {
      stop("simulate(censored = TRUE) needs one ", side,
           "-censoring point shared by every censored row (type-I ",
           "censoring); got ", length(p), " distinct points. With ",
           "row-varying censoring times an uncensored row's censoring ",
           "point is unknown, so the mechanism cannot be applied to ",
           "its draws", call. = FALSE)
    }
    p
  }
  i_r <- which(cen == 1)
  i_l <- which(cen == -1)
  list(lo = if (length(i_l)) point(i_l, "left") else -Inf,
       hi = if (length(i_r)) point(i_r, "right") else Inf)
}

# Apply the censoring mechanism to one draw of the latent response:
# every draw outside the observation window is recorded at the window's
# edge, exactly as the observed data were.
apply_censoring <- function(y, win) {
  pmin(pmax(y, win$lo), win$hi)
}

#' Simulate responses from a frmtmb fit
#'
#' A `trunc()`ed response simulates by rejection within its bounds, so
#' every draw lies in `[lb, ub]` and posterior-predictive checks
#' ([dharma_residuals()], `pp_check()`) see the same support the
#' likelihood was normalized on.
#'
#' @section Censored responses:
#' On a `cens()` fit the default draws the LATENT, uncensored response:
#' the model describes the latent distribution, and censoring is a
#' property of the observation process, not of the response. This
#' matches brms, whose `posterior_predict()` also ignores `cens()` (and
#' whose `pp_check()` therefore drops the censored rows). The draws are
#' then not comparable with the observed values on censored rows, which
#' is why `dharma_residuals()` and `residuals(type = "osa")` refuse or
#' skip them.
#'
#' `censored = TRUE` applies the fitted censoring mechanism to each
#' draw instead, so the draws are directly comparable with the observed
#' data: every draw is recorded at the edge of the observation window
#' it falls outside, capped above by the right-censoring point and
#' below by the left-censoring point. Those points are the response
#' values of the censored rows, and they must be the same for every
#' censored row on a side (type-I censoring): with row-varying
#' censoring times an uncensored row's censoring point is unknown, so
#' the mechanism cannot be applied to its draws and the call is
#' refused. Interval censoring has no single-value representation and
#' is refused too.
#'
#' @param object A `frmtmb_fit`.
#' @param nsim Number of simulated response vectors.
#' @param seed Optional RNG seed. Follows the [stats::simulate()]
#'   contract: the global RNG state is restored afterwards, and the
#'   seed used is attached as the `"seed"` attribute.
#' @param re.form `NULL` (default) conditions on the estimated random
#'   effects; `NA` redraws them from their estimated distribution
#'   (marginal simulation).
#' @param censored Apply the fitted `cens()` mechanism to the draws
#'   (see Censored responses). Ignored without `cens()`.
#' @param ... Unused.
#' @return A data frame with `nsim` columns and a `"seed"` attribute.
#' @export
simulate.frmtmb_fit <- function(object, nsim = 1, seed = NULL,
                                re.form = NULL, censored = FALSE, ...) {
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
  cwin <- NULL
  if (isTRUE(censored)) {
    if (is.null(av$cens)) {
      stop("simulate(censored = TRUE) needs a cens() response",
           call. = FALSE)
    }
    cwin <- cens_window(av, object$frame$y[[rspec$resp_name]])
  }
  for (s in seq_len(nsim)) {
    b_use <- if (marginal && length(object$frame$re_blocks)) {
      draw_b(object)
    } else {
      object$estimates[["b"]]
    }
    dp <- with_cs_offsets(object, rspec, eval_dpars(object,
                                                    b = b_use))
    dp <- dp[[rspec$resp_name]]
    out[[s]] <- if (is.null(mg)) {
      sim_response(fam, dp, av, n, extra = fit_extras(object))
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
    if (!is.null(cwin)) out[[s]] <- apply_censoring(out[[s]], cwin)
  }
  names(out) <- paste0("sim_", seq_len(nsim))
  out <- lapply(out, function(v) sim_restore_type(object, rspec, v))
  out <- sim_as_data_frame(out)
  attr(out, "seed") <- rng_state
  out
}

# Extra (non-dpar) parameters of a fit, in the shape the family's
# simulator expects: ordinal thresholds, category-specific coefficients.
fit_extras <- function(fit) {
  nms <- fit$frame$extra_names %||% character(0)
  if (!length(nms)) return(NULL)
  fit$estimates[nms]
}

# cs(x) contributes an n x (K-1) matrix of threshold-specific offsets.
# The objective builds it on the tape; eval_dpars() has no reason to,
# so add it here for the ordinal simulators that consume it.
with_cs_offsets <- function(fit, rspec, dpv) {
  for (lp in fit$frame$linpreds) {
    if (!length(lp$cs %||% list())) next
    CS <- 0
    for (ct in lp$cs) {
      CS <- CS + outer(ct$vals, fit$estimates[[ct$par]])
    }
    dpv[[lp$resp]][[".cs"]] <- CS
  }
  dpv
}

# simulate() hands draws back in the response's own type: an ordered
# factor for an ordinal fit (the 1..K codes mean nothing without the
# levels) and a matrix for a matrix response. na.exclude fits pad back
# to the original row count, the same contract fitted() and residuals()
# keep. [glmmTMB test-simulate.R; lme4#737]
sim_restore_type <- function(fit, rspec, v) {
  lv <- fit$frame$y_levels[[rspec$resp_name]]
  if (!is.null(lv)) {
    v <- factor(lv[v], levels = lv, ordered = TRUE)
  } else if (is.matrix(v)) {
    yv <- fit$frame$y[[rspec$resp_name]]
    if (is.matrix(yv) && !is.null(colnames(yv))) colnames(v) <- colnames(yv)
  }
  napred(fit, v)
}

# A matrix response needs a data frame whose COLUMNS are matrices (the
# lme4 convention); the default as.data.frame() would flatten each draw
# into one column per category.
sim_as_data_frame <- function(out) {
  if (!is.matrix(out[[1L]])) return(as.data.frame(out))
  df <- data.frame(row.names = seq_len(nrow(out[[1L]])))
  for (nm in names(out)) df[[nm]] <- out[[nm]]
  df
}

# Drop the rows napredict() padded back in, for the internal consumers
# that work in fitted-row space (bootstrap refits, DHARMa, pp_check).
na_unpad <- function(fit, x) {
  na <- fit$frame$na_action
  if (is.null(na) || !inherits(na, "exclude")) return(x)
  idx <- unclass(na)
  if (is.matrix(x) || is.data.frame(x)) x[-idx, , drop = FALSE] else x[-idx]
}
