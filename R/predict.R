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
    Q <- RTMB::sdreport(fit$obj, getJointPrecision = TRUE)$jointPrecision
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
  X <- stats::model.matrix(tt, mfp, contrasts.arg = lp$contrasts)
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

  if (!use_re) {
    # smooth wiggly parts are part of the curve, not group-level effects:
    # they stay in even for population-level predictions
    return(list(X = X, off = off, re_parts = list(), sm_parts = sm_parts))
  }

  re_parts <- list()
  for (bk in fit$frame$re_blocks) {
    if (bk$covstruct == "smooth") next
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
  list(X = X, off = off, re_parts = re_parts, sm_parts = sm_parts)
}

# RE contribution to eta for one linear predictor, given the full b.
re_eta <- function(re_parts, b, n) {
  eta <- numeric(n)
  for (rp in re_parts) {
    bk <- rp$bk
    B <- t(matrix(b[bk$b_idx], nrow = bk$dim))  # levels x D
    cols <- rp$comp$offset + seq_len(rp$comp$dim)
    contrib <- rowSums(rp$mm * B[rp$j, cols, drop = FALSE])
    contrib[is.na(contrib)] <- 0                # new levels: population
    eta <- eta + contrib
  }
  eta
}

# Smooth wiggly contribution to eta for one linear predictor.
sm_eta <- function(sm_parts, b) {
  eta <- 0
  for (sp in sm_parts) {
    eta <- eta + drop(sp$Xr %*% b[sp$bk$b_idx])
  }
  eta
}

# Numeric dpar values at the estimates (optionally with a supplied b),
# for the training data. Nested: out[[resp]][[dpar]].
eval_dpars <- function(fit, b = fit$estimates[["b"]]) {
  est <- fit$estimates
  out <- list()
  for (lp in fit$frame$linpreds) {
    if (!is.null(lp$nl_body)) {
      ev <- c(out[[lp$resp]], lp$data_list)
      eta <- eval(lp$nl_body, ev, lp$nl_env)
      out[[lp$resp]][[lp$dpar]] <- lp$link$linkinv(eta)
      next
    }
    eta <- if (ncol(lp$X)) {
      drop(as.matrix(lp$X %*% est[[lp$par]][lp$idx]))
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

# Sparse RE design (n x length(b)) for the newdata delta method,
# columns at global b positions.
re_design_matrix <- function(re_parts, n, q) {
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  for (rp in re_parts) {
    bk <- rp$bk
    D <- bk$dim
    ok <- which(!is.na(rp$j))
    for (k in seq_len(rp$comp$dim)) {
      ii <- c(ii, ok)
      jj <- c(jj, bk$b_idx[(rp$j[ok] - 1L) * D + rp$comp$offset + k])
      xx <- c(xx, rp$mm[ok, k])
    }
  }
  Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, q))
}

# The single response of a univariate fit, or an informative error.
uni_resp <- function(fit, what) {
  if (length(fit$spec$responses) > 1) {
    stop(what, " is not supported yet for multivariate fits",
         call. = FALSE)
  }
  fit$spec$responses[[1]]
}

#' Predictions from a frmtmb fit
#'
#' @param object A `frmtmb_fit`.
#' @param newdata Optional data frame to predict on. Defaults to the
#'   training data.
#' @param type `"link"` for the linear predictor, `"response"` for the
#'   dpar on its natural scale.
#' @param dpar Which distributional parameter to predict; defaults to the
#'   family's first location parameter (`"mu"` for most families).
#' @param resp For multivariate fits: which response to predict (defaults
#'   to the first).
#' @param re.form `NULL` (default) includes random effects; `NA` or `~0`
#'   gives population-level predictions.
#' @param se.fit If `TRUE`, return a list with elements `fit` and `se.fit`
#'   (delta-method standard errors accounting for fixed-effect and
#'   random-effect uncertainty).
#' @param allow_new_levels Predict unseen grouping-factor levels at the
#'   population level instead of erroring.
#' @param ... Unused.
#' @return A numeric vector, or a list when `se.fit = TRUE`.
#' @export
predict.frmtmb_fit <- function(object, newdata = NULL,
                               type = c("link", "response"),
                               dpar = NULL, resp = NULL, re.form = NULL,
                               se.fit = FALSE,
                               allow_new_levels = FALSE, ...) {
  if (...length()) {
    warning("ignoring unknown arguments to predict(): ",
            paste(...names(), collapse = ", "), call. = FALSE)
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
  sm_blocks <- Filter(function(bk) {
    bk$covstruct == "smooth" &&
      any(vapply(bk$components, function(cp) cp$lp_key == key, TRUE))
  }, object$frame$re_blocks)

  if (is.null(newdata)) {
    X <- lp$X
    off <- lp$offset
    n <- object$frame$n_obs
    eta <- drop(as.matrix(X %*% est[[lp$par]][lp$idx]))
    if (!is.null(lp$Z)) {
      if (use_re) {
        eta <- eta + as.numeric(lp$Z %*% est[["b"]])
      } else if (length(sm_blocks)) {
        # population-level: drop group effects but keep the smooth curve
        for (bk in sm_blocks) {
          eta <- eta + as.numeric(lp$Z[, bk$b_idx, drop = FALSE] %*%
                                    est[["b"]][bk$b_idx])
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
    n <- nrow(X)
    eta <- drop(as.matrix(X %*% est[[lp$par]][lp$idx]))
    if (use_re && length(re_parts)) {
      eta <- eta + re_eta(re_parts, est[["b"]], n)
    }
    if (length(sm_parts)) {
      eta <- eta + sm_eta(sm_parts, est[["b"]])
    }
  }
  if (!is.null(off)) eta <- eta + off

  if (!se.fit) {
    out <- if (type == "response") lp$link$linkinv(eta) else eta
    return(if (is.null(newdata)) napred(object, out) else out)
  }

  # Delta method: var(eta) = A V A' over the estimated coefficients (and b
  # when random effects are included).
  jc <- get_joint_cov(object)
  rn <- jc$names
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
    if (is.null(newdata)) {
      if (use_re && !is.null(lp$Z)) {
        A <- Matrix::cbind2(A, lp$Z)
        coef_pos <- c(coef_pos, b_pos)
      } else if (!use_re && length(sm_blocks) && !is.null(lp$Z)) {
        for (bk in sm_blocks) {
          A <- Matrix::cbind2(A, lp$Z[, bk$b_idx, drop = FALSE])
          coef_pos <- c(coef_pos, b_pos[bk$b_idx])
        }
      }
    } else {
      if (use_re && length(re_parts)) {
        Zn <- re_design_matrix(re_parts, n, length(est[["b"]]))
        A <- Matrix::cbind2(A, Zn)
        coef_pos <- c(coef_pos, b_pos)
      }
      for (sp in sm_parts) {
        A <- Matrix::cbind2(A, sp$Xr)
        coef_pos <- c(coef_pos, b_pos[sp$bk$b_idx])
      }
    }
  }
  V <- jc$V[coef_pos, coef_pos, drop = FALSE]
  A <- as.matrix(A)
  se_eta <- sqrt(pmax(rowSums((A %*% V) * A), 0))

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
  if (!"mu" %in% names(dp)) {
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
    if (bk$covstruct == "gr_cov") {
      # correlation is across levels, not within them
      Sigma <- exp(th[bk$theta_idx])^2 * bk$aux_A
      b[bk$b_idx] <- drop(crossprod(chol(Sigma),
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

#' Simulate responses from a frmtmb fit
#'
#' @param object A `frmtmb_fit`.
#' @param nsim Number of simulated response vectors.
#' @param seed Optional RNG seed.
#' @param re.form `NULL` (default) conditions on the estimated random
#'   effects; `NA` redraws them from their estimated distribution
#'   (marginal simulation).
#' @param ... Unused.
#' @return A data frame with `nsim` columns.
#' @export
simulate.frmtmb_fit <- function(object, nsim = 1, seed = NULL,
                                re.form = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  rspec <- uni_resp(object, "simulate()")
  fam <- rspec$family
  if (is.null(fam$sim)) {
    stop("Family '", fam$family, "' has no simulator yet", call. = FALSE)
  }
  marginal <- !is.null(re.form) && !inherits(re.form, "formula") &&
    is.na(re.form)
  n <- stats::nobs(object)
  out <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    b_use <- if (marginal && length(object$frame$re_blocks)) {
      draw_b(object)
    } else {
      object$estimates[["b"]]
    }
    dp <- eval_dpars(object, b = b_use)[[rspec$resp_name]]
    out[[s]] <- fam$sim(dp, object$frame$aterm_values[[rspec$resp_name]], n)
  }
  names(out) <- paste0("sim_", seq_len(nsim))
  as.data.frame(out)
}
