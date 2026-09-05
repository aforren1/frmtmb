# Prediction, fitted values, residuals, and simulation.
#
# Every linear predictor's Z spans the full b vector, so Z column indices
# equal b indices; |ID|-merged blocks need no special casing on the
# in-sample path. For newdata, each block contribution is rebuilt from the
# block's components (one component per contributing linear predictor).

#' Memoized joint covariance of all estimated parameters (fixed + random),
#' with row/col component names. Cached in `fit$cache`.
#'
#' @noRd
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

#' The three argument checks `predict()` and `frm_lp_basis()` share.
#'
#' They take the same `newdata`, `re.form` and `resp`, and two functions
#' describing one argument two ways is how a user learns that the second
#' one is a different argument. One template each, which is also the
#' property `test-message-uniqueness.R` asserts.
#'
#' @noRd
check_re_form <- function(re.form) {
  if (!is.null(re.form) && !inherits(re.form, "formula") &&
        !(length(re.form) == 1L && is.na(re.form))) {
    stop("`re.form` must be NULL to keep every random effect, NA to drop ",
         "them all, or a one-sided formula naming the ones to keep, not ",
         arg_desc(re.form), call. = FALSE)
  }
  invisible(NULL)
}

#' Whether `re.form` keeps the random effects.
#'
#' `NULL` keeps them, `NA` drops them, and a one-sided formula keeps
#' them unless it is `~0`. Factored out beside the message templates for
#' the same reason those were: `predict()` and `frm_lp_basis()` take the
#' same argument, and the whole point of section 3 of the review is that
#' the two return the same numbers, which they can only do while they
#' read the argument the same way.
#'
#' @noRd
re_form_keeps <- function(re.form) {
  if (is.null(re.form)) return(TRUE)
  if (!inherits(re.form, "formula")) return(FALSE)
  !identical(deparse1(re.form[[2]]), "0")
}

#' @noRd
stop_unknown_response <- function(object, resp) {
  stop("Unknown response: '", resp, "'. Available: ",
       paste(names(object$spec$responses), collapse = ", "), call. = FALSE)
}

#' @noRd
stop_newdata_missing <- function(v) {
  stop("Variable '", v, "' missing from newdata", call. = FALSE)
}

#' Numeric coefficient-space vector for a fitted model (rr factors
#' expanded through the loadings; identity otherwise).
#'
#' @noRd
coef_b <- function(fit, b = fit$estimates[["b"]]) {
  if (is.null(b)) return(b)
  expand_b(fit$frame, b, fit$estimates[["theta"]])
}

#' Numeric `mo()` column values: D times the cumulative simplex at the
#' category codes, evaluated at the current simplex estimates.
#'
#' @noRd
mo_col_values <- function(fit, mi, codes = mi$codes) {
  zeta <- exp(c(0, fit$estimates[[mi$zeta]]))
  zeta <- zeta / sum(zeta)
  cz0 <- c(0, cumsum(zeta))
  mi$D * cz0[codes + 1L]
}

#' Category codes of a `mo()` variable in new data, validated against the
#' fitted range.
#'
#' @noRd
mo_codes <- function(fit, lp, mi, newdata) {
  env <- fit$spec$responses[[lp[["resp"]]]]$formula_env
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

#' Observed-or-latent values of a `mi()` response at the estimates.
#'
#' @noRd
mi_values <- function(fit, vn) {
  xv <- fit$frame[["y"]][[vn]]
  mm_ <- fit$frame[["mi_map"]][[vn]]
  if (!is.null(mm_)) xv[mm_$rows] <- fit$estimates[["miss"]][mm_$idx]
  xv
}

#' Fill the zero placeholder columns of a stored design matrix with the
#' `mo()` and `mi()` values at the current estimates.
#'
#' @noRd
patch_mo_cols <- function(fit, lp, X) {
  for (mi in lp[["mo"]] %||% list()) {
    v <- mo_col_values(fit, mi)
    if (!is.null(mi$mult)) v <- v * mi$mult
    X[, mi$col] <- v
  }
  for (mt in lp[["mi"]] %||% list()) {
    v <- mi_values(fit, mt$var)
    if (!is.null(mt$mult)) v <- v * mt$mult
    X[, mt$col] <- v
  }
  X
}

#' `xlevels` restricted to the variables a terms object actually uses;
#' extra entries make `model.frame` warn.
#'
#' @noRd
xlev_for <- function(xlevels, tt) {
  xlevels[intersect(names(xlevels), all.vars(tt))]
}

#' The grouping factor a smooth basis is indexed by, or `NULL` when the
#' smooth is a population term.
#'
#' A smooth's wiggly part is a random-effect block in the fitted
#' objective, but that is an implementation fact, not a statement about
#' the model: `s(t)` is a population effect that happens to be penalized.
#' What separates the two is whether the basis gives every level of a
#' GROUPING factor its own exchangeable deviation. Three bases do:
#' `bs = "fs"` (one curve per level, class `fs.interaction`, which names
#' the factor in `$fterm`), `bs = "re"` (class `random.effect`, a plain
#' random intercept or slope), and a tensor product with an `re` margin
#' (`t2(t, g, bs = c("cr", "re"))`, the gamm4 spelling of the same
#' thing), which is reached by recursing into `$margin`.
#'
#' `bs = "sz"` is deliberately NOT one of them. It carries the same
#' `$fterm` field, but its level curves are contrasts against a
#' reference level, which is mgcv's spelling for a factor whose levels
#' are fixed effects, so the term stays population-level.
#'
#' The classification is asked of the smooth object itself rather than
#' of the `bs` string, so an alias or a user basis that constructs one
#' of these classes is caught too. `mf` is the fitting model frame,
#' which is what says whether a `random.effect` term is a factor (a
#' numeric term there is a single shrunken coefficient, not a group).
#'
#' @noRd
smooth_group_var <- function(sm, mf = NULL) {
  if (inherits(sm, "sz.interaction")) return(NULL)
  if (inherits(sm, "fs.interaction")) return(sm$fterm)
  if (inherits(sm, "random.effect")) {
    if (is.null(mf)) return(NULL)
    fv <- Filter(function(v) {
      !is.null(mf[[v]]) && (is.factor(mf[[v]]) || is.character(mf[[v]]))
    }, sm$term)
    return(if (length(fv)) fv[[length(fv)]] else NULL)
  }
  for (mg in sm$margin %||% list()) {
    g <- smooth_group_var(mg, mf)
    if (!is.null(g)) return(g)
  }
  NULL
}

#' The `re_blocks` indices of one linear predictor's GROUP-indexed
#' smooths, which is what `re.form = NA` drops.
#'
#' @noRd
smooth_group_block_ids <- function(lp) {
  ids <- integer(0)
  for (si in lp[["smooths"]] %||% list()) {
    if (!is.null(si$group_var)) ids <- c(ids, si$block_ids)
  }
  ids
}

#' The data columns one smooth's `PredictMat()` reads. `$term` holds
#' variable names, but a basis built on an expression keeps the
#' expression, so the names come out of the parse; `$by` is the string
#' `"NA"` when there is no `by` variable.
#'
#' @noRd
smooth_pred_vars <- function(sm) {
  vs <- as.character(sm$term)
  by <- as.character(sm$by %||% "NA")
  if (length(by) == 1L && !identical(by, "NA")) vs <- c(vs, by)
  unique(unlist(lapply(vs, function(v) {
    tryCatch(all.vars(str2lang(v)), error = function(e) v)
  }), use.names = FALSE))
}

#' Refuse a `newdata` one smooth cannot be rebuilt on, before
#' `PredictMat()` reports it as an internal length mismatch.
#'
#' Two faults meet here and they have different fixes. A missing
#' GROUPING column is only needed because the prediction is conditional
#' on the group, so `re.form = NA` is a way out of it; any other missing
#' column is simply absent data. The unseen-level check restates for a
#' factor-smooth term what the ordinary random-effect blocks already
#' promise: a level the fit never saw errors unless it is allowed
#' explicitly. mgcv's own basis matches levels by LABEL and returns a
#' zero row for one it does not know, so allowing them predicts the
#' population curve there, exactly as a new level of `(1 | g)` does.
#'
#' @noRd
smooth_newdata_check <- function(si, newdata, use_re, allow_new_levels) {
  gv <- si$group_var
  miss <- setdiff(smooth_pred_vars(si$sm), names(newdata))
  if (length(miss)) {
    if (!is.null(gv) && gv %in% miss) {
      stop("predict(newdata = ) for the factor-smooth term ", si$label,
           " needs the grouping column `", gv, "`: the term holds one ",
           "curve per level of `", gv, "`, so a prediction conditional ",
           "on it has to say which level each row belongs to. Add the ",
           "column to newdata, or ask for the population curve with ",
           "re.form = NA, which drops the term and needs no level",
           call. = FALSE)
    }
    stop("predict(newdata = ) for the smooth term ", si$label,
         " needs the column(s) ", paste0("`", miss, "`", collapse = ", "),
         ", which newdata does not have", call. = FALSE)
  }
  # fitted levels: fs smooths carry sm$flev; factor bs = "re" smooths
  # do not, so the frame records group_levels for both (older fits
  # without the field fall back to flev and keep their old behavior)
  lev <- si$group_levels %||% si$sm$flev
  if (use_re && !is.null(gv) && !is.null(lev)) {
    new <- setdiff(unique(as.character(newdata[[gv]])),
                   as.character(lev))
    if (length(new) && !allow_new_levels) {
      stop("New levels in the factor-smooth term ", si$label, ": ",
           paste(new, collapse = ", "), ". The term has no curve for ",
           "them. Use allow_new_levels = TRUE to predict them at the ",
           "population level, or re.form = NA for the population curve ",
           "at every row", call. = FALSE)
    }
    if (length(new) && is.null(si$sm$flev)) {
      # an fs basis zero-rows an unknown level; a factor bs = "re"
      # basis has one design column per fitted level and nothing else,
      # so there is no population row to hand back
      stop("allow_new_levels = TRUE cannot predict the new level(s) ",
           paste(new, collapse = ", "), " of the bs = \"re\" smooth ",
           "term ", si$label, ": its design has one column per fitted ",
           "level and no zero row for a new one. Use re.form = NA for ",
           "the population curve, which is what a new level would ",
           "receive anyway", call. = FALSE)
    }
  }
  invisible(NULL)
}

#' Rebuild the design pieces of one linear predictor for new data:
#' dense X (parametric + smooth null-space columns), per-block RE
#' component designs with level indices, and the offset.
#'
#' @noRd
pred_design <- function(fit, lp, newdata, allow_new_levels = FALSE,
                        use_re = TRUE) {
  env <- fit$spec$responses[[lp[["resp"]]]]$formula_env
  tt <- patch_predvars(lp[["terms"]], fit$frame[["predvar_map"]])
  mfp <- stats::model.frame(tt, newdata, na.action = stats::na.pass,
                            xlev = xlev_for(lp[["xlevels"]], tt))
  # sparse_x fits keep newdata designs sparse too; NA rows must stay NA
  # in eta, and sparse.model.matrix silently zeroes NA factor rows, so
  # frames with NAs fall back to the dense builder
  X <- if (isTRUE(fit$frame[["sparse_x"]]) && !anyNA(mfp)) {
    sparse_mm(tt, mfp, contrasts.arg = lp[["contrasts"]])
  } else {
    stats::model.matrix(tt, mfp, contrasts.arg = lp[["contrasts"]])
  }
  # Estimability against a rank-deficient fit. A prediction is a linear
  # functional of beta, so it is identified only when the new design row
  # is orthogonal to every direction the fit could not resolve. Testing
  # the frozen null space (not merely "a dropped column is nonzero")
  # keeps rows that restate a kept column - x2 = 2 * x, or a cell whose
  # aliased indicator is implied by the kept ones - exact.
  nonest <- rep(FALSE, nrow(X))
  if (!is.null(lp[["alias_null"]])) {
    Xa <- as.matrix(X[, rownames(lp[["alias_null"]]), drop = FALSE])
    scl <- max(abs(Xa[is.finite(Xa)]), 1)
    # NA rows already predict NA; keep them out of the estimability vote
    nonest <- rowSums(abs(Xa %*% lp[["alias_null"]]) > 1e-8 * scl,
                      na.rm = TRUE) > 0
  }
  # frozen intercept-drop / rank-deficiency column set from fit time
  X <- X[, lp[["param_colnames"]], drop = FALSE]
  off <- extract_offset(tt, mfp, env)

  # Smooths: rebuild the (wiggly, fixed) split of the basis for newdata.
  # Either way the result has the wiggly columns first (in rand order)
  # and the null-space columns last, matching the fitted X layout.
  sm_parts <- list()
  for (si in lp[["smooths"]] %||% list()) {
    # A smooth indexed by a grouping factor holds that factor's own
    # deviations, so a population-level prediction drops it the way it
    # drops (1 | g). smooth2random() leaves such a term no null space,
    # so dropping it usually removes the term entirely, and with it the
    # need for the grouping column in newdata. A group smooth
    # that DOES have unpenalized columns still has them rebuilt, because
    # they sit in X and the coefficient vector is not reindexed here.
    drop_grp <- !use_re && !is.null(si$group_var)
    if (drop_grp && si$nf == 0L) next
    smooth_newdata_check(si, newdata, use_re, allow_new_levels)
    M <- mgcv::PredictMat(si$sm, newdata)
    if (is.null(si$U)) {
      # t2(): smooth2random() gives no rotation, only pen.ind, so the
      # split is a trans.D scaling plus the frozen column permutation
      # (see smooth_pen_order()). NULL `ord` means that identity did not
      # hold at fit time; refusing beats returning wrong numbers.
      if (is.null(si$ord)) {
        stop("predict(newdata = ) is not supported for the smooth ",
             si$label, ": mgcv reported a random-effect split this ",
             "version cannot invert. In-sample fitted()/predict() work; ",
             "predict at the observed rows instead", call. = FALSE)
      }
      M <- sweep(M, 2, si$D, `*`)[, si$ord, drop = FALSE]
    } else {
      M <- sweep(M %*% si$U, 2, si$D, `*`)
    }
    pos <- 0L
    for (r in seq_along(si$nr)) {
      Xr_new <- M[, pos + seq_len(si$nr[r]), drop = FALSE]
      pos <- pos + si$nr[r]
      if (drop_grp) next
      sm_parts[[length(sm_parts) + 1L]] <- list(
        bk = fit$frame[["re_blocks"]][[si$block_ids[r]]],
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
  for (gi in lp[["gps"]] %||% list()) {
    Xc <- do.call(cbind, lapply(gi$exprs, function(ex) {
      as.numeric(eval(ex, newdata, env))
    }))
    bk <- fit$frame[["re_blocks"]][[gi$block_id]]
    extra_var <- NULL
    if (gi$type == "hsgp") {
      # dmax/center/L are the fitted scaling, so an in-sample newdata
      # row rebuilds its fitted basis row bit for bit
      Xr <- hsgp_basis(sweep(Xc / gi$dmax, 2, gi$center), gi$omega, gi$L)
    } else {
      pos <- gi$positions
      j <- match(pos_rowkey(Xc), pos_rowkey(pos))
      if (anyNA(j)) {
        th <- fit$estimates[["theta"]][bk[["theta_idx"]]]
        K <- unname(covstruct_registry[["gp"]]$vcov(th, bk))
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
  for (mi in lp[["mo"]] %||% list()) {
    v <- mo_col_values(fit, mi, mo_codes(fit, lp, mi, newdata)) *
      nd_mult(mi$mult_expr)
    X <- cbind(X, matrix(v, ncol = 1, dimnames = list(NULL, mi$label)))
  }
  for (mt in lp[["mi"]] %||% list()) {
    v <- newdata[[mt$var]]
    if (is.null(v) || anyNA(v)) {
      stop("mi(", mt$var, "): newdata must supply complete values",
           call. = FALSE)
    }
    X <- cbind(X, matrix(as.numeric(v) * nd_mult(mt$mult_expr),
                         ncol = 1, dimnames = list(NULL, mt$label)))
  }

  if (!use_re) {
    # a POPULATION smooth's wiggly part is part of the curve, not a
    # group-level effect, so it stays in; the group-indexed ones were
    # already left out of sm_parts above
    return(list(X = X, off = off, re_parts = list(), sm_parts = sm_parts,
                nonest = nonest))
  }

  re_parts <- list()
  for (bk in fit$frame[["re_blocks"]]) {
    if (bk[["covstruct"]] %in% c("smooth", "gp", "hsgp")) next
    for (comp in bk[["components"]]) {
      if (comp$lp_key != linpred_key(lp[["resp"]], lp[["dpar"]])) next
      if (!is.null(comp$mm)) {
        # One re_part per MEMBER, its design already scaled by that
        # member's weight. re_eta() and re_design_matrix() both
        # accumulate over re_parts, so the weighted sum falls out with
        # no multi-membership branch of their own, and a member level
        # that is new in newdata drops to the population value the same
        # way a new level of an ordinary factor does.
        re_parts <- c(re_parts,
                      mm_newdata_parts(comp, bk, newdata, env,
                                       lp[["xlevels"]],
                      fit$frame[["predvar_map"]],
                                       allow_new_levels))
        next
      }
      tt2 <- stats::terms(stats::as.formula(call("~", comp$bar[[2]]),
                                            env = env))
      tt2 <- patch_predvars(tt2, fit$frame[["predvar_map"]])
      mf2 <- stats::model.frame(tt2, newdata, na.action = stats::na.pass,
                                xlev = xlev_for(lp[["xlevels"]], tt2))
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
      gv <- if (bk[["covstruct"]] == "spde") spde_node_labels(gvr) else {
        as.character(gvr)
      }
      j <- match(gv, bk[["levels"]])
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

#' Rebuild one multi-membership component's design on newdata, as one
#' `re_parts` entry per member.
#'
#' @noRd
mm_newdata_parts <- function(comp, bk, newdata, env, xlevels,
                             predvar_map, allow_new_levels) {
  mms <- comp$mm
  iw <- mm_index_weights(mms, newdata, env, bk[["levels"]])
  tt2 <- stats::terms(stats::as.formula(call("~", mms$lhs), env = env))
  md <- mm_member_designs(mms, newdata, env, iw$n_members,
                          predvar_map = predvar_map,
                          xlev = xlev_for(xlevels, tt2),
                          use_model_frame = TRUE)
  if (!identical(md$cnms, comp$cnms)) {
    stop("Multi-membership design for `", comp$label, "` does not match ",
         "the fitted model (columns: ",
         paste(md$cnms, collapse = ", "), " vs ",
         paste(comp$cnms, collapse = ", "), ")", call. = FALSE)
  }
  gv <- mm_member_values(mms, newdata, env)
  if (anyNA(iw$J) && !allow_new_levels) {
    new <- unique(unlist(lapply(gv, function(v) {
      setdiff(as.character(v), bk[["levels"]])
    }), use.names = FALSE))
    stop("New levels in multi-membership factor `", mms$label, "`: ",
         paste(new, collapse = ", "),
         ". Use allow_new_levels = TRUE to predict those memberships ",
         "at the population level; the row's remaining members still ",
         "contribute their fitted effects", call. = FALSE)
  }
  lapply(seq_len(iw$n_members), function(k) {
    # new_key names WHICH unseen level this member landed on, because a
    # row whose members carry the SAME unseen label loads one draw of
    # the block, not two independent ones (see extra_var_blocks())
    list(bk = bk, comp = comp,
         mm = md$designs[[k]] * iw$W[, k], j = iw$J[, k],
         new_key = as.character(gv[[k]]))
  })
}

#' RE contribution to eta for one linear predictor, given the full
#' coefficient-space vector (see `coef_b`).
#'
#' @noRd
re_eta <- function(re_parts, cvec, n) {
  eta <- numeric(n)
  for (rp in re_parts) {
    bk <- rp$bk
    B <- t(matrix(cvec[bk[["c_idx"]]], nrow = bk[["dim"]]))  # levels x D
    cols <- rp$comp$offset + seq_len(rp$comp$dim)
    contrib <- rowSums(rp$mm * B[rp$j, cols, drop = FALSE])
    # only unmatched LEVELS predict at the population value; an NA in
    # the RE design data itself must propagate, not silently zero
    contrib[is.na(rp$j)] <- 0
    eta <- eta + contrib
  }
  eta
}

#' Smooth wiggly contribution to eta for one linear predictor.
#'
#' @noRd
sm_eta <- function(sm_parts, cvec) {
  eta <- 0
  for (sp in sm_parts) {
    eta <- eta + drop(sp$Xr %*% cvec[sp$bk[["c_idx"]]])
  }
  eta
}

#' @rdname frmtmb-extension-api
#' @export
eval_dpars <- function(fit, b = fit$estimates[["b"]]) {
  # the chokepoint every prediction, fitted value, residual and
  # simulated draw passes through, so one guard covers them all. A
  # posterior draw has real parameters written in and drops the marker
  # (draws_fit_at), so the draws surface is unaffected.
  require_fitted(fit, "predict() / fitted() / residuals() / simulate()")
  est <- fit$estimates
  if (!is.null(b)) b <- expand_b(fit$frame, b, est[["theta"]])
  out <- list()
  for (lp in fit$frame[["linpreds"]]) {
    if (!is.null(lp[["nl_body"]])) {
      ev <- c(out[[lp[["resp"]]]][c(lp[["nl_pars"]], lp[["nl_dpar_refs"]])],
              lp[["data_list"]], ps_env(lp, est, b))
      eta <- eval(lp[["nl_body"]], ev, ad_overload_env(lp[["nl_env"]],
        lp[["nl_body"]]))
      out[[lp[["resp"]]]][[lp[["dpar"]]]] <- lp[["link"]]$linkinv(eta)
      next
    }
    eta <- if (ncol(lp[["X"]])) {
      drop(as.matrix(patch_mo_cols(fit, lp, lp[["X"]]) %*%
                       est[[lp[["par"]]]][lp[["idx"]]]))
    } else {
      numeric(fit$frame[["n_obs"]])
    }
    # `b = NULL` is documented to drop the random-effect contribution,
    # and dropping it means not forming the product at all: a sparse
    # matrix times NULL sends Matrix's S4 dispatch into a recursion that
    # ends in "evaluation nested too deeply", on every model with a
    # block rather than only on unusual ones
    if (!is.null(lp[["Z"]]) && !is.null(b)) {
      eta <- eta + as.numeric(lp[["Z"]] %*% b)
    }
    if (!is.null(lp[["offset"]])) eta <- eta + lp[["offset"]]
    out[[lp[["resp"]]]][[lp[["dpar"]]]] <- lp[["link"]]$linkinv(eta)
  }
  out
}

#' Sparse RE design (`n x n_c`) for the newdata delta method, columns at
#' global coefficient-space positions.
#'
#' @noRd
re_design_matrix <- function(re_parts, n, q) {
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  for (rp in re_parts) {
    bk <- rp$bk
    D <- bk[["dim"]]
    ok <- which(!is.na(rp$j))
    for (k in seq_len(rp$comp$dim)) {
      ii <- c(ii, ok)
      jj <- c(jj, bk[["c_idx"]][(rp$j[ok] - 1L) * D + rp$comp$offset + k])
      xx <- c(xx, rp$mm[ok, k])
    }
  }
  Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, q))
}

#' Jacobians of the coefficient-space expansion for rr fits: `d cvec/d b`
#' (sparse; identity except the rr blocks' loadings) and `d cvec/d theta`
#' for the rr loading parameters (finite differences on `expand_b`).
#'
#' @noRd
rr_jacobians <- function(fit) {
  frame <- fit$frame
  est <- fit$estimates
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  th_cols <- list()
  for (bk in frame[["re_blocks"]]) {
    if (bk[["covstruct"]] == "rr") {
      L <- rr_loadings(est[["theta"]][bk[["theta_idx"]]], bk[["dim"]],
        bk[["rank"]])
      for (l in seq_len(bk[["n_levels"]])) {
        rows <- bk[["c_idx"]][(l - 1L) * bk[["dim"]] + seq_len(bk[["dim"]])]
        cols <- bk[["b_idx"]][(l - 1L) * bk[["rank"]] + seq_len(bk[["rank"]])]
        ii <- c(ii, rep(rows, times = bk[["rank"]]))
        jj <- c(jj, rep(cols, each = bk[["dim"]]))
        xx <- c(xx, as.vector(L))
      }
      for (j in bk[["theta_idx"]]) {
        h <- 1e-6 * max(1, abs(est[["theta"]][j]))
        tp <- est[["theta"]]; tp[j] <- tp[j] + h
        tn <- est[["theta"]]; tn[j] <- tn[j] - h
        dvec <- (expand_b(frame, est[["b"]], tp) -
                   expand_b(frame, est[["b"]], tn)) / (2 * h)
        th_cols[[length(th_cols) + 1L]] <- list(j = j, dvec = dvec)
      }
    } else {
      ii <- c(ii, bk[["c_idx"]])
      jj <- c(jj, bk[["b_idx"]])
      xx <- c(xx, rep(1, length(bk[["b_idx"]])))
    }
  }
  list(Jb = Matrix::sparseMatrix(i = ii, j = jj, x = xx,
                                 dims = c(frame[["n_c"]],
                                          length(est[["b"]]))),
       th_cols = th_cols)
}

#' @rdname frmtmb-extension-api
#' @export
single_response <- function(fit, what) {
  if (length(fit$spec$responses) > 1) {
    stop(what, " is not supported yet for multivariate fits",
         call. = FALSE)
  }
  fit$spec$responses[[1]]
}

#' Whether the family's response mean is just the mu dpar. All the
#' built-in identity-mean families spell their `mean_fn` as `dpars$mu`, so
#' the structural check is exact for them and conservative (general
#' path) for custom families.
#'
#' @noRd
mean_is_mu <- function(fam) {
  # a family that declares no mean has no mean: reading its absence as
  # "the mean is mu" reported a race model's drift as its fitted value
  !is.null(fam[["post"]]$mean_fn) &&
    identical(body(fam[["post"]]$mean_fn), quote(dpars[["mu"]]))
}

#' Whether a response carries `trunc()` bounds, from the spec (the stored
#' expressions) rather than from evaluated values.
#'
#' @noRd
has_trunc <- function(rspec) {
  any(c("trunc_lb", "trunc_ub") %in% names(rspec$aterms))
}

# Addition-term values (trials, se, trunc bounds, ...) re-evaluated on
# new data, for the expected-response prediction path. Terms the family
# mean cannot use (censoring, structural flags) are skipped; trials, se
# and the truncation bounds must evaluate because omitting them silently
# changes the mean.
#' An addition term as the user wrote it, for error messages.
#'
#' @noRd
aterm_label <- function(nm, ex) {
  switch(nm,
    trunc_lb = paste0("trunc(lb = ", deparse1(ex), ")"),
    trunc_ub = paste0("trunc(ub = ", deparse1(ex), ")"),
    # vint(a, b) is stored one argument per aterm as vint1, vint2, ...
    paste0(sub("[0-9]+$", "", nm), "(", deparse1(ex), ")")
  )
}

#' A custom family's `lpdf`/`mean_fn` reads its `vint()`/`vreal()` payload
#' out of aterms, so an omitted one is not a missing covariate but a
#' missing argument: the family returns a zero-length prediction.
#'
#' A term registered by another package is the same thing under another
#' name, so it is required here for the same reason.
#'
#' @noRd
is_custom_data_aterm <- function(nm) {
  grepl("^v(int|real)[0-9]+$", nm) || !is.null(registered_aterm_of(nm))
}

#' Re-evaluate a response's addition terms on new data, for the
#' expected-response prediction path. A term the family mean needs must
#' evaluate, because dropping it silently changes the mean; any other
#' term is dropped with a warning.
#'
#' @noRd
aterms_for_newdata <- function(rspec, newdata) {
  skip <- c("cens", "cens_y2", "se_sigma", "mi", "mi_sd", "weights")
  need <- c("trials", "se", "trunc_lb", "trunc_ub")
  nd_n <- nrow(newdata)
  av <- list()
  for (nm in setdiff(names(rspec$aterms), skip)) {
    ex <- rspec$aterms[[nm]]
    # the same coercion the frame applied, or newdata's factor would
    # reach the density as level codes where training data reached it as
    # whatever the contributing package meant
    reg_at <- registered_aterm_of(nm)
    v <- tryCatch(
      as.numeric(if (is.null(reg_at)) {
        eval(ex, newdata, rspec$formula_env)
      } else {
        reg_at$coerce(eval(ex, newdata, rspec$formula_env))
      }),
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
  if (!is.null(rspec$aterms[["se_sigma"]])) {
    av[["se_sigma"]] <- rspec$aterms[["se_sigma"]]
  }
  av
}

#' Every dpar of one response on the scale the DENSITY consumes, which
#' is the link inverse of its linear predictor.
#'
#' Not the same thing as `predict(type = "response", dpar = )` any more:
#' a family may declare a REPORTING scale for a dpar (a mixture's
#' mixing weights are a softmax over the component predictors, and the
#' predictor itself is not a probability), and feeding a reported value
#' back to `lpdf`, `sim` or `mean_fn` would apply the transform twice.
#' Every consumer that hands dpar values to the family goes through
#' here; only the user-facing `predict()` surface reports.
#'
#' @noRd
dpars_natural <- function(fit, rspec, newdata, re.form,
                          allow_new_levels = FALSE) {
  rn <- rspec$resp_name
  dp <- list()
  for (dnm in names(rspec$dpars)) {
    lp <- fit$frame[["linpreds"]][[linpred_key(rn, dnm)]]
    eta <- predict(fit, newdata = newdata, dpar = dnm, resp = rn,
                   re.form = re.form, type = "link",
                   allow_new_levels = allow_new_levels)
    dp[[dnm]] <- as.vector(lp[["link"]]$linkinv(eta))
  }
  dp
}

#' `y | se(s)` without `sigma = TRUE`: the residual standard deviation
#' beyond the known `s` is not in the density at all - `resid_sd()`
#' returns `s` alone and the dpar is mapped out at the link-scale zero.
#' `sigma()` has always reported that as 0; `predict(dpar = "sigma")`
#' reported the log link's inverse of the mapped-out coefficient, 1,
#' which reads as an estimate of a parameter the model does not have.
#'
#' @noRd
se_unused_sigma <- function(rspec, dpar) {
  identical(dpar, "sigma") && !is.null(rspec$aterms[["se"]]) &&
    !isTRUE(rspec$aterms[["se_sigma"]])
}

#' A family's reporting transform for one dpar, or NULL when the dpar
#' reports on its own natural scale (which is every dpar of every
#' built-in family except a mixture's `theta`).
#'
#' @noRd
dpar_report_hook <- function(fam, dpar, rspec = NULL) {
  if (!is.null(rspec) && se_unused_sigma(rspec, dpar)) {
    zero <- function(dpars, dnm) 0 * dpars[[dnm]]
    return(list(dpars = dpar, value = zero, deriv = zero))
  }
  h <- fam[["post"]][["dpar_response"]]
  if (is.null(h) || !dpar %in% h[["dpars"]]) return(NULL)
  h
}

#' Expected response over all dpars: the family mean at predicted dpar
#' values (`fitted()`'s convention, extended to newdata and `re.form`).
#'
#' @noRd
predict_mean_response <- function(fit, rspec, newdata, re.form,
                                  allow_new_levels) {
  fam <- rspec$family
  rn <- rspec$resp_name
  if (is.null(newdata) && is.null(re.form)) {
    # exactly fitted(): dpars at the estimates, conditional on the modes
    dp <- eval_dpars(fit)[[rn]]
    out <- response_mean(fam, dp, fit$frame[["aterm_values"]][[rn]])
    return(napred(fit, out))
  }
  dp <- dpars_natural(fit, rspec, newdata, re.form, allow_new_levels)
  av <- if (is.null(newdata)) {
    # in-sample dpar predictions come back napredict-ed; pad the
    # per-observation aterm values the same way (a no-op under na.omit)
    lapply(fit$frame[["aterm_values"]][[rn]], function(v) {
      if (is.numeric(v) && length(v) == fit$frame[["n_obs"]]) {
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
#' @section Ordinal responses:
#' `cumulative()`, `sratio()`, `cratio()` and `acat()` have no mean on
#' the response scale, so `type = "response"` (and its alias
#' `type = "conditional"`) returns an `n x K` matrix of category
#' probabilities instead of a vector - the brms `fitted()` convention -
#' with the response's own factor levels as column names. The rows sum
#' to one. `cs()` category-specific terms are honored: they enter each
#' threshold separately and are re-evaluated on `newdata`.
#'
#' [fitted()] returns the same matrix, so the usual
#' `predict(type = "response") == fitted()` identity holds here too.
#'
#' `type = "link"` (the default) and `dpar = "mu"` still give the latent
#' linear predictor, which is where the fixed-effect coefficients live
#' and where `se.fit` is available. `se.fit` on the response scale is
#' refused: the prediction is a K-vector per row, not one number.
#' `emmeans` and `insight::get_predicted()` stay on that latent scale,
#' which is the `mode = "latent"` convention for `clm`-like models.
#'
#' `type = "conditional"` is glmmTMB's name for the conditional MEAN, so
#' it gives the category probabilities here too rather than the linear
#' predictor: an ordinal response has no mean, and answering a question
#' about a mean with a latent predictor is the confusion this section
#' exists to remove. Ask for the predictor by name (`type = "link"`, or
#' `dpar = "mu"`) when that is what you want.
#' @param resp For multivariate fits: which response to predict (defaults
#'   to the first).
#' @param re.form `NULL` (default) includes random effects; `NA` or `~0`
#'   gives population-level predictions. See
#'   *What `re.form = NA` drops* for what that means when the model has
#'   smooths.
#' @section What `re.form = NA` drops:
#' `re.form = NA` (equivalently `~0`) asks for the POPULATION-level
#' prediction. Every `(x | g)` block is dropped, and so is any smooth
#' whose basis gives each level of a grouping factor its own curve.
#' Everything else stays.
#'
#' Dropped:
#'
#' * `(1 | g)`, `(x | g)`, and the structured spellings of them
#'   (`gr()`, `cs()`, `ar()`, `mm()`, `car()`, `spde()`, ...).
#' * `s(t, g, bs = "fs")`, the factor-smooth interaction: one curve per
#'   level of `g`, so the curves ARE the group deviations.
#' * `s(g, bs = "re")` and `s(x, g, bs = "re")`, which are a random
#'   intercept and a random slope written as a smooth.
#' * `t2(t, g, bs = c("cr", "re"))` and any other tensor product with an
#'   `re` margin, which is the same random smooth in a different
#'   spelling.
#'
#' Kept:
#'
#' * `s(t)`, `s(t, by = x)`, `te()`, `t2()`, and every other population
#'   smooth. A smooth's wiggly part is stored as a random-effect block
#'   because that is how a penalty is written as a mixed model, but the
#'   term is a population effect and the population prediction is the
#'   fitted curve, not the null-space line through it.
#' * `gp()` and `hsgp()` terms.
#'
#' The test is what the basis MEANS, not the `bs` string: `bs = "sz"`
#' names a factor the way `bs = "fs"` does, but writes the level curves
#' as contrasts against a reference level, which is mgcv's spelling for
#' a factor whose levels are fixed effects, so it would count as
#' population-level. (`sz` has no random-effect representation, so it is
#' not fittable here at all; the classification is stated for
#' completeness.)
#'
#' The result is `mgcv::predict.gam(exclude = )` on the factor-smooth
#' term, and `tests/testthat/test-smooth-population.R` asserts the two
#' agree to 1e-6 on a shared fit.
#'
#' This deliberately follows mgcv rather than brms: brms stores every
#' smooth's wiggly part as population parameters, so its
#' `re_formula = NA` KEEPS factor-smooth curves. A ported brms call
#' with a `bs = "fs"` term therefore returns different numbers here,
#' on purpose: the retained per-level curve is not a population
#' quantity, and mgcv, the authority frmtmb's smooth estimation
#' already follows, drops it too.
#'
#' A dropped factor-smooth term needs nothing from `newdata`, so the
#' grouping column may be left out entirely when `re.form = NA`. It is
#' required for a conditional prediction, and its absence is reported by
#' name rather than by an mgcv internal message.
#'
#' `conditional_effects()` draws its curves at the population level, so
#' it follows this rule too: on a model with a factor-smooth term the
#' displayed curve is the population smooth, and the grouping factor is
#' not offered as an effect to plot.
#' @param se.fit If `TRUE`, return a list with elements `fit` and `se.fit`
#'   (delta-method standard errors accounting for fixed-effect and
#'   random-effect uncertainty). Exact `gp()` terms predict unseen
#'   positions by kriging: the conditional mean at the fitted kernel,
#'   with the GP conditional variance added to the standard errors.
#' @section Standard errors of the expected response:
#' For a family whose mean is the `mu` dpar, `se.fit` on
#' A dpar whose response scale is not its own link inverse, such as a
#' `mixture()` mixing weight reporting the softmax, takes the delta
#' method through that transform with respect to its OWN predictor. That is exact for a two-component mixture and
#' conservative for three or more; see `?mixture`.
#'
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
#'   population level instead of erroring. A factor-smooth term
#'   (`bs = "fs"`) follows the same rule: a level it never saw
#'   contributes nothing, which leaves the population curve.
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
#' @return A numeric vector, or a list when `se.fit = TRUE`. For an
#'   ordinal family with `type = "response"`, an `n x K` matrix of
#'   category probabilities.
#'
#' @srrstats {G2.3,G2.3a} `type` is a univariate character parameter and
#'   is restricted with `match.arg()` to the documented set, so an
#'   unexpected value errors and names the permitted ones. `dpar` is
#'   checked against the family's own parameter names, and an unknown
#'   `resp` errors with the available responses listed.
#' @srrstats {G3.0} A computed floating-point quantity is never compared
#'   for equality. The estimability test for a rank-deficient design is
#'   a relative tolerance of `1e-8` against the null space of the fitted
#'   design, the same test [stats::predict.lm()] uses, and the
#'   documentation states the tolerance and its consequence for
#'   near-aliased designs. Exact equality is used only where the
#'   compared value is a data category and not a computed quantity:
#'   `y == 0` and `y == 1` are the structural-zero and one-inflation
#'   indicators of the zero-inflated, hurdle and zero-one-inflated
#'   families, where the point mass is defined at that exact value, and
#'   guards against a degenerate input, such as a zero standard
#'   deviation or a zero row-sum of weights, reject on the boundary
#'   rather than accept there.
#' @srrstats {RE1.3} Output structures retain the relevant aspects of the
#'   input. Predictions, `fitted()` values, and residuals carry the row
#'   names of the data they were computed from; `vcov()`, `confint()`,
#'   and `fixef()` share one coefficient naming scheme;
#'   `tests/testthat/test-methods-audit.R` asserts that
#'   `vcov(full = TRUE)` carries exactly the row names of `confint()`
#'   and that `vcov()` is its coefficient sub-block. The stored model
#'   frame keeps the input row names.
#' @srrstats {RE4.9} Modelled values of the response are returned by
#'   `fitted()`, and by `predict(type = "response")`, which is asserted
#'   to equal `fitted()` on the training data wherever both are defined:
#'   the fuzz harness checks the invariant across the family grid, and
#'   `tests/testthat/test-methods-audit.R` checks it on the cases where
#'   the two could plausibly diverge (zero-inflated, `trials()`
#'   binomial). A family with no mean on the response scale, such as
#'   `cox()`, refuses both alike. On an ordinal family the modelled
#'   response is a category
#'   distribution rather than a mean, so both return the same `n x K`
#'   matrix of category probabilities (the brms convention), named by the
#'   response's own levels; the latent linear predictor stays reachable
#'   as `predict(type = "link")`.
#' @srrstats {RE4.14} Uncertainty is available away from the observed
#'   data. `se.fit = TRUE` returns delta-method standard errors that
#'   include fixed-effect and random-effect uncertainty; unseen grouping
#'   levels add their block's marginal variance, and exact `gp()` terms
#'   add the Gaussian-process conditional (kriging) variance, so the
#'   reported error grows with distance from the observed positions.
#' @srrstats {RE4.16} New groups can be submitted to `predict()`. Levels
#'   of a grouping factor that were not in the training data error by
#'   default, naming the offending levels, and are predicted at the
#'   population level under `allow_new_levels = TRUE` (the lme4 spelling
#'   `allow.new.levels` is accepted as well).
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
#' fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
#'
#' # the link scale by default; "response" is what fitted() returns
#' head(predict(fit))
#' max(abs(predict(fit, type = "response") - fitted(fit)))
#'
#' # re.form = NA drops the random effects: the population prediction
#' nd <- data.frame(x = c(-1, 0, 1), g = factor(1, levels = levels(dd$g)))
#' predict(fit, newdata = nd, re.form = NA, type = "response")
#'
#' # delta-method standard errors, on whichever scale was asked for
#' p <- predict(fit, newdata = nd, se.fit = TRUE)
#' cbind(fit = p$fit, se = p$se.fit)
#'
#' # a level the fit never saw errors unless it is allowed explicitly,
#' # in which case it is predicted at the population level
#' nd_new <- data.frame(x = 0, g = factor("new"))
#' try(predict(fit, newdata = nd_new))
#' predict(fit, newdata = nd_new, allow_new_levels = TRUE)
#'
#' # a distributional parameter instead of the mean
#' fit2 <- frm(bf(y ~ x, sigma ~ x) + gaussian(), data = dd)
#' head(predict(fit2, dpar = "sigma", type = "response"))
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
  require_fitted(object, "predict()")
  # re.form is read three lines down as "NULL, a formula, or anything
  # else drops the random effects", so a typo used to return the
  # POPULATION prediction and say nothing. newdata that is not
  # rectangular reached the design builder and failed there on
  # "non-conformable arrays", which names neither the argument nor the
  # fault.
  check_flag(se.fit, "se.fit")
  check_flag(allow_new_levels, "allow_new_levels")
  if (!is.null(newdata) && !is.data.frame(newdata)) {
    stop("`newdata` must be a data frame, or NULL to predict on the ",
         "training data, not ", arg_desc(newdata), call. = FALSE)
  }
  check_re_form(re.form)
  type <- match.arg(type)
  use_re <- re_form_keeps(re.form)

  resp <- resp %||% names(object$spec$responses)[1]
  rspec <- object$spec$responses[[resp]]
  if (is.null(rspec)) {
    stop_unknown_response(object, resp)
  }
  # An ordinal response has no mean on the response scale: what
  # "response" means there is the category distribution, one row of K
  # probabilities per observation (the brms fitted()/epred convention).
  # The mu predictor is still reachable as type = "link" or dpar = "mu".
  if (identical(rspec$family[["type"]], "ordinal") && is.null(dpar) &&
      type %in% c("response", "conditional")) {
    if (se.fit) {
      stop("se.fit is not supported on the response scale for an ",
           "ordinal family: the prediction is a K-vector of category ",
           "probabilities per row, not one number, and the thresholds ",
           "enter every one of them. Use type = \"link\" for the ",
           "standard error of the latent predictor", call. = FALSE)
    }
    return(predict_ordinal(object, rspec, newdata, use_re,
                           allow_new_levels))
  }
  # A nominal response has no mean either: "response" is the K-vector of
  # category probabilities, the same convention the ordinal families and
  # brms both follow.
  if (identical(rspec$family[["type"]], "categorical") && is.null(dpar) &&
      type %in% c("response", "conditional")) {
    if (se.fit) {
      stop("se.fit is not supported on the response scale for a ",
           "categorical family: the prediction is a K-vector of category ",
           "probabilities per row, not one number. Use type = \"link\" ",
           "with dpar = for the standard error of one category's latent ",
           "predictor", call. = FALSE)
    }
    return(predict_categorical(object, rspec, newdata, use_re,
                               allow_new_levels))
  }
  # A structured family's expected value may condition on the WHOLE
  # observed response rather than on the row: an HMM's is the
  # occupancy-weighted mean, sum_k P(S_t = k | y) mu_k(x_t). No per-row
  # family mean can produce that, and without this branch `type =
  # "response"` falls through to the first location dpar and reports
  # state 1's mean at every row, silently. A family whose per-row mean
  # IS rowwise (a group-level mixture) leaves `fitted_mean` empty and
  # falls through here on purpose.
  st <- fam_structure(rspec$family)
  if (!is.null(st) && is.null(dpar) &&
      type %in% c("response", "conditional")) {
    fam_ <- rspec$family
    if (se.fit) {
      structure_gate(st, "se_fit_response",
                     structure_generic(fam_, "se.fit on the response scale"))
    }
    if (!is.null(newdata)) {
      structure_gate(st, "newdata_response",
                     structure_generic(
                       fam_, "predict(newdata =, type = \"response\")"))
    }
    if (!is.null(re.form)) {
      structure_gate(st, "re_form",
                     structure_generic(
                       fam_, "re.form = on the response scale"))
    }
    fm <- st[["fitted_mean"]]
    if (!is.null(fm)) {
      return(napred(object,
                    fm(object, frame_block_of(object$frame,
                                              rspec$resp_name))))
    }
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
           " parameter; family '", rspec$family[["family"]], "' has none",
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
  lp <- object$frame[["linpreds"]][[key]]
  if (is.null(lp)) {
    stop("Unknown dpar: '", dpar, "' for response '", resp,
         "'. Available: ", paste(names(rspec$dpars), collapse = ", "),
         call. = FALSE)
  }

  if (!is.null(lp[["nl_body"]])) {
    if (se.fit) {
      stop("se.fit is not supported for the nonlinear predictor yet; ",
           "request the nonlinear parameters (dpar = '",
           rspec$nlpars[1], "', ...) instead", call. = FALSE)
    }
    # only the parameters this body names, each on its own linear-
    # predictor scale. A parameter that carries a body of its own comes
    # back through this same branch, so a chain of nlf() formulas
    # unwinds by recursion in dependency order.
    vals <- list()
    for (np in lp[["nl_pars"]]) {
      vals[[np]] <- predict(object, newdata = newdata, dpar = np,
                            resp = resp, re.form = re.form,
                            allow_new_levels = allow_new_levels)
    }
    # a reference to another dpar reads its VALUE, so it comes back
    # through that parameter's link inverse - the NATURAL scale, which
    # is what the body computes with, not a reporting scale
    for (dr in lp[["nl_dpar_refs"]] %||% character(0)) {
      lpr <- object$frame[["linpreds"]][[linpred_key(resp, dr)]]
      vals[[dr]] <- lpr[["link"]]$linkinv(
        predict(object, newdata = newdata, dpar = dr,
                type = "link", resp = resp, re.form = re.form,
                allow_new_levels = allow_new_levels))
    }
    dl <- if (is.null(newdata)) {
      lp[["data_list"]]
    } else {
      lapply(stats::setNames(names(lp[["data_list"]]),
             names(lp[["data_list"]])),
             function(v) {
               if (is.null(newdata[[v]])) {
                 stop_newdata_missing(v)
               }
               newdata[[v]]
             })
    }
    # `check` only on the newdata path: in sample the fit-end report has
    # already said whatever there is to say about the knot span, and
    # saying it again on every fitted() call would be noise.
    eta <- eval(lp[["nl_body"]],
                c(vals, dl, ps_env(lp, object$estimates, coef_b(object),
                                   check = !is.null(newdata))),
                ad_overload_env(lp[["nl_env"]], lp[["nl_body"]]))
    out <- if (type == "response") lp[["link"]]$linkinv(eta) else eta
    return(if (is.null(newdata)) napred(object, out) else out)
  }

  ed <- lp_eta_design(object, lp, newdata, use_re, allow_new_levels)
  eta <- ed[["eta"]]
  n <- ed[["n"]]
  # a dpar whose RESPONSE scale is not its own link inverse: a mixture's
  # mixing weight is the softmax over the component predictors, so the
  # identity-link predictor it used to report was not a probability and
  # was not even bounded by one. The transform is the family's, and the
  # density never sees it (dpars_natural()).
  hook <- if (type == "response") {
    dpar_report_hook(rspec$family, dpar, rspec)
  }
  # dpars_natural() reads through predict(), so an in-sample value comes
  # back ALREADY padded for the rows na.action dropped, while eta and
  # its standard error are on the kept rows. The padding comes off here
  # and goes back on once at the end, rather than twice.
  hook_val <- function(what) {
    nat <- dpars_natural(object, rspec, newdata, re.form,
                         allow_new_levels)
    v <- hook[[what]](nat, dpar)
    if (is.null(newdata) && length(v) != n) {
      v <- v[!is.na(napred(object, rep(1, n)))]
    }
    v
  }

  if (!se.fit) {
    out <- if (!is.null(hook)) {
      hook_val("value")
    } else if (type == "response") {
      lp[["link"]]$linkinv(eta)
    } else {
      eta
    }
    return(if (is.null(newdata)) napred(object, out) else out)
  }

  # The delta method is written as a two-line consumer of
  # frm_lp_basis() on purpose: if predict() could not be written in
  # terms of the exported seam, the seam would be the wrong shape.
  has_rr <- isTRUE(object$frame[["has_rr"]])
  rrj <- if (has_rr) rr_jacobians(object)
  jc <- get_joint_cov(object)
  da <- lp_delta_A(object, lp, ed, newdata, use_re, jc, has_rr, rrj)
  lb <- lp_basis_out(object, jc, ed[["eta"]], as.matrix(da$A), da$coef_pos,
                     lp_extra_var_vec(object, ed, use_re), ed[["nonest"]])
  var_eta <- pmax(rowSums((lb$A %*% lb$V) * lb$A), 0) + lb$extra_var
  se_eta <- sqrt(var_eta)
  # the kept columns still have a finite variance, but it is not the
  # standard error of anything the fit estimates
  if (any(ed[["nonest"]])) se_eta[ed[["nonest"]]] <- NA_real_

  out <- if (!is.null(hook)) {
    # the delta method through the reporting transform, with respect to
    # this dpar's own predictor: the same one-predictor rule the link
    # inverse gets below
    list(fit = hook_val("value"),
         se.fit = abs(hook_val("deriv")) * se_eta)
  } else if (type == "response") {
    list(fit = lp[["link"]]$linkinv(eta),
         se.fit = abs(lp[["link"]]$mu_eta(eta)) * se_eta)
  } else {
    list(fit = eta, se.fit = se_eta)
  }
  if (is.null(newdata)) {
    out$fit <- napred(object, out$fit)
    out$se.fit <- napred(object, out$se.fit)
  }
  out
}

#' eta and the design pieces of one linear predictor, in sample or on
#' newdata. Shared by `predict()` and by the joint delta method for the
#' expected response, so both see exactly the same eta.
#'
#' @noRd
lp_eta_design <- function(object, lp, newdata, use_re, allow_new_levels) {
  est <- object$estimates
  key <- linpred_key(lp[["resp"]], lp[["dpar"]])
  re_parts <- list()
  sm_parts <- list()
  nonest <- FALSE
  sm_ids <- which(vapply(object$frame[["re_blocks"]], function(bk) {
    bk[["covstruct"]] %in% c("smooth", "gp", "hsgp") &&
      any(vapply(bk[["components"]], function(cp) cp$lp_key == key, TRUE))
  }, TRUE))
  # the blocks a population-level prediction keeps: gp()/hsgp() curves
  # and the POPULATION smooths, but not a smooth whose basis is indexed
  # by a grouping factor (see smooth_group_var())
  if (!use_re) sm_ids <- setdiff(sm_ids, smooth_group_block_ids(lp))
  sm_blocks <- object$frame[["re_blocks"]][sm_ids]

  if (is.null(newdata)) {
    X <- patch_mo_cols(object, lp, lp[["X"]])
    off <- lp[["offset"]]
    n <- object$frame[["n_obs"]]
    eta <- drop(as.matrix(X %*% est[[lp[["par"]]]][lp[["idx"]]]))
    if (!is.null(lp[["Z"]])) {
      cvec <- coef_b(object)
      if (use_re) {
        eta <- eta + as.numeric(lp[["Z"]] %*% cvec)
      } else if (length(sm_blocks)) {
        # population-level: drop group effects but keep the smooth curve
        for (bk in sm_blocks) {
          eta <- eta + as.numeric(lp[["Z"]][, bk[["c_idx"]], drop = FALSE] %*%
                                    cvec[bk[["c_idx"]]])
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
    eta <- drop(as.matrix(X %*% est[[lp[["par"]]]][lp[["idx"]]]))
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
            paste(lp[["dropped_colnames"]], collapse = ", "),
            "; predicting NA there", call. = FALSE)
    eta[nonest] <- NA_real_
  }
  list(eta = eta, X = X, off = off, re_parts = re_parts,
       sm_parts = sm_parts, sm_blocks = sm_blocks, nonest = nonest, n = n)
}

#' A standard error that leaves out a variance component is only honest
#' if it says so, so this is a warning rather than a note in the docs.
#'
#' @noRd
warn_modes_conditional_se <- function() {
  warning("The fitted objective marginalizes the random effects ",
          "(quadrature = TRUE), so its covariance carries no ",
          "random-effect block: se.fit is conditional on the ",
          "conditional modes and omits random-effect uncertainty. ",
          "Refit with quadrature = FALSE for the full delta method",
          call. = FALSE)
  invisible(NULL)
}

#' Delta method: `var(eta) = A V A'` over the estimated coefficients (and b
#' when random effects are included). For rr fits the Z matrices span the
#' coefficient space, so the b columns go through the Jacobian of the
#' loadings expansion, and the rr loading parameters (theta) contribute
#' their own columns. Returns A and the positions of its columns in the
#' joint covariance, so several linear predictors can be combined.
#'
#' @noRd
lp_delta_A <- function(object, lp, ed, newdata, use_re, jc, has_rr, rrj) {
  est <- object$estimates
  rn <- jc$names
  X <- ed[["X"]]
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
  if (lp[["par"]] == "beta") {
    coef_pos <- which(rn == "beta")[lp[["idx"]]]
    A <- X
  } else {
    tpl_len <- length(object$frame[["par_template"]][["betad"]])
    est_rank <- match(lp[["idx"]],
                      setdiff(seq_len(tpl_len),
                              object$frame[["betad_fixed_idx"]]))
    keep <- !is.na(est_rank)
    coef_pos <- which(rn == "betad")[est_rank[keep]]
    A <- X[, keep, drop = FALSE]
  }
  if (length(object$frame[["re_blocks"]])) {
    b_pos <- which(rn == "b")
    th_pos <- which(rn == "theta")
    would_add <- if (is.null(newdata)) {
      !is.null(lp[["Z"]]) && (use_re || length(ed[["sm_blocks"]]) > 0L)
    } else {
      (use_re && length(ed[["re_parts"]]) > 0L) || length(ed[["sm_parts"]]) > 0L
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
      if (use_re && !is.null(lp[["Z"]])) {
        upd <- add_b_cols(A, coef_pos, lp[["Z"]], b_pos, th_pos)
        A <- upd$A
        coef_pos <- upd$coef_pos
      } else if (!use_re && length(ed[["sm_blocks"]]) && !is.null(lp[["Z"]])) {
        for (bk in ed[["sm_blocks"]]) {
          A <- Matrix::cbind2(A, lp[["Z"]][, bk[["c_idx"]], drop = FALSE])
          coef_pos <- c(coef_pos, b_pos[bk[["b_idx"]]])
        }
      }
    } else {
      if (use_re && length(ed[["re_parts"]])) {
        Zn <- re_design_matrix(ed[["re_parts"]], ed[["n"]],
                               object$frame[["n_c"]] %||% length(est[["b"]]))
        upd <- add_b_cols(A, coef_pos, Zn, b_pos, th_pos)
        A <- upd$A
        coef_pos <- upd$coef_pos
      }
      for (sp in ed[["sm_parts"]]) {
        A <- Matrix::cbind2(A, sp$Xr)
        coef_pos <- c(coef_pos, b_pos[sp$bk[["b_idx"]]])
      }
    }
  }
  list(A = A, coef_pos = coef_pos)
}

#' Variance sources that are not coefficient uncertainty. A new grouping
#' level (`allow_new_levels`) contributes its block's marginal variance,
#' the population-prediction-interval convention; for `|ID|`-merged blocks
#' that is the JOINT block's slice for this component. An exact `gp()` at
#' an unseen position contributes the GP's own conditional variance
#' (zero at observed positions).
#'
#' @noRd
lp_extra_var <- function(object, ed, use_re) {
  nl <- list()
  if (use_re && length(ed[["re_parts"]])) {
    th <- object$estimates[["theta"]]
    for (rp in ed[["re_parts"]]) {
      nas <- which(is.na(rp$j))
      if (!length(nas)) next
      bk <- rp$bk
      # the levels ARE the structure there, so there is no marginal
      # variance to hand an unseen one
      if (bk[["covstruct"]] %in% c("gr_cov", "gr_prec", "car", "spde")) next
      S <- covstruct_registry[[bk[["covstruct"]]]]$vcov(th[bk[["theta_idx"]]],
        bk)
      # a Student-t block's registry vcov() is the SCALE matrix; the
      # variance an unseen level contributes is nu/(nu-2) times it. The
      # interval is still built as a gaussian one, so it is the right
      # variance around a heavier-tailed truth, not the right quantile
      if (is_student_block(bk)) S <- S * student_var_factor(bk[["dist_nu"]])
      nl[[length(nl) + 1L]] <- list(
        bk = bk,
        S = S,
        cols = rp$comp$offset + seq_len(rp$comp$dim),
        mm = rp$mm, nas = nas, new_key = rp$new_key
      )
    }
  }
  gp <- list()
  for (sp in ed[["sm_parts"]]) {
    if (!is.null(sp$extra_var)) gp[[length(gp) + 1L]] <- sp$extra_var
  }
  list(new_levels = nl, gp = gp)
}

#' Collect the new-level variance sources into INDEPENDENT DRAWS.
#'
#' Two `lp_extra_var()` entries load the SAME draw exactly when they
#' name the same block and, on that row, the same level of it. Summing
#' their design rows first and taking one quadratic form is then the
#' right answer, and it is not the same as adding the two quadratic
#' forms: the shared draw contributes `(w1 + w2)^2 S`, two distinct
#' draws contribute `w1^2 S + w2^2 S`.
#'
#' Three cases meet here, and the grouping key settles all three.
#' Components of one `|ID|`-merged block, and one block appearing in
#' several linear predictors, always name one level per row, so they
#' share a key and their cross-covariance is kept. A multi-membership
#' term is the case where one block contributes SEVERAL entries per row
#' - one per member - and whether two of them are one draw depends on
#' the data: `allow_new_levels = TRUE` on a row whose two members carry
#' the same unseen label is one draw, two different unseen labels are
#' two. `new_key` carries the label that decides it; entries without
#' one (every single-membership term) key on a constant, which
#' reproduces the one-level-per-row grouping.
#'
#' @noRd
extra_var_blocks <- function(nl, n, weights = NULL) {
  out <- list()
  for (idx in seq_along(nl)) {
    e <- nl[[idx]]
    if (!length(e$nas)) next
    kv <- e$new_key %||% rep(".", n)
    w <- if (is.null(weights)) rep(1, n) else weights[[idx]]
    bkey <- as.character(e$bk[["c_idx"]][1L])
    for (lev in unique(kv[e$nas])) {
      rows <- e$nas[!is.na(kv[e$nas]) & kv[e$nas] == lev]
      if (!length(rows)) next
      key <- paste0(bkey, "\r", lev)
      B <- out[[key]] %||% list(S = e$S, M = matrix(0, n, e$bk[["dim"]]),
                                rows = integer(0))
      B$M[rows, e$cols] <- B$M[rows, e$cols] +
        w[rows] * e$mm[rows, , drop = FALSE]
      B$rows <- union(B$rows, rows)
      out[[key]] <- B
    }
  }
  out
}

#' Central-difference gradient of the expected response with respect to
#' one dpar's linear predictor. Analytic gradients exist for the simple
#' `mean_fn` forms, but the family set (and the truncated means) is wide
#' enough that one differencing rule beats a table of hand derivatives;
#' every `mean_fn` is elementwise, so a whole column of the Jacobian costs
#' two evaluations. The step is relative so it survives both tiny and
#' large etas.
#'
#' @noRd
mean_eta_grad <- function(fam, dp, av, dnm, link, eta) {
  h <- 1e-5 * pmax(1, abs(eta))
  dp_hi <- dp
  dp_lo <- dp
  dp_hi[[dnm]] <- link$linkinv(eta + h)
  dp_lo[[dnm]] <- link$linkinv(eta - h)
  (response_mean(fam, dp_hi, av) - response_mean(fam, dp_lo, av)) / (2 * h)
}

#' Delta-method SEs for the expected response of a family whose mean is
#' not the mu dpar (zero-inflation, hurdles, lognormal, trials-binomial,
#' or any truncated response). The mean runs through EVERY dpar's linear
#' predictor, so the gradient row stacks `dm/deta_k` times each predictor's
#' own A matrix and the quadratic form is taken over the JOINT
#' coefficient covariance: the cross-dpar covariances (and the shared b
#' block) are part of the answer, not an afterthought. Random effects
#' enter conditional on their modes, the same convention `predict()` uses
#' for eta.
#'
#' @noRd
predict_mean_se <- function(object, rspec, newdata, use_re,
                            allow_new_levels) {
  fam <- rspec$family
  rnm <- rspec$resp_name
  jc <- get_joint_cov(object)
  has_rr <- isTRUE(object$frame[["has_rr"]])
  rrj <- if (has_rr) rr_jacobians(object)
  dnames <- names(rspec$dpars)
  eds <- list()
  das <- list()
  evs <- list()
  dp <- list()
  for (dnm in dnames) {
    lp <- object$frame[["linpreds"]][[linpred_key(rnm, dnm)]]
    if (!is.null(lp[["nl_body"]])) {
      stop("se.fit is not supported on the response scale for a ",
           "nonlinear predictor yet; request the nonlinear parameters ",
           "(dpar = '", rspec$nlpars[1], "', ...) instead", call. = FALSE)
    }
    ed <- lp_eta_design(object, lp, newdata, use_re, allow_new_levels)
    eds[[dnm]] <- ed
    das[[dnm]] <- lp_delta_A(object, lp, ed, newdata, use_re, jc,
                             has_rr, rrj)
    evs[[dnm]] <- lp_extra_var(object, ed, use_re)
    dp[[dnm]] <- lp[["link"]]$linkinv(ed[["eta"]])
  }
  n <- eds[[1L]]$n
  av <- if (is.null(newdata)) {
    object$frame[["aterm_values"]][[rnm]]
  } else {
    aterms_for_newdata(rspec, newdata)
  }
  m <- response_mean(fam, dp, av)
  if (is.null(m) || !is.numeric(m) || !is.null(dim(m)) || length(m) != n) {
    stop("se.fit is not available for the expected response of family '",
         fam[["family"]], "': its mean is not one number per observation",
         call. = FALSE)
  }

  grad <- list()
  for (dnm in dnames) {
    lp <- object$frame[["linpreds"]][[linpred_key(rnm, dnm)]]
    grad[[dnm]] <- mean_eta_grad(fam, dp, av, dnm, lp[["link"]], eds[[dnm]]$eta)
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
  nl_all <- list()
  nl_grad <- list()
  for (dnm in dnames) {
    for (nl in evs[[dnm]]$new_levels) {
      nl_all[[length(nl_all) + 1L]] <- nl
      nl_grad[[length(nl_grad) + 1L]] <- grad[[dnm]]
    }
  }
  for (B in extra_var_blocks(nl_all, n, nl_grad)) {
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
  formula$frame[["data_frame"]]
}

#' Reinsert NAs for `na.exclude` fits (`napredict` is a no-op for
#' `na.omit`).
#'
#' @noRd
napred <- function(fit, x) {
  stats::napredict(fit$frame[["na_action"]], x)
}

#' Fitted values
#'
#' The modelled response at the estimates, conditional on the
#' random-effect modes. Equal to `predict(object, type = "response")` on
#' the training data, for every family.
#'
#' @param object A `frmtmb_fit`.
#' @param ... Unused.
#' @return A numeric vector of expected responses; for an ordinal family
#'   (`cumulative()`, `sratio()`, `cratio()`, `acat()`) an `n x K` matrix
#'   of category probabilities.
#' @section Ordinal responses:
#' An ordinal response has no mean, so `fitted()` returns the `n x K`
#' matrix of category probabilities, with the response's own factor
#' levels as column names and rows summing to one - the brms `fitted()`
#' convention. `cs()` terms are honored. The latent linear predictor,
#' which is where the coefficients live and where `se.fit` is available,
#' is `predict(object, type = "link")`.
#' @seealso [predict.frmtmb_fit()], [residuals.frmtmb_fit()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100))
#' dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x))
#' fit <- frm(bf(y ~ x) + poisson(), data = dd)
#' max(abs(fitted(fit) - predict(fit, type = "response")))
#' @export
fitted.frmtmb_fit <- function(object, ...) {
  rspec <- single_response(object, "fitted()")
  # a structured family whose row mean conditions on the whole observed
  # response supplies it here; everything else is rowwise
  fm <- fam_structure(rspec$family)[["fitted_mean"]]
  if (!is.null(fm)) {
    return(napred(object,
                  fm(object, frame_block_of(object$frame,
                                            rspec$resp_name))))
  }
  if (identical(rspec$family[["type"]], "ordinal")) {
    # no mean exists; the modelled response IS the category
    # distribution, and predict(type = "response") agrees by construction
    return(predict_ordinal(object, rspec, NULL, TRUE, FALSE))
  }
  if (identical(rspec$family[["type"]], "categorical")) {
    # same reasoning: the modelled response IS the category distribution
    return(predict_categorical(object, rspec, NULL, TRUE, FALSE))
  }
  dp <- eval_dpars(object)[[rspec$resp_name]]
  if (!"mu" %in% names(dp) && is.null(rspec$family[["post"]]$mean_fn)) {
    stop("fitted() is not defined for family '", rspec$family[["family"]], "'",
         call. = FALSE)
  }
  fam <- rspec$family
  out <- response_mean(fam, dp,
                       object$frame[["aterm_values"]][[rspec$resp_name]])
  napred(object, out)
}

#' Number of ordinal categories, from the threshold vector rather than
#' from the data: the top category may be unobserved.
#'
#' @noRd
ordinal_ncat <- function(fit) {
  raw <- fit$estimates[["tau_raw"]]
  if (is.null(raw)) {
    rspec <- single_response(fit, "residuals()")
    return(max(fit$frame[["y"]][[rspec$resp_name]]))
  }
  length(raw) + 1L
}

#' The `n x (K-1)` matrix of threshold-specific offsets a `cs()` term
#' contributes, or NULL when the predictor has none. In sample the
#' column values were kept at frame time; on newdata the term has to be
#' re-evaluated, which is what `label` carries (it is `"cs"` followed by
#' the deparsed expression).
#'
#' @noRd
ord_cs_values <- function(object, lp, newdata, n) {
  cst <- lp[["cs"]] %||% list()
  if (!length(cst)) return(list())
  env <- object$spec$responses[[lp[["resp"]]]]$formula_env
  lapply(cst, function(ct) {
    v <- if (is.null(newdata)) {
      ct$vals
    } else {
      ex <- ct$expr %||% str2lang(sub("^cs", "", ct$label))
      as.numeric(eval(ex, newdata, env))
    }
    if (length(v) == 1L) v <- rep(v, n)
    if (length(v) != n) {
      stop("cs() term '", ct$label, "' evaluated to ", length(v),
           " value(s) on ", n, " rows of newdata", call. = FALSE)
    }
    list(par = ct$par, vals = v, label = ct$label)
  })
}

#' Category-specific offsets for an ordinal linear predictor.
#'
#' Returns the `n x (K - 1)` matrix that `cs()` terms add to the
#' thresholds, or `NULL` when the model has none. It is separate from
#' the threshold arithmetic because the offset is the only part that
#' depends on `newdata`.
#'
#' @noRd
ord_cs_offsets <- function(object, lp, newdata, n, K1) {
  cv <- ord_cs_values(object, lp, newdata, n)
  if (!length(cv)) return(NULL)
  CS <- matrix(0, n, K1)
  for (ct in cv) CS <- CS + outer(ct$vals, object$estimates[[ct$par]])
  CS
}

#' `n x K` category probabilities of an ordinal fit.
#'
#' The probabilities come out of the family's OWN log-density, one
#' category at a time: `P(y = k) = exp(lpdf(k, eta, thresholds))`. The
#' four ordinal lpdfs are proper pmfs on `1..K`, so this is exact, it
#' cannot drift away from the likelihood the model was fitted with, and
#' it needs no second copy of the cumulative / sequential /
#' adjacent-category algebra (it agrees with the simulators'
#' `ord_cat_probs()` to machine precision, which the tests assert). A
#' custom ordinal family gets the same treatment for free.
#'
#' @noRd
ord_probs_from_eta <- function(fam, eta, cs, extra, K) {
  n <- length(eta)
  dp <- list(mu = eta)
  if (!is.null(cs)) dp[[".cs"]] <- cs
  P <- matrix(NA_real_, n, K)
  # ordinal lpdfs take the extras (thresholds) as a fourth argument;
  # categorical's takes three - dispatch on arity like fam_lcdf() does
  four <- length(formals(fam[["lpdf"]])) >= 4L
  for (k in seq_len(K)) {
    P[, k] <- exp(as.numeric(if (four) {
      fam[["lpdf"]](rep.int(k, n), dp, list(), extra)
    } else {
      fam[["lpdf"]](rep.int(k, n), dp, list())
    }))
  }
  # analytically the rows already sum to one; the division only removes
  # the last bit of rounding, and turns an overflowed row into NaN
  # instead of a silent zero vector
  P / rowSums(P)
}

#' `n x K` category probabilities in FITTED-row space (no `na.exclude`
#' padding), for the internal consumers that work alongside `y`.
#'
#' @noRd
ord_probs <- function(object, rspec, newdata = NULL, use_re = TRUE,
                      allow_new_levels = FALSE) {
  fam <- rspec$family
  lp <- object$frame[["linpreds"]][[linpred_key(rspec$resp_name, "mu")]]
  if (!is.null(lp[["nl_body"]])) {
    stop("type = \"response\" is not supported for an ordinal family ",
         "with a nonlinear predictor", call. = FALSE)
  }
  ed <- lp_eta_design(object, lp, newdata, use_re, allow_new_levels)
  eta <- unname(ed[["eta"]])
  n <- length(eta)
  K <- ordinal_ncat(object)
  cs <- ord_cs_offsets(object, lp, newdata, n, K - 1L)
  # the ordinal lpdfs read only `extra` (the thresholds and the cs
  # coefficients); no addition term enters a category probability
  P <- ord_probs_from_eta(fam, eta, cs, fit_extras(object), K)
  colnames(P) <- object$frame[["y_levels"]][[rspec$resp_name]] %||%
    as.character(seq_len(K))
  rn <- names(ed[["eta"]])
  if (is.null(rn) && is.null(newdata)) {
    rn <- rownames(object$frame[["data_frame"]])
  }
  if (!is.null(rn) && length(rn) == n) rownames(P) <- rn
  # a row that cannot be estimated from the retained design columns has
  # no category distribution either
  if (any(ed[["nonest"]])) P[ed[["nonest"]], ] <- NA_real_
  P
}

#' Response-scale prediction for an ordinal family.
#'
#' Wraps `ord_probs()` and returns the `n x K` matrix of category
#' probabilities, padded by `napred()` on the training data so that
#' `na.exclude` gives back the input row count. The modelled response of
#' an ordinal family is a category distribution, not a mean, which is
#' why this is not the general response path.
#'
#' @noRd
predict_ordinal <- function(object, rspec, newdata, use_re,
                            allow_new_levels) {
  P <- ord_probs(object, rspec, newdata, use_re, allow_new_levels)
  if (is.null(newdata)) P <- napred(object, P)
  P
}

#' `n x K` category probabilities of a categorical fit, in fitted-row
#' space.
#'
#' Unlike the ordinal families, which share one latent predictor, a
#' categorical fit carries K-1 of them - one per non-reference category,
#' in category order - so the softmax runs over the whole set. The
#' reference category's column is the pinned zero, and the row maximum
#' comes out before `exp()` so a wide predictor cannot overflow.
#'
#' @noRd
cat_probs <- function(object, rspec, newdata = NULL, use_re = TRUE,
                      allow_new_levels = FALSE) {
  dpn <- rspec$primary_dpars
  K <- length(dpn) + 1L
  E <- NULL
  nonest <- NULL
  rn <- NULL
  for (j in seq_along(dpn)) {
    lp <- object$frame[["linpreds"]][[linpred_key(rspec$resp_name, dpn[j])]]
    if (!is.null(lp[["nl_body"]])) {
      stop("type = \"response\" is not supported for a categorical ",
           "family with a nonlinear predictor", call. = FALSE)
    }
    ed <- lp_eta_design(object, lp, newdata, use_re, allow_new_levels)
    if (is.null(E)) {
      n <- length(ed[["eta"]])
      E <- matrix(0, n, K)
      nonest <- rep(FALSE, n)
      rn <- names(ed[["eta"]])
    }
    E[, j + 1L] <- unname(ed[["eta"]])
    nonest <- nonest | ed[["nonest"]]
  }
  P <- exp(E - apply(E, 1L, max))
  P <- P / rowSums(P)
  colnames(P) <- object$frame[["y_levels"]][[rspec$resp_name]] %||%
    as.character(seq_len(K))
  if (is.null(rn) && is.null(newdata)) {
    rn <- rownames(object$frame[["data_frame"]])
  }
  if (!is.null(rn) && length(rn) == nrow(P)) rownames(P) <- rn
  if (any(nonest)) P[nonest, ] <- NA_real_
  P
}

#' Response-scale prediction for a categorical family.
#'
#' The counterpart of `predict_ordinal()` for a nominal response: wraps
#' `cat_probs()` and pads with `napred()` on the training data. The two
#' stay separate because the ordinal probabilities come from thresholds
#' and the categorical ones from one linear predictor per category.
#'
#' @noRd
predict_categorical <- function(object, rspec, newdata, use_re,
                                allow_new_levels) {
  P <- cat_probs(object, rspec, newdata, use_re, allow_new_levels)
  if (is.null(newdata)) P <- napred(object, P)
  P
}

#' Delta-method standard errors of the `n x K` category probabilities on
#' a prediction grid.
#'
#' A category probability runs through the linear predictor AND through
#' the thresholds (and the `cs()` coefficients), so the gradient is
#' taken over all of them jointly and the quadratic form uses the joint
#' covariance: thresholds are estimated too, and pretending otherwise
#' would understate every band. The eta part reuses `lp_delta_A()`, so
#' the coefficient bookkeeping is exactly `predict(se.fit = TRUE)`'s;
#' the derivative of `p_k` with respect to eta and with respect to each
#' extra parameter is a central difference of the family's own lpdf,
#' the same differencing rule `mean_eta_grad()` uses and for the same
#' reason (a custom ordinal family gets it for free).
#'
#' @noRd
ord_prob_se <- function(object, rspec, lp, ed, newdata, use_re,
                        weights = NULL) {
  fam <- rspec$family
  K <- ordinal_ncat(object)
  eta <- unname(ed[["eta"]])
  n <- length(eta)
  extra <- fit_extras(object)
  csv <- ord_cs_values(object, lp, newdata, n)
  CS <- if (length(csv)) {
    M <- matrix(0, n, K - 1L)
    for (ct in csv) M <- M + outer(ct$vals, object$estimates[[ct$par]])
    M
  }
  probs <- function(e, cs, ex) ord_probs_from_eta(fam, e, cs, ex, K)
  P0 <- probs(eta, CS, extra)

  jc <- get_joint_cov(object)
  has_rr <- isTRUE(object$frame[["has_rr"]])
  rrj <- if (has_rr) rr_jacobians(object)
  da <- lp_delta_A(object, lp, ed, newdata, use_re, jc, has_rr, rrj)
  A <- as.matrix(da$A)
  pos <- da$coef_pos

  h <- 1e-5 * pmax(1, abs(eta))
  dPde <- (probs(eta + h, CS, extra) - probs(eta - h, CS, extra)) /
    (2 * h)

  # one n x K derivative block per estimated extra parameter, paired
  # with its row in the joint covariance
  extra_d <- list()
  extra_pos <- integer(0)
  cs_par <- vapply(csv, `[[`, "", "par")
  for (nm in names(extra)) {
    ep <- which(jc$names == nm)
    val <- extra[[nm]]
    ci <- match(nm, cs_par)
    for (j in seq_along(val)) {
      if (j > length(ep)) next
      hj <- 1e-5 * max(1, abs(val[j]))
      d <- if (is.na(ci)) {
        ehi <- extra
        elo <- extra
        ehi[[nm]][j] <- val[j] + hj
        elo[[nm]][j] <- val[j] - hj
        (probs(eta, CS, ehi) - probs(eta, CS, elo)) / (2 * hj)
      } else {
        chi <- CS
        clo <- CS
        chi[, j] <- chi[, j] + hj * csv[[ci]]$vals
        clo[, j] <- clo[, j] - hj * csv[[ci]]$vals
        (probs(eta, chi, extra) - probs(eta, clo, extra)) / (2 * hj)
      }
      extra_d[[length(extra_d) + 1L]] <- d
      extra_pos <- c(extra_pos, ep[j])
    }
  }
  n_beta <- ncol(A)
  V <- jc$V[c(pos, extra_pos), c(pos, extra_pos), drop = FALSE]
  # `weights` collapses the K columns into ONE displayed quantity,
  # sum_k w_k p_k - the expected category number for w = 1..K. The
  # weights go on the GRADIENT before the quadratic form, so the
  # covariances between the category probabilities are kept, which
  # summing K separate standard errors would throw away.
  if (!is.null(weights)) {
    stopifnot(length(weights) == K)
    P0 <- matrix(as.vector(P0 %*% weights), n, 1L)
    dPde <- matrix(as.vector(dPde %*% weights), n, 1L)
    extra_d <- lapply(extra_d, function(d) {
      matrix(as.vector(d %*% weights), n, 1L)
    })
  }
  nq <- ncol(P0)
  SE <- matrix(NA_real_, n, nq)
  G <- matrix(0, n, n_beta + length(extra_d))
  for (k in seq_len(nq)) {
    G[, seq_len(n_beta)] <- dPde[, k] * A
    for (i in seq_along(extra_d)) G[, n_beta + i] <- extra_d[[i]][, k]
    SE[, k] <- sqrt(pmax(rowSums((G %*% V) * G), 0))
  }
  if (any(ed[["nonest"]])) {
    P0[ed[["nonest"]], ] <- NA_real_
    SE[ed[["nonest"]], ] <- NA_real_
  }
  if (is.null(weights)) {
    colnames(P0) <- colnames(SE) <-
      object$frame[["y_levels"]][[rspec$resp_name]] %||%
        as.character(seq_len(K))
  }
  list(P = P0, se = SE)
}

#' Mean and variance of the CATEGORY INDEX under the fitted category
#' distribution, `E[Y] = sum_k k p_k` and `Var[Y]`, in fitted-row space.
#'
#' An ordinal response has no mean, but every consumer that needs one
#' number per observation (a residual, a residuals-versus-fitted plot,
#' DHARMa's `fittedPredictedResponse`) needs one anyway. Scoring the
#' categories by their own integer codes is the standard fallback: it is
#' what `y - E[Y]` means in brms's `residuals()` (which subtracts drawn
#' categories from the observed ones) and it is monotone in the latent
#' predictor, which is all the plots and the rank transform use it for.
#'
#' @noRd
ord_cat_moments <- function(object, rspec) {
  P <- ord_probs(object, rspec)
  k <- seq_len(ncol(P))
  m <- as.numeric(P %*% k)
  v <- as.numeric(P %*% (k^2)) - m^2
  list(mean = m, var = v, P = P)
}

#' OSA integration window and row split for a censored response, or NULL
#' when nothing is censored.
#'
#' A censored row contributes a probability MASS: on the tape its
#' likelihood no longer depends on the observation, so `oneStepPredict`
#' either inverts a singular system (`fullGaussian`) or integrates a flat
#' slice to infinity (`oneStepGeneric`). Both are real: the "observation"
#' on such a row is an event, not a value, and has no one-step CDF.
#' What IS well defined is the CDF of the uncensored rows conditional on
#' the censoring events, which is what subset/conditional buys. Under
#' type-I censoring (one censoring point per side) an uncensored row is
#' exactly a draw that landed inside the window, so its PIT renormalizes
#' on that window just as a `trunc()` fit's does.
#'
#' @noRd
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
#' @section Ordinal responses:
#' An ordinal response has no mean, so `"response"` and `"pearson"`
#' score the categories by the integer codes `1..K` the likelihood
#' itself uses: `"response"` is `y - E[Y]` with
#' `E[Y] = sum_k k * P(y = k)` taken from [fitted()]'s category
#' probabilities, and `"pearson"` divides by the standard deviation of
#' that same distribution. This is the frequentist point-estimate form
#' of what brms's `residuals()` reports on an ordinal fit (there, the
#' observed category minus a drawn one). It is a residual on a SCORE,
#' not on the ordinal scale, so read it for gross lack of fit and
#' pattern, not as a calibrated quantity: `"osa"` and
#' [dharma_residuals()] give residuals that use only the order.
#' `"deviance"` is refused, as it is for every family without a
#' standard unit deviance.
#'
#' `"osa"` uses `"oneStepGeneric"` over the discrete support `1..K`,
#' which makes the residuals randomized quantile residuals.
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
#' @section Residual correlation terms:
#' Under an `ar()`, `ma()`, `arma()`, `cosy()` or `unstr()` term (see
#' [frmtmb-autocor]) the residual covariance of a group is `D R D` with
#' `R` unit-diagonal, so the marginal residual SD of a row is still
#' `sigma` and `"pearson"` is unchanged - it divides by exactly that.
#' What `"response"` and `"pearson"` do NOT do is decorrelate: plotted
#' against time within a group they still show the fitted
#' autocorrelation, which is the intended reading. `"osa"` is refused,
#' because the taped likelihood is a joint density per group rather
#' than a product of per-observation terms; [dharma_residuals()] works,
#' since [simulate()] draws one correlated residual per group.
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
#'
#' @srrstats {G2.2} Parameters that expect a univariate response refuse a
#'   multivariate fit rather than silently using the first response. One
#'   guard serves every such method (`residuals()`, `fitted()`,
#'   `simulate()`, `dharma_residuals()`, `mixture_probs()`, `pp_check()`),
#'   and errors naming the method that is not yet multivariate.
#' @srrstats {RE4.10} Model residuals are returned by `residuals()`, in
#'   four types, with enough documentation to interpret them and to hand
#'   them to a user's own test. `"response"` and `"pearson"` are the usual
#'   forms; `"deviance"` is defined per family in a dedicated section that
#'   lists the supported families and distinguishes the residual from the
#'   model deviance reported by `deviance()`; `"osa"` gives one-step-ahead
#'   quantile residuals through [TMB::oneStepPredict()], which are
#'   standard normal under a correctly specified model whatever the
#'   family, and which stay valid for censored, truncated, and ordinal
#'   responses where Pearson residuals mislead. `dharma_residuals()`
#'   hands the simulation-based equivalent to DHARMa for the user's own
#'   tests, and `vignette("diagnostics")` works through the choice.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
#' fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
#'
#' # raw and variance-standardized residuals
#' head(residuals(fit))
#' head(residuals(fit, type = "pearson"))
#' # the usual overdispersion check for a poisson fit
#' sum(residuals(fit, type = "pearson")^2) / df.residual(fit)
#'
#' # one-step-ahead quantile residuals are standard normal under a
#' # correctly specified model, whatever the family
#' r <- residuals(fit, type = "osa")
#' qqnorm(r); qqline(r)
#' @export
residuals.frmtmb_fit <- function(object, type = c("response", "pearson",
                                                  "deviance", "osa"),
                                 osa_method = NULL, ...) {
  type <- match.arg(type)
  rspec <- single_response(object, "residuals()")
  fam <- rspec$family
  if (identical(fam[["type"]], "categorical")) {
    # nothing here is defined on a nominal scale: the categories carry
    # no order, so there is no y - E[Y] to form and no CDF to invert
    stop("residuals() is not defined for a categorical family: the ",
         "categories carry no order, so a residual has no scale to live ",
         "on. Compare fitted(fit) (the n x K category probabilities) ",
         "against the observed categories instead", call. = FALSE)
  }
  st <- fam_structure(fam)
  if (!is.null(st)) {
    # Both of these need a PER-ROW likelihood: one-step prediction needs
    # the taped density of one observation given the earlier ones, and
    # the unit deviance needs a row's likelihood to compare with its
    # saturated fit. A structured likelihood has neither.
    if (type == "osa") {
      structure_gate(st, "osa",
                     structure_generic(fam, "residuals(type = \"osa\")"))
    }
    if (type == "deviance") {
      structure_gate(st, "deviance",
                     structure_generic(fam,
                                       "residuals(type = \"deviance\")"))
    }
    fm <- st[["fitted_mean"]]
    if (!is.null(fm)) {
      blk <- frame_block_of(object$frame, rspec$resp_name)
      r <- object$frame[["y"]][[rspec$resp_name]] - fm(object, blk)
      if (type == "pearson") {
        fv <- st[["fitted_var"]]
        if (is.null(fv)) {
          stop("residuals(type = \"pearson\") needs the conditional ",
               "variance of each row given the whole response, and the ",
               "'", fam[["family"]], "' family declares no fitted_var(). Use ",
               "type = \"response\"", call. = FALSE)
        }
        r <- r / sqrt(fv(object, blk))
      }
      # a masked (NA) response has no residual, and the placeholder
      # value standing in for it on the tape must never look like one
      if (!is.null(blk[["miss"]])) r[blk[["miss"]]] <- NA_real_
      return(napred(object, r))
    }
  }
  if (type == "osa") {
    if (!is.null(object$frame[["autocor"]][[rspec$resp_name]])) {
      # oneStepPredict needs the taped density of ONE observation given
      # the previous ones; under an R-side residual the tape holds a
      # joint density per group and never registers an observation
      # vector (no OBS() call), so there is nothing to step through
      stop("residuals(type = \"osa\") is not available for a fit with a ",
           "residual correlation term (",
           object$frame[["autocor"]][[rspec$resp_name]]$label,
           "): the likelihood is a joint density per group, not a ",
           "product of per-observation terms. Use type = \"pearson\", ",
           "which divides by the marginal residual SD, or ",
           "dharma_residuals(), which uses simulate() and does draw ",
           "correlated residuals", call. = FALSE)
    }
    av0 <- object$frame[["aterm_values"]][[rspec$resp_name]]
    tb <- trunc_bounds(av0, object$frame[["n_obs"]])
    cb <- osa_cens_domain(av0, object$frame[["y"]][[rspec$resp_name]])
    ordinal <- identical(fam[["type"]], "ordinal")
    method <- osa_method %||%
      if (!is.null(tb) || !is.null(cb) || ordinal) "oneStepGeneric"
      else if (identical(fam[["family"]], "gaussian")) "fullGaussian"
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
      args$discrete <- identical(fam[["type"]], "discrete") || ordinal
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
      full <- rep(NA_real_, object$frame[["n_obs"]])
      full[cb$subset] <- r
      r <- full
    }
    return(napred(object, r))
  }
  if (identical(fam[["type"]], "ordinal") &&
      type %in% c("response", "pearson")) {
    # scored by the category codes the likelihood itself uses; see
    # ord_cat_moments()
    mom <- ord_cat_moments(object, rspec)
    r <- object$frame[["y"]][[rspec$resp_name]] - mom$mean
    if (type == "pearson") r <- r / sqrt(mom$var)
    return(napred(object, r))
  }
  dp <- eval_dpars(object)[[rspec$resp_name]]
  av <- object$frame[["aterm_values"]][[rspec$resp_name]]
  yv <- object$frame[["y"]][[rspec$resp_name]]
  if (type == "deviance") {
    return(napred(object, deviance_residuals(fam, yv, dp, av,
                                             object$frame[["n_obs"]])))
  }
  # on a truncated response the residual is against the truncated mean,
  # which is what the data were actually drawn from
  mu <- response_mean(fam, dp, av)
  r <- yv - mu
  if (type == "pearson") {
    if (is.null(fam[["post"]]$var_fn)) {
      stop("Family '", fam[["family"]], "' has no variance function; ",
           "pearson residuals are unavailable", call. = FALSE)
    }
    # the scale stays the untruncated family variance; only the centering
    # is truncation-aware, so pearson residuals on a truncated model are
    # slightly conservative
    v <- fam[["post"]]$var_fn(dp, av)
    r <- r / sqrt(v)
  }
  napred(object, r)
}

#' One draw of the full b vector from its estimated distribution
#' `N(0, Sigma)`, or the multivariate t with scale `Sigma` on a
#' `gr(dist = "student")` block.
#'
#' @noRd
draw_b <- function(fit) {
  th <- fit$estimates[["theta"]]
  b <- numeric(length(fit$estimates[["b"]] %||% numeric(0)))
  for (bk in fit$frame[["re_blocks"]]) {
    if (bk[["covstruct"]] == "rr") {
      # standard-normal factors; eval_dpars expands them via loadings
      b[bk[["b_idx"]]] <- stats::rnorm(length(bk[["b_idx"]]))
      next
    }
    if (bk[["covstruct"]] == "gr_cov") {
      # correlation is across levels, not within them
      S <- covstruct_registry[["gr_cov"]]$vcov(th[bk[["theta_idx"]]], bk)
      K <- kronecker(bk[["aux_A"]], S)
      b[bk[["b_idx"]]] <- drop(crossprod(chol(K),
                                    stats::rnorm(nrow(K))))
      next
    }
    if (bk[["covstruct"]] == "gr_prec") {
      # x = U^-1 z with U'U = Q has covariance Q^-1; for correlated
      # slopes Q is the level-major Kronecker precision Q (x) Sigma^-1
      Qb <- if (bk[["dim"]] == 1L) {
        exp(-2 * th[bk[["theta_idx"]][1]]) * bk[["aux_Q"]]
      } else {
        S <- covstruct_registry[["gr_prec"]]$vcov(th[bk[["theta_idx"]]], bk)
        Matrix::kronecker(bk[["aux_Q"]], methods::as(solve(unname(S)),
                                                "generalMatrix"))
      }
      U <- Matrix::chol(Qb)
      b[bk[["b_idx"]]] <- as.vector(
        Matrix::solve(U, stats::rnorm(length(bk[["b_idx"]]))))
      next
    }
    if (bk[["covstruct"]] == "car") {
      # the whole field is one draw from its (dense) covariance
      K <- car_cov(th[bk[["theta_idx"]]], bk)
      b[bk[["b_idx"]]] <- drop(crossprod(chol(K), stats::rnorm(nrow(K))))
      next
    }
    if (bk[["covstruct"]] == "spde") {
      U <- Matrix::chol(spde_prec(th[bk[["theta_idx"]]], bk))
      b[bk[["b_idx"]]] <- as.vector(Matrix::solve(U,
                                             stats::rnorm(bk[["n_levels"]])))
      next
    }
    V <- covstruct_registry[[bk[["covstruct"]]]]$vcov(th[bk[["theta_idx"]]], bk)
    L <- chol(V)
    U <- matrix(stats::rnorm(bk[["n_levels"]] * bk[["dim"]]),
      bk[["n_levels"]]) %*% L
    if (is_student_block(bk)) {
      # V is the SCALE matrix, and the mixing variable is per LEVEL -
      # shared across that level's coefficients, which is what makes the
      # draw a multivariate t rather than d independent ones
      nu <- bk[["dist_nu"]]
      U <- U * sqrt(nu / stats::rchisq(bk[["n_levels"]], df = nu))
    }
    b[bk[["b_idx"]]] <- as.vector(t(U))   # level-major
  }
  b
}

#' The observation window the censoring mechanism imposes: the censoring
#' point is the response value on a censored row, and it must be the
#' same for every censored row on a side (type-I censoring). With
#' row-varying censoring times the point is unknown on the rows that
#' were NOT censored (the data only say it was never reached), so the
#' mechanism cannot be applied to their draws at all.
#'
#' @noRd
cens_window <- function(av, yobs) {
  cen <- av[["cens"]]
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

#' Apply the censoring mechanism to one draw of the latent response:
#' every draw outside the observation window is recorded at the window's
#' edge, exactly as the observed data were.
#'
#' @noRd
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
#' @section Structured draws:
#' Most families draw each row on its own. Some cannot, and those go
#' through one implementation that [simulate()], `posterior_predict()`
#' and [frm_simulate()] all reach (see `sim_ctx` in
#' [frmtmb_family()]):
#' - a `mixture(groups = ~g)` draw takes one class per GROUP and then
#'   simulates each row from its group's component;
#' - a [mixture_mvn()] draw takes a class per row and then a
#'   multivariate normal with that class's own covariance;
#' - a residual correlation term (`ar()`, `ma()`, `cosy()`, ...) is one
#'   multivariate residual draw per group added to the mean predictor,
#'   so the draws carry the fitted autocorrelation;
#' - a family from an extension package draws through whatever its
#'   structure declares, by the same route.
#'
#' A structured draw covers whole sequences or groups, so `trunc()`
#' rejection cannot resample single rows within it (every structured
#' model refuses `trunc()` when the frame is assembled) and
#' `posterior_predict(newdata =)` is refused: the structure indexes the
#' rows the model was fitted on.
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
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
#' fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
#'
#' # one column per draw; the seed used is attached
#' sims <- simulate(fit, nsim = 5, seed = 42)
#' str(sims)
#' attr(sims, "seed")
#'
#' # re.form = NA redraws the group effects, which is the right choice
#' # for a parametric bootstrap over new groups
#' sims_m <- simulate(fit, nsim = 5, re.form = NA, seed = 42)
#' apply(sims_m, 2, var) > apply(sims, 2, var)
#'
#' # a posterior-predictive check by hand: does the fit reproduce the
#' # share of zeros in the data?
#' mean(dd$y == 0)
#' colMeans(simulate(fit, nsim = 20, seed = 1) == 0)
#' @export
simulate.frmtmb_fit <- function(object, nsim = 1, seed = NULL,
                                re.form = NULL, censored = FALSE, ...) {
  # nsim reaches vapply()/replicate() as a length, where a length-2 or
  # character value reports "invalid 'length' argument" and names
  # neither simulate() nor nsim
  check_count(nsim, "nsim", min = 1L)
  check_flag(censored, "censored")
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
    on.exit(assign(".Random.seed", saved_seed, envir = globalenv()),
            add = TRUE)
  }
  rspec <- single_response(object, "simulate()")
  fam <- rspec$family
  # the draw itself is the structure's own sim_ctx; what is decided here
  # is only what these two options would MEAN for a whole-block draw
  st <- fam_structure(fam)
  if (!is.null(st)) {
    if (!is.null(re.form)) {
      structure_gate(st, "re_form",
                     structure_generic(fam, "simulate(re.form =)"),
                     context = "simulate")
    }
    if (isTRUE(censored)) {
      structure_gate(st, "cens_trunc",
                     structure_generic(fam, "simulate(censored = TRUE)"),
                     context = "simulate")
    }
  }
  if (!sim_can(fam)) {
    stop("simulate(): family '", fam[["family"]], "' has no simulator yet",
         sim_note(fam), call. = FALSE)
  }
  marginal <- !is.null(re.form) && !inherits(re.form, "formula") &&
    is.na(re.form)
  n <- stats::nobs(object)
  out <- vector("list", nsim)
  av <- object$frame[["aterm_values"]][[rspec$resp_name]]
  cwin <- NULL
  if (isTRUE(censored)) {
    if (is.null(av[["cens"]])) {
      stop("simulate(censored = TRUE) needs a cens() response",
           call. = FALSE)
    }
    cwin <- cens_window(av, object$frame[["y"]][[rspec$resp_name]])
  }
  for (s in seq_len(nsim)) {
    b_use <- if (marginal && length(object$frame[["re_blocks"]])) {
      draw_b(object)
    } else {
      object$estimates[["b"]]
    }
    dp <- with_cs_offsets(object, rspec, eval_dpars(object,
                                                    b = b_use))
    dp <- dp[[rspec$resp_name]]
    out[[s]] <- sim_draw(sim_context(object, rspec, dp, aterms = av,
                                     n = n,
                                     extra = fit_extras(object)))
    if (!is.null(cwin)) out[[s]] <- apply_censoring(out[[s]], cwin)
  }
  names(out) <- paste0("sim_", seq_len(nsim))
  out <- lapply(out, function(v) sim_restore_type(object, rspec, v))
  out <- sim_as_data_frame(out)
  attr(out, "seed") <- rng_state
  out
}

#' @rdname frmtmb-extension-api
#' @export
fit_extras <- function(fit) {
  nms <- fit$frame[["extra_names"]] %||% character(0)
  if (!length(nms)) return(NULL)
  fit$estimates[nms]
}

#' @rdname frmtmb-extension-api
#' @export
dpar_linpred <- function(frame, params, resp, dpar) {
  lp <- frame[["linpreds"]][[linpred_key(resp, dpar)]]
  if (is.null(lp) || !ncol(lp[["X"]])) return(NULL)
  eta <- drop(as.matrix(lp[["X"]] %*% params[[lp[["par"]]]][lp[["idx"]]]))
  if (!is.null(lp[["offset"]])) eta <- eta + lp[["offset"]]
  eta
}

#' `cs(x)` contributes an `n x (K-1)` matrix of threshold-specific
#' offsets. The objective builds it on the tape; `eval_dpars()` has no
#' reason to, so add it here for the ordinal simulators that consume it.
#'
#' @noRd
with_cs_offsets <- function(fit, rspec, dpv) {
  for (lp in fit$frame[["linpreds"]]) {
    if (!length(lp[["cs"]] %||% list())) next
    CS <- 0
    for (ct in lp[["cs"]]) {
      CS <- CS + outer(ct$vals, fit$estimates[[ct$par]])
    }
    dpv[[lp[["resp"]]]][[".cs"]] <- CS
  }
  dpv
}

#' `simulate()` hands draws back in the response's own type: an ordered
#' factor for an ordinal fit (the 1..K codes mean nothing without the
#' levels) and a matrix for a matrix response. `na.exclude` fits pad back
#' to the original row count, the same contract `fitted()` and
#' `residuals()` keep. `[glmmTMB test-simulate.R; lme4#737]`
#'
#' @noRd
sim_restore_type <- function(fit, rspec, v) {
  lv <- fit$frame[["y_levels"]][[rspec$resp_name]]
  if (!is.null(lv)) {
    # a categorical response's levels are nominal: ordering the draws
    # would claim an order the model never used
    v <- factor(lv[v], levels = lv,
                ordered = !identical(rspec$family[["type"]], "categorical"))
  } else if (is.matrix(v)) {
    yv <- fit$frame[["y"]][[rspec$resp_name]]
    if (is.matrix(yv) && !is.null(colnames(yv))) colnames(v) <- colnames(yv)
  }
  napred(fit, v)
}

#' A matrix response needs a data frame whose COLUMNS are matrices (the
#' lme4 convention); the default `as.data.frame()` would flatten each draw
#' into one column per category.
#'
#' @noRd
sim_as_data_frame <- function(out) {
  if (!is.matrix(out[[1L]])) return(as.data.frame(out))
  df <- data.frame(row.names = seq_len(nrow(out[[1L]])))
  for (nm in names(out)) df[[nm]] <- out[[nm]]
  df
}

#' Drop the rows `napredict()` padded back in, for the internal consumers
#' that work in fitted-row space (bootstrap refits, DHARMa, `pp_check`).
#'
#' @noRd
na_unpad <- function(fit, x) {
  na <- fit$frame[["na_action"]]
  if (is.null(na) || !inherits(na, "exclude")) return(x)
  idx <- unclass(na)
  if (is.matrix(x) || is.data.frame(x)) x[-idx, , drop = FALSE] else x[-idx]
}

#' The joint covariance of the fixed and random coefficients
#'
#' The covariance of everything the fit estimates, `beta`, `betad`,
#' `theta` and the random-effect coefficients `b` together, in one
#' matrix. It is what a delta method over a fitted CURVE needs and what
#' [vcov()] cannot return: a penalized smooth's wiggly part is a
#' random-effect block even when the smooth is a population term, so a
#' covariance that stops at the fixed effects covers none of it.
#'
#' `vcov(full = TRUE)` returns the OUTER parameter vector's covariance
#' and its row names are documented to be exactly [confint()]'s, so `b`
#' is not in it and cannot be added. `frm_lp_basis()` is the accessor
#' for a design at `newdata`; this is the accessor for the covariance
#' those coefficients have.
#'
#' The result is memoized on the fit, so the joint-precision solve is
#' paid once however many curves are drawn from it. It is also the ONLY
#' route to the covariance of an AUTOSCALED fit
#' (`frmtmb_control(autoscale = TRUE)`): a fresh
#' `RTMB::sdreport(getJointPrecision = TRUE)` goes round the
#' reparameterization and returns a covariance built on the unscaled
#' Hessian, which is a different matrix and is not marked as one.
#'
#' @param object A fitted `frmtmb_fit`.
#' @return A list with
#' \describe{
#'   \item{`V`}{the `p x p` joint covariance.}
#'   \item{`names`}{length `p`; the PARAMETER COMPONENT each row belongs
#'     to (`"beta"`, `"betad"`, `"b"`, `"theta"`, ...), which is what
#'     [frm_lp_basis()]`$coef_pos` indexes.}
#'   \item{`labels`}{length `p`; one label per row,
#'     `beta.<coefficient>` for a fixed effect and `b.<block>.<level>`
#'     for a random one.}
#' }
#' A fit with no random effects has no joint precision, and `V` is then
#' the fixed-effect covariance `sdreport()$cov.fixed`; `names` says so.
#' @seealso [frm_lp_basis()], [vcov()], [frmtmb-extension-api]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(120), g = factor(rep(1:12, each = 10)))
#' dd$y <- rnorm(120, 1 + 2 * dd$x + rnorm(12, 0, 0.5)[dd$g], 0.4)
#' fit <- frm(bf(y ~ x + (1 | g)), data = dd)
#' jc <- frm_joint_cov(fit)
#' dim(jc$V)
#' table(jc$names)
#' head(jc$labels)
#' @export
frm_joint_cov <- function(object) {
  require_frmtmb_fit(object, "frm_joint_cov()")
  require_fitted(object, "frm_joint_cov()")
  jc <- get_joint_cov(object)
  jc$labels <- joint_coef_labels(object, jc)
  jc
}

#' One label per row of the joint covariance.
#'
#' The component names alone do not say WHICH coefficient a row is, and
#' a caller assembling a curve has to line its design columns up against
#' something. Built here rather than cached with `V` so that the
#' memoized object keeps the shape every internal consumer already
#' reads.
#'
#' @noRd
joint_coef_labels <- function(fit, jc = get_joint_cov(fit)) {
  rn <- jc$names
  if (is.null(rn)) return(NULL)
  tpl <- fit$frame[["par_template"]]
  out <- rn
  for (cp in unique(rn)) {
    i <- which(rn == cp)
    nm <- if (cp == "b") b_coef_labels(fit) else names(tpl[[cp]])
    if (is.null(nm) || length(nm) != length(i)) {
      nm <- as.character(seq_along(i))
    }
    out[i] <- paste0(cp, ".", nm)
  }
  out
}

#' Per-coefficient labels of the `b` vector, block by block.
#'
#' @noRd
b_coef_labels <- function(fit) {
  out <- character(0)
  for (bk in fit$frame[["re_blocks"]] %||% list()) {
    lev <- bk[["levels"]]
    nl <- bk[["n_levels"]]
    d <- if (identical(bk[["covstruct"]], "rr")) bk[["rank"]] else bk[["dim"]]
    cn <- bk[["cnms"]]
    if (length(cn) != d) cn <- as.character(seq_len(d))
    # a smooth or gp block is one unnamed "level", so it gets no level
    # segment rather than an empty one
    seg <- if (nl == 1L && is.null(lev)) rep("", d) else {
      lv <- if (length(lev) == nl) lev else as.character(seq_len(nl))
      paste0(rep(lv, each = d), ".")
    }
    out <- c(out, paste0(bk[["term_label"]], ".", seg, rep(cn, nl)))
  }
  out
}

#' The design of a linear predictor over the coefficient vector
#'
#' `predict(se.fit = TRUE)` builds a matrix `A` with one row per
#' prediction and one column per contributing coefficient, forms
#' `A V A'` and keeps only its diagonal. Every delta-method quantity
#' over a fitted curve needs the whole thing: a contrast between two
#' grids, an average marginal effect with a correct standard error, a
#' simultaneous band, a derivative, the time of a peak. This returns the
#' pieces so that an extension does not have to rebuild `A` by
#' perturbation, one `predict()` call per coefficient.
#'
#' The name follows `emmeans::emm_basis()`, which is the same idea for
#' the fixed block alone.
#'
#' @param object A fitted `frmtmb_fit`.
#' @param newdata Data to build the design at, or `NULL` for the
#'   training data.
#' @param dpar,resp The distributional parameter and response to take
#'   the linear predictor of. Both default the way [predict()] defaults
#'   them.
#' @param re.form `NULL` keeps every random effect, `NA` drops them all,
#'   a one-sided formula keeps the ones it names.
#' @param allow_new_levels Whether a grouping level the fit never saw is
#'   allowed.
#' @return A list with
#' \describe{
#'   \item{`eta`}{the linear predictor, exactly
#'     `predict(type = "link")`.}
#'   \item{`A`}{`n x p`; `d eta / d coef`.}
#'   \item{`coef_pos`}{length `p`; the rows of `V` the columns of `A`
#'     belong to, in `V`'s own order.}
#'   \item{`V`}{`p x p`; [frm_joint_cov()] subset to `coef_pos`.}
#'   \item{`coef_names`}{length `p`; the labels of those rows.}
#'   \item{`extra_var`}{length `n`; variance that is NOT coefficient
#'     uncertainty, kept separate rather than folded into `A V A'`
#'     because an exact [gp()]'s kriging variance and a new grouping
#'     level's marginal variance are not.}
#'   \item{`nonest`}{length `n`; rows that load on a direction the
#'     rank-deficient design could not identify.}
#' }
#' `var(eta)` is `rowSums((A %*% V) * A) + extra_var`, and
#' `predict(se.fit = TRUE)` is written that way.
#' @section A nonlinear body:
#' For a nonlinear predictor `A` is a JACOBIAN rather than a design, and
#' it is computed by taping the body against the coefficients it reaches
#' through. That includes a [ps()] block, whose coefficients enter the
#' body through a spline evaluated at an argument the parameters move.
#' `predict(se.fit = TRUE)` stays refused for a nonlinear predictor;
#' this is the route.
#'
#' The Jacobian is exact, and the delta method built on it is still a
#' first-order approximation, which for a warped curve is a stronger
#' assumption than it is for a linear one. `allow_new_levels = TRUE` is
#' refused there, and so is a contributing exact [gp()], because neither
#' variance has a chain rule through the body that has been measured.
#' @section A reduced-rank block:
#' An `rr()` block's loadings live in `theta`, so a design over
#' `(beta, b)` alone is INCOMPLETE. `A` carries the loading columns
#' too, through `rr_jacobians()`, and `coef_pos` names their `theta`
#' rows; a caller therefore gets the whole delta method rather than
#' discovering a missing piece.
#' @seealso [frm_joint_cov()] for the covariance alone,
#'   [frmtmb-extension-api]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(120), g = factor(rep(1:12, each = 10)))
#' dd$y <- rnorm(120, 1 + 2 * dd$x + rnorm(12, 0, 0.5)[dd$g], 0.4)
#' fit <- frm(bf(y ~ x + (1 | g)), data = dd)
#' nd <- data.frame(x = c(-1, 0, 1), g = factor(1, levels = levels(dd$g)))
#' lb <- frm_lp_basis(fit, newdata = nd, re.form = NA)
#' str(lb$A)
#' lb$coef_names
#'
#' # the covariance of the WHOLE grid, which predict() reduces to its
#' # diagonal
#' Sigma <- lb$A %*% lb$V %*% t(lb$A)
#' all.equal(sqrt(diag(Sigma)),
#'           predict(fit, newdata = nd, re.form = NA, se.fit = TRUE)$se.fit)
#' @export
frm_lp_basis <- function(object, newdata = NULL, dpar = NULL, resp = NULL,
                         re.form = NULL, allow_new_levels = FALSE) {
  require_frmtmb_fit(object, "frm_lp_basis()")
  require_fitted(object, "frm_lp_basis()")
  check_flag(allow_new_levels, "allow_new_levels")
  if (!is.null(newdata) && !is.data.frame(newdata)) {
    stop("`newdata` must be a data frame, or NULL to use the training ",
         "data, not ", arg_desc(newdata), call. = FALSE)
  }
  check_re_form(re.form)
  use_re <- re_form_keeps(re.form)
  resp <- resp %||% names(object$spec$responses)[1]
  rspec <- object$spec$responses[[resp]]
  if (is.null(rspec)) {
    stop_unknown_response(object, resp)
  }
  dpar <- dpar %||% if ("mu" %in% names(rspec$dpars)) "mu" else
    rspec$primary_dpars[1]
  lp <- object$frame[["linpreds"]][[linpred_key(resp, dpar)]]
  if (is.null(lp)) {
    stop("frm_lp_basis(): unknown dpar '", dpar, "' for response '",
         resp, "'. Available: ",
         paste(names(rspec$dpars), collapse = ", "), call. = FALSE)
  }
  jc <- get_joint_cov(object)
  if (!is.null(lp[["nl_body"]])) {
    return(lp_basis_nl(object, lp, rspec, newdata, use_re,
                       allow_new_levels, jc))
  }
  ed <- lp_eta_design(object, lp, newdata, use_re, allow_new_levels)
  has_rr <- isTRUE(object$frame[["has_rr"]])
  rrj <- if (has_rr) rr_jacobians(object)
  da <- lp_delta_A(object, lp, ed, newdata, use_re, jc, has_rr, rrj)
  lp_basis_out(object, jc, ed[["eta"]], as.matrix(da$A), da$coef_pos,
               lp_extra_var_vec(object, ed, use_re), ed[["nonest"]])
}

#' The variance sources that are not coefficient uncertainty, summed to
#' one vector per row. `predict()` adds them in place; a caller that
#' gets `A` and `V` separately needs them separately too.
#'
#' @noRd
lp_extra_var_vec <- function(object, ed, use_re) {
  n <- ed[["n"]]
  out <- numeric(n)
  ev <- lp_extra_var(object, ed, use_re)
  for (B in extra_var_blocks(ev$new_levels, n)) {
    Mr <- B$M[B$rows, , drop = FALSE]
    out[B$rows] <- out[B$rows] + rowSums((Mr %*% B$S) * Mr)
  }
  for (gv in ev$gp) out <- out + gv
  out
}

#' Assemble the return value once, so the linear and nonlinear branches
#' cannot drift into different shapes.
#'
#' @noRd
lp_basis_out <- function(object, jc, eta, A, coef_pos, extra_var, nonest) {
  lab <- joint_coef_labels(object, jc)
  list(eta = eta, A = A, coef_pos = coef_pos,
       V = jc$V[coef_pos, coef_pos, drop = FALSE],
       coef_names = if (is.null(lab)) NULL else lab[coef_pos],
       extra_var = extra_var,
       nonest = nonest %||% rep(FALSE, length(eta)))
}

#' Which component of the parameter list each row of the joint
#' covariance belongs to, and its index inside that component.
#'
#' @noRd
joint_pos_map <- function(jc) {
  rn <- jc$names
  idx <- integer(length(rn))
  for (cp in unique(rn)) {
    i <- which(rn == cp)
    idx[i] <- seq_along(i)
  }
  list(comp = rn, idx = idx)
}

#' `d eta / d coef` for a NONLINEAR predictor, by taping the body.
#'
#' `eta` is affine in the coefficients for every LINEAR dpar the body
#' names, exactly: `eta(c) = eta(chat) + A (c - chat)`, with `A` the
#' design `lp_delta_A()` already builds. So the whole body is a
#' composition of affine maps, link inverses and whatever R code the
#' body contains, and taping that composition against the coefficient
#' subvector gives the Jacobian with no perturbation and no finite
#' difference.
#'
#' A `ps()` block joins through the same tape: its closures are rebuilt
#' from the perturbed parameter list on every tape evaluation, which is
#' why they live in `ev` and not in `nl_env`.
#'
#' @noRd
lp_basis_nl <- function(object, lp, rspec, newdata, use_re,
                        allow_new_levels, jc) {
  if (isTRUE(allow_new_levels)) {
    stop("frm_lp_basis(): allow_new_levels = TRUE is refused for a ",
         "nonlinear predictor. A new level contributes its block's ",
         "marginal variance, and how that variance passes through a ",
         "nonlinear body has not been measured", call. = FALSE)
  }
  pm <- joint_pos_map(jc)
  est <- object$estimates
  parts <- list()          # per contributing dpar: eta_hat, A, coef_pos
  build <- function(nm) {
    if (!is.null(parts[[nm]])) return(invisible(NULL))
    lpk <- object$frame[["linpreds"]][[linpred_key(rspec$resp_name, nm)]]
    if (is.null(lpk)) {
      stop("frm_lp_basis(): the nonlinear body names '", nm,
           "', which is not a parameter of response '",
           rspec$resp_name, "'", call. = FALSE)
    }
    if (!is.null(lpk[["nl_body"]])) {
      for (sub in c(lpk[["nl_pars"]], lpk[["nl_dpar_refs"]])) build(sub)
      parts[[nm]] <<- list(lp = lpk, nl = TRUE)
      return(invisible(NULL))
    }
    edk <- lp_eta_design(object, lpk, newdata, use_re, FALSE)
    ex <- lp_extra_var_vec(object, edk, use_re)
    if (any(ex != 0)) {
      stop("frm_lp_basis(): the nonlinear body reaches '", nm,
           "', which contributes variance that is not coefficient ",
           "uncertainty (an exact gp() kriging variance). Its chain ",
           "rule through a nonlinear body has not been measured",
           call. = FALSE)
    }
    has_rr <- isTRUE(object$frame[["has_rr"]])
    dak <- lp_delta_A(object, lpk, edk, newdata, use_re, jc, has_rr,
                      if (has_rr) rr_jacobians(object))
    parts[[nm]] <<- list(lp = lpk, nl = FALSE, eta = unname(edk[["eta"]]),
                         A = as.matrix(dak$A), pos = dak$coef_pos,
                         nonest = edk[["nonest"]])
    invisible(NULL)
  }
  for (nm in c(lp[["nl_pars"]], lp[["nl_dpar_refs"]])) build(nm)

  pos <- integer(0)
  for (p in parts) if (!isTRUE(p$nl)) pos <- c(pos, p$pos)
  # a ps() block reaches the body through its own coefficients rather
  # than through any design, so its rows are added by name
  for (lpk in c(list(lp), lapply(parts, `[[`, "lp"))) {
    for (pt in lpk[["ps_terms"]] %||% list()) {
      pos <- c(pos, which(pm$comp == pt[["par"]])[pt[["beta_idx"]]],
               which(pm$comp == "b")[pt[["b_idx"]]])
    }
  }
  pos <- sort(unique(pos))
  if (!length(pos)) {
    stop("frm_lp_basis(): the nonlinear predictor '", lp[["dpar"]],
         "' reaches no estimated coefficient, so it has no design",
         call. = FALSE)
  }
  chat <- vapply(pos, function(k) est[[pm$comp[k]]][pm$idx[k]], 0)

  dl <- lp_basis_nl_data(object, newdata)
  eta_fun <- function(cc) {
    "[<-" <- RTMB::ADoverload("[<-")
    pars <- est
    for (k in seq_along(pos)) {
      pars[[pm$comp[pos[k]]]][pm$idx[pos[k]]] <- cc[k]
    }
    bvec <- expand_b(object$frame, pars[["b"]], pars[["theta"]])
    vals <- list()
    ev_one <- function(lpk, nm) {
      p <- parts[[nm]]
      if (isTRUE(p$nl)) {
        e <- eval(lpk[["nl_body"]],
                  c(vals[c(lpk[["nl_pars"]], lpk[["nl_dpar_refs"]])],
                    dl[[nm]], ps_env(lpk, pars, bvec)),
                  ad_overload_env(lpk[["nl_env"]], lpk[["nl_body"]]))
      } else {
        e <- p$eta
        for (k in seq_along(p$pos)) {
          j <- match(p$pos[k], pos)
          e <- e + p$A[, k] * (cc[j] - chat[j])
        }
      }
      e
    }
    for (nm in names(parts)) {
      p <- parts[[nm]]
      e <- ev_one(p$lp, nm)
      # the objective stores every dpar as linkinv(eta) and reads the
      # body's names out of that store, so this does the same
      vals[[nm]] <- p$lp[["link"]]$linkinv(e)
    }
    eval(lp[["nl_body"]],
         c(vals[c(lp[["nl_pars"]], lp[["nl_dpar_refs"]])],
           dl[[lp[["dpar"]]]], ps_env(lp, pars, bvec)),
         ad_overload_env(lp[["nl_env"]], lp[["nl_body"]]))
  }
  tp <- RTMB::MakeTape(eta_fun, chat)
  eta <- as.numeric(tp(chat))
  A <- tp$jacobian(chat)
  nonest <- rep(FALSE, length(eta))
  for (p in parts) {
    if (!isTRUE(p$nl) && !is.null(p$nonest)) nonest <- nonest | p$nonest
  }
  lp_basis_out(object, jc, eta, A, pos, numeric(length(eta)), nonest)
}

#' The raw data columns each nonlinear body reads, at `newdata` or in
#' sample, one list per dpar that has a body.
#'
#' @noRd
lp_basis_nl_data <- function(object, newdata) {
  out <- list()
  for (lpk in object$frame[["linpreds"]]) {
    if (is.null(lpk[["nl_body"]])) next
    dl <- lpk[["data_list"]]
    if (!is.null(newdata)) {
      dl <- lapply(stats::setNames(names(dl), names(dl)), function(v) {
        if (is.null(newdata[[v]])) {
          stop_newdata_missing(v)
        }
        newdata[[v]]
      })
    }
    out[[lpk[["dpar"]]]] <- dl
  }
  out
}
