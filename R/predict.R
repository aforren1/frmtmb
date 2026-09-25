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
  Q <- joint_precision(fit)
  if (is.null(Q)) {
    V <- sdr_of(fit)$cov.fixed
    rn <- rownames(V)
  } else {
    # same degradation vcov() uses: a singular joint precision gives NaN
    # standard errors and one warning naming diagnose(), not a raw
    # LAPACK message from deep inside predict()
    V <- as.matrix(solve_joint_precision(Q, cache, fit))
    rn <- rownames(Q)
  }
  cache$Vjoint <- list(V = V, names = rn)
  cache$Vjoint
}

#' The joint precision of every estimated parameter, fixed and random,
#' or `NULL` when the objective has no inner parameters. Memoized in
#' `fit$cache`: `get_joint_cov()` inverts it, and `predict()` draws the
#' group effects from its conditional blocks.
#'
#' @noRd
joint_precision <- function(fit) {
  # the stored sdreport's own matrix wins whenever it has one, so a
  # caller that replaces the sdreport is read, not a copy from before;
  # only the separate joint-precision sdreport an ML fit needs is kept
  Q <- sdr_of(fit)$jointPrecision
  if (!is.null(Q)) return(Q)
  cache <- fit$cache
  if (!is.null(cache$Qjoint)) return(cache$Qjoint)
  if (!is.null(sdr_of(fit)$par.random) && length(sdr_of(fit)$par.random)) {
    Q <- autoscale_sdreport(fit, jp = TRUE)$jointPrecision
  }
  if (is.environment(cache)) cache$Qjoint <- Q
  Q
}

#' The three argument checks `predict()` and `frm_lp_basis()` share.
#'
#' They take the same `newdata`, `re_formula` and `resp`, and two functions
#' describing one argument two ways is how a user learns that the second
#' one is a different argument. One template each, which is also the
#' property `test-message-uniqueness.R` asserts.
#'
#' @noRd
check_re_form <- function(re_formula) {
  if (!is.null(re_formula) && !inherits(re_formula, "formula") &&
        !(length(re_formula) == 1L && is.na(re_formula))) {
    frm_stop(
      "`re_formula` must be NULL to keep every random effect, NA to drop ",
      "them all, or a one-sided formula naming the ones to keep, not ",
      arg_desc(re_formula), call. = FALSE)
  }
  invisible(NULL)
}

#' Whether `re_formula` keeps the random effects.
#'
#' `NULL` keeps them, `NA` drops them, and a one-sided formula keeps
#' them when it names at least one group-level term: `~0` and `~1` name
#' none, and brms reads both as "no group-level effects". WHICH terms a
#' formula keeps is `re_resolve()`'s question; this answers only whether
#' any are. Factored out beside the message templates for
#' the same reason those were: `predict()` and `frm_lp_basis()` take the
#' same argument, and the whole point of section 3 of the review is that
#' the two return the same numbers, which they can only do while they
#' read the argument the same way.
#'
#' @noRd
re_form_keeps <- function(re_formula) {
  if (is.null(re_formula)) return(TRUE)
  if (!inherits(re_formula, "formula")) return(FALSE)
  length(re_formula_bars(re_formula)) > 0L
}

#' @noRd
stop_unknown_response <- function(object, resp) {
  frm_stop("Unknown response: '", resp, "'. Available: ",
           paste(names(object$spec$responses), collapse = ", "), call. = FALSE)
}

#' @noRd
stop_newdata_missing <- function(v) {
  frm_stop("Variable '", v, "' missing from newdata", call. = FALSE)
}

#' Refuse, by name, the two newdata faults model.frame() would otherwise
#' report in its own words: a variable the design needs that neither
#' `newdata` nor the formula environment holds, and a factor level the
#' fit never saw. brms refuses both as a brms_error.
#'
#' The lookup is model.frame()'s own, `newdata` first and then the
#' environment of `tt`, so a variable it would find is never refused.
#'
#' @noRd
check_newdata_frame <- function(tt, newdata, xlev) {
  env <- environment(tt) %||% globalenv()
  nms <- names(newdata)
  for (v in setdiff(all.vars(tt), nms)) {
    if (!exists(v, envir = env)) stop_newdata_missing(v)
  }
  for (v in intersect(names(xlev), nms)) {
    x <- newdata[[v]]
    if (!is.factor(x) && !is.character(x)) next
    lev <- unique(as.character(x[!is.na(x)]))
    new <- lev[!lev %in% xlev[[v]]]
    if (length(new)) {
      frm_stop("newdata has ", if (length(new) == 1L) "a level" else
                 "levels", " of `", v, "` that the fit did not see: ",
               paste0("'", new, "'", collapse = ", "), ". The fitted ",
               "levels are ", paste0("'", xlev[[v]], "'", collapse = ", "),
               ". A population-level factor cannot take a new level",
               call. = FALSE)
    }
  }
  invisible(NULL)
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
      frm_stop("mo(): new data contain unknown categories", call. = FALSE)
    }
    return(as.integer(v) - 1L)
  }
  codes <- as.integer(round(v))   # grids may land between categories
  if (any(codes < 0 | codes > mi$D)) {
    frm_stop("mo(): new data outside the fitted range 0..", mi$D,
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
#' smooths, which is what `re_formula = NA` drops.
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
#' on the group, so `re_formula = NA` is a way out of it; any other missing
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
      frm_stop("predict(newdata = ) for the factor-smooth term ", si$label,
               " needs the grouping column `", gv, "`: the term holds one ",
               "curve per level of `", gv, "`, so a prediction conditional ",
               "on it has to say which level each row belongs to. Add the ",
               "column to newdata, or ask for the population curve with ",
               "re_formula = NA, which drops the term and needs no level",
               call. = FALSE)
    }
    frm_stop("predict(newdata = ) for the smooth term ", si$label,
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
      frm_stop("New levels in the factor-smooth term ", si$label, ": ",
               paste(new, collapse = ", "), ". The term has no curve for ",
               "them. Use allow_new_levels = TRUE to predict them at the ",
               "population level, or re_formula = NA for the population curve ",
               "at every row", call. = FALSE)
    }
    if (length(new) && is.null(si$sm$flev)) {
      # an fs basis zero-rows an unknown level; a factor bs = "re"
      # basis has one design column per fitted level and nothing else,
      # so there is no population row to hand back
      frm_stop("allow_new_levels = TRUE cannot predict the new level(s) ",
               paste(new, collapse = ", "), " of the bs = \"re\" smooth ",
               "term ", si$label, ": its design has one column per fitted ",
               "level and no zero row for a new one. Use re_formula = NA for ",
               "the population curve, which is what a new level would ",
               "receive anyway", call. = FALSE)
    }
  }
  invisible(NULL)
}

#' The grouping columns one linear predictor's random effects need, and
#' what to do when `newdata` does not carry them.
#'
#' brms's rule, from `validate_newdata()`: "grouping factors do not need
#' to be specified by the user if new levels are allowed". It fills the
#' absent ones with `NA`, so every row is a level the fit never saw and
#' the block contributes its population value and, in an interval, its
#' marginal variance. frmtmb reached base R's
#' `eval(comp$bar[[3]], newdata, env)` instead and stopped at "object 'g'
#' not found", which names neither the argument nor the rule
#' (dev/adefects-findings.md, D6).
#'
#' Without `allow_new_levels` brms stops too, and on the same base R
#' error; measured on `brmsfit_example1` with `visit` removed. The
#' refusal here says which column and which argument, and is classed,
#' which item 2.6e asks of every refusal frmtmb makes.
#'
#' Only the ordinary grouping blocks are read. A factor smooth has its
#' own named refusal in `smooth_newdata_check()`, and an `spde()` block's
#' levels are mesh row numbers, where a new level means nothing; both
#' keep the behavior they had.
#'
#' @noRd
fill_new_group_vars <- function(fit, lp, newdata, allow_new_levels, env) {
  key <- linpred_key(lp[["resp"]], lp[["dpar"]])
  need <- character(0)
  for (bk in fit$frame[["re_blocks"]]) {
    if (bk[["covstruct"]] %in% c("smooth", "gp", "hsgp", "spde")) next
    for (comp in bk[["components"]]) {
      if (comp$lp_key != key) next
      need <- c(need, if (is.null(comp$mm)) all.vars(comp$bar[[3L]]) else
                        comp$mm$gvars)
    }
  }
  # the environment lookup is model.frame()'s own, the one
  # check_newdata_frame() honors, so a variable R would have found is
  # neither refused nor overwritten with NA here
  need <- setdiff(unique(need), names(newdata))
  need <- need[!vapply(need, exists, NA, envir = env)]
  if (!length(need)) return(newdata)
  if (!allow_new_levels) {
    frm_stop("newdata has no column ",
             paste0("`", need, "`", collapse = ", "),
             ", the grouping factor", if (length(need) > 1L) "s" else "",
             " of a random effect. Add ",
             if (length(need) > 1L) "them" else "it",
             ", or say what the prediction should do without ",
             if (length(need) > 1L) "them" else "it",
             ": allow_new_levels = TRUE treats every row as an unseen ",
             "level (the population value, plus that block's variance in ",
             "an interval), and re_formula = NA drops the random effects ",
             "altogether", call. = FALSE)
  }
  for (v in need) newdata[[v]] <- NA
  newdata
}

#' Rebuild the design pieces of one linear predictor for new data:
#' dense X (parametric + smooth null-space columns), per-block RE
#' component designs with level indices, and the offset.
#'
#' @noRd
pred_design <- function(fit, lp, newdata, allow_new_levels = FALSE,
                        use_re = TRUE) {
  env <- fit$spec$responses[[lp[["resp"]]]]$formula_env
  if (use_re) {
    newdata <- fill_new_group_vars(fit, lp, newdata, allow_new_levels, env)
  }
  tt <- patch_predvars(lp[["terms"]], fit$frame[["predvar_map"]])
  check_newdata_frame(tt, newdata, xlev_for(lp[["xlevels"]], tt))
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
        frm_stop("predict(newdata = ) is not supported for the smooth ",
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
      frm_stop("Interaction multiplier '", deparse1(mult_expr),
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
      frm_stop("mi(", mt$var, "): newdata must supply complete values",
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
      check_newdata_frame(tt2, newdata, xlev_for(lp[["xlevels"]], tt2))
      mf2 <- stats::model.frame(tt2, newdata, na.action = stats::na.pass,
                                xlev = xlev_for(lp[["xlevels"]], tt2))
      mm <- stats::model.matrix(tt2, mf2)
      if (!identical(colnames(mm), comp$cnms)) {
        frm_stop("Random-effect design for `", comp$label, "` does not match ",
                 "the fitted model (columns: ",
                 paste(colnames(mm), collapse = ", "), " vs ",
                 paste(comp$cnms, collapse = ", "), ")", call. = FALSE)
      }
      # a partial re_formula keeps some columns of this term (re_view())
      if (!is.null(comp$keep_cols)) mm[, !comp$keep_cols] <- 0
      gvr <- eval(comp$bar[[3]], newdata, env)
      # an spde block's levels are mesh ROW NUMBERS, so the node has to
      # be read as a number here too: as.character() on a double would
      # spell node 100000 as "1e+05" and lose the column
      gv <- if (bk[["covstruct"]] == "spde") spde_node_labels(gvr) else {
        as.character(gvr)
      }
      j <- match(gv, bk[["levels"]])
      if (anyNA(j) && !allow_new_levels) {
        frm_stop("New levels in grouping factor `", deparse1(comp$bar[[3]]),
                 "`: ", paste(unique(gv[is.na(j)]), collapse = ", "),
                 ". Use allow_new_levels = TRUE to predict them at the ",
                 "population level", call. = FALSE)
      }
      # new_key is the level label: two rows at the SAME unseen level
      # load one draw of its effect and two DIFFERENT unseen levels load
      # independent ones. Without it every unseen level of a block
      # shared one key, which a per-row variance cannot see and a joint
      # draw can: predict(summary = FALSE) gave three distinct new
      # groups the same effect in every replicate
      re_parts[[length(re_parts) + 1L]] <- list(bk = bk, comp = comp,
                                                mm = mm, j = j,
                                                new_key = gv)
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
    frm_stop("Multi-membership design for `", comp$label, "` does not match ",
             "the fitted model (columns: ",
             paste(md$cnms, collapse = ", "), " vs ",
             paste(comp$cnms, collapse = ", "), ")", call. = FALSE)
  }
  gv <- mm_member_values(mms, newdata, env)
  if (anyNA(iw$J) && !allow_new_levels) {
    new <- unique(unlist(lapply(gv, function(v) {
      setdiff(as.character(v), bk[["levels"]])
    }), use.names = FALSE))
    frm_stop("New levels in multi-membership factor `", mms$label, "`: ",
             paste(new, collapse = ", "),
             ". Use allow_new_levels = TRUE to predict those memberships ",
             "at the population level; the row's remaining members still ",
             "contribute their fitted effects", call. = FALSE)
  }
  lapply(seq_len(iw$n_members), function(k) {
    mmk <- md$designs[[k]] * iw$W[, k]
    # a partial re_formula keeps some columns of this term (re_view())
    if (!is.null(comp$keep_cols)) mmk[, !comp$keep_cols] <- 0
    # new_key names WHICH unseen level this member landed on, because a
    # row whose members carry the SAME unseen label loads one draw of
    # the block, not two independent ones (see extra_var_blocks())
    list(bk = bk, comp = comp,
         mm = mmk, j = iw$J[, k],
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

#' Jacobians of the coefficient-space expansion: `d cvec/d b` (sparse)
#' and `d cvec/d theta` for the rr loading parameters (finite
#' differences on `expand_b`).
#'
#' `d cvec/d b` is the identity for most blocks, but NOT for the two
#' whose coefficients are a function of the parameters rather than the
#' parameters themselves. An `rr` block expands through its loadings.
#' An `esicar` block is centered per connected component, so its entry
#' is the projection `P` that `car_center_jacobian()` builds; pairing
#' `Z` with `b` through the identity there put the inert component
#' means into every standard error, exactly `con_sd^2` of variance that
#' the field the predictor sees does not have.
#'
#' The name is historical: rr was the first block that needed this.
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
    } else if (block_is_esicar(bk)) {
      # d(P b)/d b = P, per connected component. The same predicate
      # expand_b() branches on, so the Jacobian cannot describe a
      # different expansion than the one it differentiates.
      pj <- car_center_jacobian(bk[["aux_car"]])
      ii <- c(ii, bk[["c_idx"]][pj$i])
      jj <- c(jj, bk[["b_idx"]][pj$j])
      xx <- c(xx, pj$x)
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
    frm_stop(what, " is not supported yet for multivariate fits",
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
        frm_stop("Addition term ", label, " could not be evaluated on ",
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
      frm_warning("Addition term ", label, " could not be evaluated on ",
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
#' Not the same thing as `frm_linpred(type = "response", dpar = )` any more:
#' a family may declare a REPORTING scale for a dpar (a mixture's
#' mixing weights are a softmax over the component predictors, and the
#' predictor itself is not a probability), and feeding a reported value
#' back to `lpdf`, `sim` or `mean_fn` would apply the transform twice.
#' Every consumer that hands dpar values to the family goes through
#' here; only the user-facing `predict()` surface reports.
#'
#' @noRd
dpars_natural <- function(fit, rspec, newdata, re_formula,
                          allow_new_levels = FALSE) {
  rn <- rspec$resp_name
  dp <- list()
  for (dnm in names(rspec$dpars)) {
    lp <- fit$frame[["linpreds"]][[linpred_key(rn, dnm)]]
    eta <- frm_linpred(fit, newdata = newdata, dpar = dnm, resp = rn,
                       re_formula = re_formula, type = "link",
                       allow_new_levels = allow_new_levels)
    dp[[dnm]] <- as.vector(lp[["link"]]$linkinv(eta))
  }
  dp
}

#' `y | se(s)` without `sigma = TRUE`: the residual standard deviation
#' beyond the known `s` is not in the density at all - `resid_sd()`
#' returns `s` alone and the dpar is mapped out at the link-scale zero.
#' `sigma()` has always reported that as 0; `frm_linpred(dpar = "sigma")`
#' reported the log link's inverse of the mapped-out coefficient, 1,
#' which reads as an estimate of a parameter the model does not have.
#'
#' The dpar this applies to is the one the FAMILY says `se()` replaces,
#' so a family that declares `se_dpar = "tau"` reports `tau` the same
#' way, and one that declares `se_dpar = NA` has no such dpar.
#'
#' @noRd
se_unused_sigma <- function(rspec, dpar) {
  repl <- family_se_dpar(rspec$family)
  !is.null(repl) && !is.na(repl) && identical(dpar, repl) &&
    !is.null(rspec$aterms[["se"]]) &&
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
#' values (`fitted()`'s convention, extended to newdata and `re_formula`).
#'
#' @noRd
predict_mean_response <- function(fit, rspec, newdata, re_formula,
                                  allow_new_levels) {
  fam <- rspec$family
  rn <- rspec$resp_name
  if (is.null(newdata) && is.null(re_formula)) {
    # exactly fitted(): dpars at the estimates, conditional on the modes
    dp <- eval_dpars(fit)[[rn]]
    out <- response_mean(fam, dp, fit$frame[["aterm_values"]][[rn]])
    return(napred(fit, out))
  }
  dp <- dpars_natural(fit, rspec, newdata, re_formula, allow_new_levels)
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

#' Arguments `predict()` used to take, or that brms takes and a point
#' estimate cannot answer.
#'
#' The two retired spellings are named here rather than left to the
#' generic unknown-argument message because both were live in 0.57.0:
#' one of them silently changed the answer, and the other silently did
#' not.
#'
#' @noRd
predict_retired <- c(
  re.form = paste("lme4's spelling, no longer accepted. brms is the",
                  "tiebreaker on a name, so this setting is",
                  "`re_formula` here as it is in brms. Pass",
                  "re_formula ="),
  allow.new.levels = paste("lme4's spelling, no longer accepted. brms",
                           "spells it `allow_new_levels`, and so does",
                           "this. Pass allow_new_levels ="),
  ndraws = paste("frm_linpred() evaluates one parameter vector, so",
                 "there is nothing to thin. predict() SIMULATES and",
                 "takes ndraws; on posterior draws,",
                 "frmtmb.sample's posterior_linpred() takes it too"),
  draw_ids = paste("a maximum likelihood fit has no draws to index.",
                   "frmtmb.sample's posterior_linpred() takes draw_ids"),
  sort = "there is no draws dimension here to sort against",
  summary = paste("the return value is one number per row and there is",
                  "nothing behind it to summarize. fitted() is the",
                  "expected response in brms's summary shape"),
  robust = "a median over draws needs draws",
  probs = paste("quantiles of a linear predictor need draws. Use",
                "se.fit = TRUE here, or fitted() for brms's",
                "quantile columns"),
  nlpar = paste("a non-linear parameter is reached through `dpar` in",
                "this package"),
  transform = paste("the scale is chosen with `type` here:",
                    "type = \"response\" is brms's default scale and",
                    "type = \"link\" is this one's")
)

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
#' `frm_linpred(type = "response") == fitted()` identity holds here too.
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
#' @param re_formula Which group-level terms enter the prediction.
#'   `NULL` (default) keeps all of them; `NA` keeps none, which is the
#'   population-level prediction. A one-sided formula keeps the terms it
#'   names; see *A one-sided `re_formula`*. See
#'   *What `re_formula = NA` drops* for what that means when the model has
#'   smooths.
#' @section A one-sided `re_formula`:
#' A formula keeps the group-level terms it names and drops the others,
#' with brms's rule (`update_re_terms()`). A term in the formula is
#' matched to a term of the fit that has the same grouping factor and
#' whose columns include the formula term's columns. So on a fit with
#' `(1 + x | g) + (1 | h)`:
#'
#' * `~ (1 | g)` keeps the intercept of `g` and drops its slope and `h`.
#' * `~ (0 + x | g)` keeps the slope of `g` alone.
#' * `~ (1 + x | g) + (1 | h)` keeps everything, and the answer is
#'   identical to `re_formula = NULL`.
#' * `~0` and `~1` name no group-level term, so they are `NA`.
#'
#' The same term is kept in every distributional parameter that has it.
#' brms's spellings are read too: an id (`(1 | p | g)`), a `gr()`
#' wrapper, `||` and a nested group (`a/b`). A term outside the bars is
#' ignored.
#'
#' A dropped term is dropped everywhere: from the estimate, from the
#' standard error, and from `predict()`'s draws. Its grouping column is
#' not needed in `newdata`, and a level of it the fit never saw is not
#' an error.
#'
#' Two cases are refused. A formula term that matches no term of the
#' fit is an error that names it, because a misspelled grouping factor
#' would otherwise change the answer with nothing said; brms drops such
#' a term silently. And a formula that keeps SOME terms is refused on a
#' fit that also has group-level content a formula cannot name, such as
#' a factor-smooth term: `NA` drops that content and `NULL` keeps it,
#' and a partial formula cannot say which.
#' @section What `re_formula = NA` drops:
#' `re_formula = NA` (equivalently `~0`) asks for the POPULATION-level
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
#' grouping column may be left out entirely when `re_formula = NA`. It is
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
#' The estimate is at the random-effect modes, and the standard error
#' carries their uncertainty, through the joint covariance, as `se.fit`
#' does for the linear predictor. Unseen grouping levels
#' (`allow_new_levels = TRUE`) add their block's marginal variance,
#' propagated through the same gradients.
#' @param allow_new_levels Predict unseen grouping-factor levels at the
#'   population level instead of erroring. A factor-smooth term
#'   (`bs = "fs"`) follows the same rule: a level it never saw
#'   contributes nothing, which leaves the population curve.
#'
#'   It also makes the grouping COLUMN optional, as it does in brms
#'   (`validate_newdata()`: "grouping factors do not need to be
#'   specified by the user if new levels are allowed"). A column
#'   `newdata` does not carry is filled with `NA`, so every row is an
#'   unseen level. Without `allow_new_levels` a missing grouping column
#'   is refused by name, and the refusal offers this argument and
#'   `re_formula = NA`, which drops the random effects instead.
#' @param ... Refused. An argument this method does not have is an
#'   error naming it, and the two lme4 spellings that were live in
#'   0.57.0 (`re.form`, `allow.new.levels`) are refused by name with
#'   the brms spelling that replaced them.
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
#'   input. `frm_linpred()` and `model.frame()` carry the row names of
#'   the data they were computed from, and
#'   `tests/testthat/test-edgecases.R` asserts it, including across a
#'   row dropped for missingness. `fitted()`, `residuals()` and
#'   `predict()` are brms's summary matrices and leave their row
#'   dimnames NULL as brms does, because brms names those rows after
#'   the draws matrix's columns, which are unnamed; the same test
#'   asserts THAT, so neither convention can drift. `vcov()`,
#'   `confint()` and `fixef()` share one coefficient naming scheme, and
#'   so does `coef()`;
#'   `tests/testthat/test-methods-audit.R` asserts that
#'   `vcov(full = TRUE)` carries exactly the row names of `confint()`
#'   and that `vcov()` is its population-level sub-block, matched
#'   through `brms_coef_table()` because `vcov()` names its rows as
#'   brms does. The stored model frame keeps the input row names.
#' @srrstats {RE4.9} Modelled values of the response are returned by
#'   `fitted()`, and by `frm_linpred(type = "response")`, which is asserted
#'   to equal `fitted()` on the training data wherever both are defined:
#'   the fuzz harness checks the invariant across the family grid, and
#'   `tests/testthat/test-methods-audit.R` checks it on the cases where
#'   the two could plausibly diverge (zero-inflated, `trials()`
#'   binomial). A family with no mean on the response scale, such as
#'   `cox()`, refuses both alike. On an ordinal family the modelled
#'   response is a category
#'   distribution rather than a mean, so both return the same `n x K`
#'   matrix of category probabilities (the brms convention), which
#'   `fitted()` names `P(Y = k)` as brms does and `frm_linpred()` names
#'   by the response's own levels; the latent linear predictor stays
#'   reachable as `frm_linpred(type = "link")`.
#' @srrstats {RE4.14} Uncertainty is available away from the observed
#'   data, and in three shapes: `fitted()` reports brms's `Est.Error`
#'   and `Q` columns around the expected response,
#'   [predict.frmtmb_fit()] simulates the PREDICTIVE distribution, which
#'   is much wider because it carries the observation noise, and
#'   `se.fit = TRUE` here returns delta-method standard errors that
#'   include fixed-effect and random-effect uncertainty; unseen grouping
#'   levels add their block's marginal variance, and exact `gp()` terms
#'   add the Gaussian-process conditional (kriging) variance, so the
#'   reported error grows with distance from the observed positions.
#' @srrstats {RE4.16} New groups can be submitted here, and to
#'   `fitted()` and `predict()`. Levels
#'   of a grouping factor that were not in the training data error by
#'   default, naming the offending levels, and are predicted at the
#'   population level under `allow_new_levels = TRUE`, which is brms's
#'   spelling and the only one: lme4's `allow.new.levels` is refused
#'   and the refusal names the replacement.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
#' fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
#'
#' # the link scale by default; "response" is what fitted() estimates
#' head(frm_linpred(fit))
#' max(abs(frm_linpred(fit, type = "response") -
#'           fitted(fit)[, "Estimate"]))
#'
#' # re_formula = NA drops the random effects: the population prediction
#' nd <- data.frame(x = c(-1, 0, 1), g = factor(1, levels = levels(dd$g)))
#' frm_linpred(fit, newdata = nd, re_formula = NA, type = "response")
#'
#' # delta-method standard errors, on whichever scale was asked for
#' p <- frm_linpred(fit, newdata = nd, se.fit = TRUE)
#' cbind(fit = p$fit, se = p$se.fit)
#'
#' # a level the fit never saw errors unless it is allowed explicitly,
#' # in which case it is predicted at the population level
#' nd_new <- data.frame(x = 0, g = factor("new"))
#' try(frm_linpred(fit, newdata = nd_new))
#' frm_linpred(fit, newdata = nd_new, allow_new_levels = TRUE)
#'
#' # a distributional parameter instead of the mean
#' fit2 <- frm(bf(y ~ x, sigma ~ x) + gaussian(), data = dd)
#' head(frm_linpred(fit2, dpar = "sigma", type = "response"))
#' @seealso [fitted.frmtmb_fit()] for the same expected response in
#'   brms's four-column shape, [predict.frmtmb_fit()] for the predictive
#'   summary, and [frmtmb-scales], which states which scale every method
#'   reports. The default here is the LINK scale.
#' @export
frm_linpred <- function(object, newdata = NULL,
                        type = c("link", "response", "conditional",
                                 "zprob", "zlink", "disp"),
                        dpar = NULL, resp = NULL, re_formula = NULL,
                        se.fit = FALSE, allow_new_levels = FALSE, ...) {
  # A warning here was not enough: `predict(fit, re_form = NA)` warned
  # and then returned the CONDITIONAL prediction, and a warning in a
  # loop or under suppressWarnings() is a wrong number with no record.
  frm_check_dots(..., .unsupported = predict_retired)
  require_fitted(object, "frm_linpred()")
  # re_formula is read three lines down as "NULL, a formula, or anything
  # else drops the random effects", so a typo used to return the
  # POPULATION prediction and say nothing. newdata that is not
  # rectangular reached the design builder and failed there on
  # "non-conformable arrays", which names neither the argument nor the
  # fault.
  check_flag(se.fit, "se.fit")
  check_flag(allow_new_levels, "allow_new_levels")
  if (!is.null(newdata) && !is.data.frame(newdata)) {
    frm_stop("`newdata` must be a data frame, or NULL to predict on the ",
             "training data, not ", arg_desc(newdata), call. = FALSE)
  }
  check_re_form(re_formula)
  type <- frm_match_arg(type)
  # the structure gate below asks whether the caller set re_formula at
  # all, which the resolution would hide by turning a partial formula
  # into NULL on a reduced design
  re_asked <- !is.null(re_formula)
  rr <- re_resolve(object, re_formula, "frm_linpred()")
  object <- rr$fit
  re_formula <- rr$re_formula
  use_re <- re_form_keeps(re_formula)

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
      frm_stop("se.fit is not supported on the response scale for an ",
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
      frm_stop("se.fit is not supported on the response scale for a ",
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
    if (re_asked) {
      structure_gate(st, "re_form",
                     structure_generic(
                       fam_, "re_formula = on the response scale"))
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
      frm_stop("type = '", type, "' selects its own dpar; drop dpar =",
               call. = FALSE)
    }
    dpar <- switch(type,
      zprob = ,
      zlink = intersect(c("zi", "hu"), names(rspec$dpars))[1],
      disp = intersect(c("sigma", "shape", "phi"),
                       names(rspec$dpars))[1]
    )
    if (is.na(dpar)) {
      frm_stop("type = '", type, "' needs a family with a ",
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
    return(predict_mean_response(object, rspec, newdata, re_formula,
                                 allow_new_levels))
  }
  dpar <- dpar %||% if ("mu" %in% names(rspec$dpars)) "mu" else
    rspec$primary_dpars[1]
  key <- linpred_key(resp, dpar)
  lp <- object$frame[["linpreds"]][[key]]
  if (is.null(lp)) {
    frm_stop("Unknown dpar: '", dpar, "' for response '", resp,
             "'. Available: ", paste(names(rspec$dpars), collapse = ", "),
             call. = FALSE)
  }

  if (!is.null(lp[["nl_body"]])) {
    if (se.fit) {
      frm_stop("se.fit is not supported for the nonlinear predictor yet; ",
               "request the nonlinear parameters (dpar = '",
               rspec$nlpars[1], "', ...) instead", call. = FALSE)
    }
    # only the parameters this body names, each on its own linear-
    # predictor scale. A parameter that carries a body of its own comes
    # back through this same branch, so a chain of nlf() formulas
    # unwinds by recursion in dependency order.
    vals <- list()
    for (np in lp[["nl_pars"]]) {
      vals[[np]] <- frm_linpred(object, newdata = newdata, dpar = np,
                                resp = resp, re_formula = re_formula,
                                allow_new_levels = allow_new_levels)
    }
    # a reference to another dpar reads its VALUE, so it comes back
    # through that parameter's link inverse - the NATURAL scale, which
    # is what the body computes with, not a reporting scale
    for (dr in lp[["nl_dpar_refs"]] %||% character(0)) {
      lpr <- object$frame[["linpreds"]][[linpred_key(resp, dr)]]
      vals[[dr]] <- lpr[["link"]]$linkinv(
        frm_linpred(object, newdata = newdata, dpar = dr,
                    type = "link", resp = resp, re_formula = re_formula,
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
    nat <- dpars_natural(object, rspec, newdata, re_formula,
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
    frm_warning("Rank-deficient fit: ", sum(nonest), " row(s) of newdata are ",
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
  frm_warning("The fitted objective marginalizes the random effects ",
              "(quadrature = TRUE), so its covariance carries no ",
              "random-effect block: se.fit is conditional on the ",
              "conditional modes and omits random-effect uncertainty. ",
              "Refit with quadrature = FALSE for the full delta method",
              call. = FALSE)
  invisible(NULL)
}

#' Delta method: `var(eta) = A V A'` over the estimated coefficients (and b
#' when random effects are included). Returns A and the positions of its
#' columns in the joint covariance, so several linear predictors can be
#' combined.
#'
#' The Z matrices span COEFFICIENT space, so their b columns are only
#' `Z` itself while `b` IS the coefficient vector. Two block types make
#' it something else: an `rr` block expands through its loadings (whose
#' `theta` parameters then contribute columns of their own), and an
#' `esicar` block is centered per connected component. Both go through
#' `d cvec/d b` from `rr_jacobians()`, and the need for one is derived
#' from the frame the way `expand_b()` derives it, with the caller's
#' `has_rr` only as a fast path: a cached flag can go stale, and pairing
#' `Z` with `b` through the identity is silent when it is wrong.
#'
#' @noRd
lp_delta_A <- function(object, lp, ed, newdata, use_re, jc, has_rr, rrj) {
  est <- object$estimates
  rn <- jc$names
  X <- ed[["X"]]
  if (is.null(rrj) && (has_rr || frame_needs_expand(object$frame))) {
    rrj <- rr_jacobians(object)
  }
  use_jac <- !is.null(rrj)
  add_b_cols <- function(A, coef_pos, Zc, b_pos, th_pos) {
    if (use_jac) {
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
    # a missing level label (conditional_effects() grids write NA for
    # "a group the fit did not see") is its own unseen level per row;
    # left NA, the row matched no key below and lost its variance
    if (anyNA(kv)) kv[is.na(kv)] <- paste0(".na.", which(is.na(kv)))
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
#' block) are part of the answer, not an afterthought. The mean is at
#' the random-effect modes and the b block of the joint covariance puts
#' their uncertainty in the standard error, as for eta.
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
      frm_stop("se.fit is not supported on the response scale for a ",
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
    frm_stop("se.fit is not available for the expected response of family '",
             fam[["family"]], "': its mean is not one number per observation",
             call. = FALSE,
             package = frm_family_package(fam))
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
  frm_check_dots(...)
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
#' brms's `fitted()`: a summary of the expected response, in the columns
#' `Estimate`, `Est.Error`, `Q2.5` and `Q97.5`. The estimate is the
#' modelled response at the estimates, conditional on the random-effect
#' modes, which is `frm_linpred(object, type = "response")` on the
#' training data for every family.
#'
#' A maximum-likelihood fit has no draws to summarize, so `Est.Error` is
#' the delta-method standard error [frm_linpred()] reports for the same
#' quantity, and the `Q` columns are the Wald interval at those
#' probabilities. The estimate is at the modes, and the standard error
#' carries the uncertainty in them: at a grouping level the fit saw,
#' the level's group effect enters with its covariance taken jointly
#' with the fixed effects from the joint precision, the frequentist
#' analogue of the posterior of that effect in brms's
#' `posterior_epred()`. The interval covers the expected response at a
#' known level with the nominal coverage averaged over the groups, not
#' for one group's realized effect. The interval is around the EXPECTED
#' response and
#' carries no observation noise; [predict.frmtmb_fit()] is the
#' predictive interval that does.
#'
#' A predictor whose standard error this package cannot produce, such
#' as a nonlinear body or a structured likelihood, reports `NA` in
#' `Est.Error` and in the `Q` columns rather than a number it does not
#' have. The estimate is unaffected.
#'
#' @param object A `frmtmb_fit`.
#' @param newdata Optional data frame to evaluate on. Defaults to the
#'   training data.
#' @param re_formula Which group-level terms enter the answer: `NULL`
#'   (default) keeps all of them, `NA` keeps none, and a one-sided
#'   formula keeps the terms it names, as in brms (see [frm_linpred()]).
#'   brms's spelling, and the only one: lme4's `re.form` is not
#'   accepted here.
#' @param scale `"response"` (default) for the modelled response, or
#'   `"linear"` for the linear predictor. brms's spelling of what
#'   [predict.frmtmb_fit()] calls `type`.
#' @param resp For multivariate fits: which response (defaults to the
#'   first).
#' @param dpar Which distributional parameter to report instead of the
#'   mean.
#' @param nlpar brms's name for a non-linear parameter, which is a
#'   distributional parameter here: a synonym for `dpar`, and giving
#'   both is an error.
#' @param ndraws,draw_ids,sort,summary,robust brms's arguments, in
#'   brms's positions so that a positional brms call asks the same
#'   question. Each needs posterior draws, and a maximum-likelihood fit
#'   has none, so each is refused by name with the reason and with the
#'   place it does work: `frmtmb.sample`'s `posterior_epred()`. The
#'   default of each is accepted and changes nothing.
#' @param probs Probabilities of the two quantile columns. brms takes
#'   the posterior quantiles there; here they are the ends of the Wald
#'   interval at those probabilities.
#' @param allow_new_levels Predict unseen grouping-factor levels at the
#'   population level instead of erroring, as in
#'   [frm_linpred()].
#' @param ... Refused. An argument this method does not have is an
#'   error naming it, because a swallowed `re_formula` returned the
#'   conditional fit and said nothing.
#' @return An `n x 4` matrix with the columns `Estimate`, `Est.Error`
#'   and one per entry of `probs`. For an ordinal or categorical family
#'   an `n x 4 x K` array, the third dimension named `P(Y = k)`, which
#'   is brms's shape. The ROW dimnames are `NULL`, as brms's are; the
#'   data's row names are on `frm_linpred()` and `model.frame()`.
#' @section Ordinal responses:
#' An ordinal response has no mean, so `fitted()` summarizes the `K`
#' category probabilities, with rows of `Estimate` summing to one,
#' which is the brms `fitted()` convention. `cs()` terms are honored.
#' The standard error of a category probability is the
#' finite-difference delta method
#' over the whole outer parameter vector, because the probability
#' depends on the thresholds and the `cs()` coefficients as well as on
#' the linear predictor. On a mixed ordinal fit the group effects join
#' the differenced vector, with their covariance taken jointly with the
#' parameters from the joint precision, so the standard error carries
#' their uncertainty as the scalar route does. The estimates are
#' unaffected. The latent linear predictor, which is where the
#' coefficients live, is `frm_linpred(object, type = "link")`.
#' @seealso [predict.frmtmb_fit()] for the predictive interval,
#'   [frm_linpred()] for the linear predictor,
#'   [residuals.frmtmb_fit()], and [frmtmb-scales] for which scale each
#'   method reports
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100))
#' dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x))
#' fit <- frm(bf(y ~ x) + poisson(), data = dd)
#' head(fitted(fit))
#' max(abs(fitted(fit)[, "Estimate"] -
#'           frm_linpred(fit, type = "response")))
#' @export
fitted.frmtmb_fit <- function(object, newdata = NULL, re_formula = NULL,
                              scale = c("response", "linear"),
                              resp = NULL, dpar = NULL, nlpar = NULL,
                              ndraws = NULL, draw_ids = NULL, sort = FALSE,
                              summary = TRUE, robust = FALSE,
                              probs = c(0.025, 0.975), ...,
                              allow_new_levels = FALSE) {
  frm_check_dots(..., .unsupported = fitted_no_draws)
  scale <- frm_match_arg(scale)
  fitted_refuse_draws_args("fitted()", ndraws, draw_ids, sort, summary,
                           robust)
  check_flag(allow_new_levels, "allow_new_levels")
  # brms reaches a non-linear parameter with `nlpar`; in this package an
  # nlf() parameter IS a dpar, so the two names address one thing
  if (!is.null(nlpar)) {
    if (!is.null(dpar)) {
      frm_stop("fitted(): give `dpar` or `nlpar`, not both. A non-linear ",
               "parameter is a distributional parameter in this package, ",
               "so the two names address the same thing", call. = FALSE)
    }
    dpar <- nlpar
  }
  # frm_linpred() defaults an unnamed multivariate response to the
  # first; fitted() refuses instead, as it always has, because the
  # caller who did not name one is asking for all of them
  if (is.null(resp)) single_response(object, "fitted()")
  check_re_form(re_formula)
  rs <- object$spec$responses[[resp %||% names(object$spec$responses)[1L]]]
  if (is.null(rs) || is.null(fam_structure(rs$family))) {
    # once, so the finite-difference route below perturbs the reduced
    # design rather than resolving the formula again at every step; a
    # structured family keeps the formula, which frm_linpred() gates
    rr <- re_resolve(object, re_formula, "fitted()")
    object <- rr$fit
    re_formula <- rr$re_formula
  }
  est <- fitted_point(object, newdata, re_formula, scale, resp, dpar,
                      allow_new_levels)
  se <- fitted_point_se(object, newdata, re_formula, scale, resp, dpar,
                        allow_new_levels, est)
  if (is.matrix(est)) {
    # a category distribution: one summary layer per category, brms's
    # n x 4 x K array, named by the response's own categories the way
    # brms's posterior_epred() names them
    out <- brms_summary_array(est, se, probs,
                              third = brms_category_labels(colnames(est),
                                                           ncol(est)))
    dn <- dimnames(out)
    dn[1L] <- list(NULL)
    dimnames(out) <- dn
    return(prob_clamp_quantiles(out))
  }
  out <- brms_summary_matrix(est, se, probs, rownames = NULL)
  out
}

#' Hold the quantile columns of a category-probability summary between
#' 0 and 1.
#'
#' The `Q` columns of a maximum-likelihood summary are a Wald interval,
#' `est + qnorm(p) * se`, which is unbounded; brms's are quantiles OF
#' probabilities and cannot leave the range. Measured on a three-
#' category ordinal fit, 15 of 450 bounds were below 0 and 5 above 1
#' (`dev/reviews/20260918-shapes.md`, m1). Clamping keeps the reported
#' bound inside the range the quantity lives in; it does NOT widen the
#' interval, so a bound that was clamped is a sign that the normal
#' approximation is poor there.
#'
#' @noRd
prob_clamp_quantiles <- function(out) {
  q <- setdiff(seq_len(dim(out)[2L]), 1:2)
  if (length(q)) {
    out[, q, ] <- pmin(pmax(out[, q, , drop = FALSE], 0), 1)
  }
  out
}

#' The fitted value itself, which is what every internal caller wants.
#'
#' One implementation, not two: every documented identity between
#' `fitted()` and `frm_linpred(type = "response")` was previously a
#' claim about two bodies that happened to agree, and the ordinal,
#' categorical and structured branches were written out twice.
#'
#' @noRd
fitted_point <- function(object, newdata = NULL, re_formula = NULL,
                         scale = "response", resp = NULL, dpar = NULL,
                         allow_new_levels = FALSE) {
  frm_linpred(object, newdata = newdata,
              type = if (scale == "response") "response" else "link",
              dpar = dpar, resp = resp, re_formula = re_formula,
              allow_new_levels = allow_new_levels)
}

#' The standard error of a fitted value, or `NULL` where none is
#' available.
#'
#' Two routes, and neither is new arithmetic. A scalar fitted value has
#' the delta-method standard error `frm_linpred(se.fit = TRUE)` already
#' computes, which carries the fixed-effect block, the random-effect
#' block and the `gp()` kriging variance. A CATEGORY DISTRIBUTION has
#' none, because a category probability depends on the thresholds and on
#' the `cs()` coefficients as well as on the linear predictor, so it
#' takes the finite-difference delta method over the whole outer
#' parameter vector.
#'
#' A predictor `se.fit` cannot answer, such as a nonlinear body or a
#' structured likelihood, reports `NA` rather than a number it does
#' not have.
#'
#' @noRd
fitted_point_se <- function(object, newdata, re_formula, scale, resp, dpar,
                            allow_new_levels, est) {
  if (is.matrix(est)) {
    f <- function(fit) {
      fitted_point(fit, newdata, re_formula, scale, resp, dpar,
                   allow_new_levels)
    }
    # a category probability moves with the group effects of the levels
    # the fit saw, so they join the differenced vector: the scalar route
    # below carries them through the joint covariance, and this route
    # carried the outer parameters alone.
    # Only the levels THESE ROWS load are differenced. Every other kept
    # level has a derivative of exactly zero here and would cost two
    # model evaluations, which on a fit with many levels is the whole
    # cost. `NULL` from re_used_b() means the design could not be
    # rebuilt, and then every kept level is differenced, as before.
    b_idx <- if (re_form_keeps(re_formula)) {
      gov <- re_governed_b(object)
      used <- re_used_b(object, newdata, resp, allow_new_levels)
      if (is.null(used)) gov else intersect(gov, used)
    }
    # In sample every level IS loaded, so the bound above cannot help
    # there; what does is that the levels of ONE block can be perturbed
    # together and attributed by row, exactly (re_b_batches()).
    bt <- if (length(b_idx)) {
      re_b_batches(object, newdata, resp, allow_new_levels, b_idx)
    }
    return(fit_fd_se(object, f, b_idx = b_idx, b_batch = bt))
  }
  p <- tryCatch(suppressWarnings(
    frm_linpred(object, newdata = newdata,
                type = if (scale == "response") "response" else "link",
                dpar = dpar, resp = resp, re_formula = re_formula,
                allow_new_levels = allow_new_levels, se.fit = TRUE)),
    error = function(e) NULL)
  if (is.null(p)) NULL else p$se.fit
}

#' The draws-only arguments of `fitted()` and `residuals()`, refused
#' when they are set to anything but their default.
#'
#' @noRd
fitted_refuse_draws_args <- function(what, ndraws, draw_ids, sort, summary,
                                     robust) {
  check_flag(summary, "summary")
  check_flag(robust, "robust")
  if (!is.null(ndraws)) {
    frm_stop(what, " cannot honor `ndraws`: ", fitted_no_draws[["ndraws"]],
             call. = FALSE)
  }
  if (!is.null(draw_ids)) {
    frm_stop(what, " cannot honor `draw_ids`: ",
             fitted_no_draws[["draw_ids"]], call. = FALSE)
  }
  if (!isFALSE(sort)) {
    frm_stop(what, " cannot honor `sort`: ", fitted_no_draws[["sort"]],
             call. = FALSE)
  }
  if (!summary) {
    frm_stop(what, " cannot honor summary = FALSE: ",
             fitted_no_draws[["summary"]], call. = FALSE)
  }
  if (robust) {
    frm_stop(what, " cannot honor robust = TRUE: ",
             fitted_no_draws[["robust"]], call. = FALSE)
  }
  invisible(NULL)
}

#' Arguments `fitted.brmsfit()` has that a point estimate cannot answer.
#'
#' Named rather than reported as unknown: `ndraws` is a real brms
#' argument, and "there is no such argument" would send the caller
#' looking for a typo that is not there.
#'
#' @noRd
fitted_no_draws <- c(
  ndraws = paste("a maximum likelihood fit carries one estimate, not a",
                 "posterior, so there is nothing to thin. Draw from the",
                 "fit with frmtmb.sample::frm_sample() and use",
                 "posterior_epred(ndraws = ) there"),
  draw_ids = paste("a maximum likelihood fit has no draws to index.",
                   "frmtmb.sample's posterior_epred() takes draw_ids"),
  sort = paste("rows come back in the order of the data, always: there",
               "is no draws dimension here to sort against.",
               "frmtmb.sample's posterior_epred() takes sort"),
  summary = paste("brms returns the posterior draws there, and a",
                  "maximum-likelihood fit has none: the Estimate column",
                  "IS the fitted value and there is nothing behind it.",
                  "frmtmb.sample's posterior_epred() takes summary"),
  robust = paste("a median over draws needs draws.",
                 "frmtmb.sample's posterior_epred() takes robust")
)

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
      frm_stop("cs() term '", ct$label, "' evaluated to ", length(v),
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
    frm_stop("type = \"response\" is not supported for an ordinal family ",
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
      frm_stop("type = \"response\" is not supported for a categorical ",
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
#' the coefficient bookkeeping is exactly `frm_linpred(se.fit = TRUE)`'s;
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
#' `discrete` refuses rather than shifts. The window this builds is the
#' continuous one, `lo < Y < hi` written as `[lo, hi]` because the ends
#' carry no mass there. A count censored at `hi` under the inclusive
#' convention (see `row_lpdf()`) can BE `hi` and still be censored, so
#' an uncensored row's support is `[lo + 1, hi - 1]` and its PIT
#' renormalizes on `F(hi - 1) - F(lo)`. Shifting the window is one line
#' and verifying it is not, so it is refused until something measures
#' it.
#'
#' @noRd
osa_cens_domain <- function(av, y, discrete = FALSE) {
  cen <- av[["cens"]]
  if (is.null(cen) || !any(cen != 0)) return(NULL)
  if (discrete) {
    frm_stop("residuals(type = \"osa\") is not supported on a cens() fit ",
             "with a discrete family. A discrete censoring bound is ",
             "INCLUSIVE (right censoring at k is Y >= k), so an uncensored ",
             "row's support is [lo + 1, hi - 1] rather than the [lo, hi] ",
             "the one-step window is built on, and no reference has ",
             "measured the shifted window. Every other residual type ",
             "works, and dharma_residuals() covers the same ground through ",
             "simulate(censored = TRUE)", call. = FALSE)
  }
  if (any(cen == 2)) {
    frm_stop("residuals(type = \"osa\") does not support interval censoring ",
             "(cens code 2): an interval-censored row observes an event, not ",
             "a value, and the uncensored rows' observation window is then ",
             "not a single interval", call. = FALSE)
  }
  i_obs <- which(cen == 0)
  if (!length(i_obs)) {
    frm_stop("residuals(type = \"osa\") needs at least one uncensored ",
             "observation", call. = FALSE)
  }
  point <- function(idx, side) {
    p <- unique(y[idx])
    if (length(p) > 1L) {
      frm_stop("residuals(type = \"osa\") on a cens() fit needs one ",
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
    frm_stop("residuals(type = \"osa\") on a cens() fit found uncensored ",
             "responses outside the censoring window [", lo, ", ", hi,
             "]; the censoring is not type-I and the one-step CDF has no ",
             "well-defined domain", call. = FALSE)
  }
  list(lo = lo, hi = hi, subset = i_obs, conditional = c(i_l, i_r))
}

#' Arguments `residuals.brmsfit()` has that this one does not.
#'
#' `newdata` and `re_formula` are not refusals of principle: a residual
#' needs an observed response, and evaluating one away from the fitted
#' rows is a feature this package does not have yet. They are named so
#' that a ported brms call is told which it is.
#'
#' @noRd
residuals_unsupported <- c(
  newdata = paste("a residual needs the observed response, and",
                  "residuals() here reads the fitted rows only.",
                  "Compute it yourself from predict(newdata = )"),
  re_formula = paste("residuals() here is conditional on the",
                     "random-effect modes, always. For the population",
                     "residual take y - predict(re_formula = NA,",
                     "type = \"response\")"),
  method = paste("brms's `method` chooses the predictive distribution;",
                 "this package has one. The one-step-ahead algorithm is",
                 "chosen with `osa_method`"),
  resp = paste("residuals() is not supported for multivariate fits at",
               "all yet, so there is no response to choose")
)

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
#' censoring are refused, and so is a DISCRETE family. A discrete
#' censoring bound is inclusive (right censoring at `k` is `Y >= k`;
#' see [frmtmb_family()]), so an uncensored count's support is
#' `[lo + 1, hi - 1]` rather than the `[lo, hi]` this window is built
#' on, and no reference has measured the shifted window. Every other
#' residual type works there.
#' `dharma_residuals()` is not a substitute on a
#' censored fit, because [simulate.frmtmb_fit()] draws the latent
#' uncensored response by default (as brms's `posterior_predict()`
#' does) and those draws are not comparable with the observed censored
#' values; `simulate(censored = TRUE)` makes them comparable, but the
#' resulting point mass at each censoring point is not a distribution
#' DHARMa's rank transform can use.
#'
#' @section Ordinal and other category responses:
#' `"response"` (`"ordinary"`) and `"pearson"` are refused for an
#' ordinal family (`cumulative()`, `sratio()`, `cratio()`, `acat()`), a
#' `categorical()` family and a `multinomial()` family. brms refuses
#' the same two types for the same families ("Predictive errors are not
#' defined for ordinal or categorical models"): the response is a
#' category or a vector of counts over categories, so `y - E[Y]` has no
#' scale to be read on.
#'
#' On an ordinal fit, `"osa"` gives a residual that uses only the
#' order. It uses `"oneStepGeneric"` over the discrete support `1..K`,
#' which makes the residuals randomized quantile residuals.
#' [dharma_residuals()] is the simulation-based alternative.
#' `"deviance"` is refused, as it is for every family without a
#' standard unit deviance.
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
#' @param type `"response"` (brms spells the same thing `"ordinary"`,
#'   and both are accepted), `"pearson"`, `"deviance"`, or `"osa"`.
#'   brms's `residuals()` takes `newdata` in this position; this one
#'   has always taken `type` there, and `newdata` is refused by name.
#' @param osa_method Method for [TMB::oneStepPredict()]; defaults to
#'   `"fullGaussian"` for gaussian models and `"oneStepGeneric"`
#'   otherwise. A truncated, censored or ordinal response always uses
#'   `"oneStepGeneric"` (a truncated gaussian is not gaussian) with the
#'   integration domain and discrete support taken from the `trunc()`
#'   bounds or the censoring window, which must then be the same for
#'   every row.
#' @param ... For `type = "osa"`: passed to [TMB::oneStepPredict()],
#'   and checked against that function's own formals. For every other
#'   type: refused, naming the argument. The `residuals.brmsfit()`
#'   arguments this one does not have (`newdata`, `re_formula`,
#'   `method`, `resp`, `ndraws`, `draw_ids`, `sort`, `summary`,
#'   `robust`, `probs`) are refused with the reason rather than
#'   reported as unknown names.
#' @param ndraws,draw_ids,sort,summary,robust brms's arguments. Each
#'   needs posterior draws and a maximum-likelihood fit has none, so
#'   each is refused by name with the reason. The default of each is
#'   accepted and changes nothing.
#' @param probs Probabilities of the two quantile columns, the ends of
#'   the Wald interval at those probabilities.
#' @return brms's summary matrix: `n` rows with `NULL` dimnames, as
#'   brms's are, and the columns `Estimate`, `Est.Error` and one per
#'   entry of `probs`, `NA` on censored rows.
#'   The observed response is fixed, so `Est.Error` is the standard
#'   error of the fitted value; a `"deviance"` or `"osa"` residual has
#'   none and reports `NA` there.
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
#' pr <- residuals(fit, type = "pearson")[, "Estimate"]
#' sum(pr^2) / df.residual(fit)
#'
#' # one-step-ahead quantile residuals are standard normal under a
#' # correctly specified model, whatever the family. This block is
#' # \donttest{} because oneStepPredict() is 5.9 s of this example's
#' # 8.9 s, which puts the examples phase over R CMD check's 5 s
#' # threshold; --run-donttest still runs it.
#' \donttest{
#' r <- residuals(fit, type = "osa")[, "Estimate"]
#' qqnorm(r); qqline(r)
#' }
#' @seealso [frmtmb-scales] for which scale each type is on.
#' @export
residuals.frmtmb_fit <- function(object, type = c("response", "ordinary",
                                                  "pearson", "deviance",
                                                  "osa"),
                                 osa_method = NULL, ...,
                                 ndraws = NULL, draw_ids = NULL,
                                 sort = FALSE, summary = TRUE,
                                 robust = FALSE, probs = c(0.025, 0.975)) {
  type <- frm_match_arg(type)
  # brms spells the raw residual "ordinary"; this package has always
  # spelled it "response", and both reach the same branch. The
  # spelling the CALLER used is what a refusal has to quote back
  written <- type
  if (identical(type, "ordinary")) type <- "response"
  fitted_refuse_draws_args("residuals()", ndraws, draw_ids, sort, summary,
                           robust)
  # the dots are checked HERE, not in residual_values(): the message and
  # the S3 contract table (`na.rm` on residuals()) are both keyed by the
  # name the caller typed, and residual_values() is not it
  frm_check_dots(..., .unsupported = residuals_unsupported,
                 .allow = if (identical(type, "osa")) {
                   names(formals(TMB::oneStepPredict))
                 })
  if (type %in% c("response", "pearson")) {
    # brms refuses its two residual types for every polytomous family
    # (brms 2.23.0 .predictive_error(), dev/correct-log/brms-resid.txt).
    # A categorical fit is left to residual_values(), which refuses
    # every type there, "osa" included, with its own reason.
    fam <- single_response(object, "residuals()")$family
    if (fam_is_polytomous(fam) &&
          !identical(fam[["type"]], "categorical")) {
      frm_stop("residuals(type = \"", written, "\") is not defined for ",
               "the '", fam[["family"]], "' family, as in brms: the ",
               "response is a category or a set of counts over ",
               "categories, so ",
               "y - E[Y] has no scale to be read on. On an ordinal fit, ",
               "type = \"osa\" gives randomized quantile residuals that ",
               "use only the order, and dharma_residuals() is the ",
               "simulation-based check", call. = FALSE)
    }
  }
  r <- residual_values(object, type = type, osa_method = osa_method, ...)
  se <- residual_point_se(object, type, r)
  if (is.matrix(r)) {
    # a matrix-valued response (multinomial counts, the mvn mixture's
    # columns, lca item codes) has one residual per CELL, so the summary
    # is the same n x 4 x K array fitted() returns there
    return(brms_summary_array(r, se, probs))
  }
  # brms's row dimnames are NULL here (measured,
  # dev/shapes-rev-brmsref.rds); carrying the data's row names made
  # rownames(residuals(fit)) differ from rownames(residuals(brmsfit))
  brms_summary_matrix(r, se, probs, rownames = NULL)
}

#' The standard error of a residual, or `NULL` where none is defined.
#'
#' The observed response is fixed, so the spread of `y - mu` is the
#' spread of `mu`, which is the quantity [fitted()] reports. A Pearson
#' residual divides by the row's residual standard deviation, so its
#' standard error divides by the same number; that divisor is read off
#' the two residual types rather than recomputed, and a row whose raw
#' residual is exactly zero carries no ratio and reports `NA`.
#'
#' A deviance or one-step-ahead residual is a nonlinear transform of the
#' whole row likelihood, and neither has a standard error here.
#'
#' @noRd
residual_point_se <- function(object, type, r) {
  if (!type %in% c("response", "pearson")) return(NULL)
  se_mu <- tryCatch(fitted_point_se(object, NULL, NULL, "response", NULL,
                                    NULL, FALSE, r),
                    error = function(e) NULL)
  if (is.null(se_mu)) return(NULL)
  if (is.matrix(r)) {
    # a cell of a matrix-valued response has no Pearson scale to divide
    # by, so only the raw residual carries one
    return(if (identical(type, "response") && is.matrix(se_mu)) {
      se_mu
    } else NULL)
  }
  if (is.matrix(se_mu)) return(NULL)
  if (identical(type, "response")) return(se_mu)
  r0 <- residual_values(object, type = "response")
  sc <- r0 / r
  sc[!is.finite(sc)] <- NA_real_
  se_mu / sc
}

#' @noRd
residual_values <- function(object, type = c("response", "pearson",
                                             "deviance", "osa"),
                            osa_method = NULL, ...) {
  type <- frm_match_arg(type)
  # The dots reach TMB::oneStepPredict() and ONLY on the osa branch.
  # They are CHECKED by residuals(), the public method, whose name is
  # what the refusal has to print.
  rspec <- single_response(object, "residuals()")
  fam <- rspec$family
  if (identical(fam[["type"]], "categorical")) {
    # nothing here is defined on a nominal scale: the categories carry
    # no order, so there is no y - E[Y] to form and no CDF to invert
    frm_stop("residuals() is not defined for a categorical family: the ",
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
    if (type == "deviance" && !is.null(st[["loglik_row"]]) &&
        is.null(st[["fitted_mean"]])) {
      # only a structure that supplies the magnitude through loglik_row
      # is told the sign is missing; a rowwise family with neither falls
      # through to the rowwise path and its own refusal
      # the structured deviance path runs inside the fitted_mean branch
      # below, so without one the family fell through to the rowwise
      # path and was told it had no unit deviance, which is false: the
      # magnitude is in loglik_row(); only the sign is missing
      frm_stop("residuals(type = \"deviance\") needs the sign of each row's ",
               "departure from its conditional mean, and the '",
               fam[["family"]], "' family declares no fitted_mean(). The ",
               "magnitude is available from loglik_row(); the sign is not. ",
               "Declare fitted_mean() in frmtmb_structure() to enable it.",
               call. = FALSE,
               package = frm_family_package(fam))
    }
    fm <- st[["fitted_mean"]]
    if (!is.null(fm)) {
      blk <- frame_block_of(object$frame, rspec$resp_name)
      r <- object$frame[["y"]][[rspec$resp_name]] - fm(object, blk)
      if (type == "pearson") {
        fv <- st[["fitted_var"]]
        if (is.null(fv)) {
          frm_stop("residuals(type = \"pearson\") needs the conditional ",
                   "variance of each row given the whole response, and the ",
                   "'", fam[["family"]],
                   "' family declares no fitted_var(). Use ",
                   "type = \"response\"", call. = FALSE,
                   package = frm_family_package(fam))
        }
        r <- r / sqrt(fv(object, blk))
      }
      if (type == "deviance") {
        # sign from the conditional mean, magnitude from the family's
        # own per-row log-density. Reached only for a family that
        # declares supports$deviance, which the gate above has already
        # checked; what is checked here is that it can produce the two
        # halves of a unit deviance.
        # The FAMILY applies the row weights: the core passes them into
        # the slot and does not reapply them, which is the convention
        # the importance correction already uses on the same slot
        # (R/importance.R). Multiplying here as well made a family that
        # follows the documented instruction report deviance residuals
        # exactly sqrt(w) too large.
        r <- sign(r) * sqrt(structure_unit_deviance(object, rspec, st, blk))
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
      frm_stop("residuals(type = \"osa\") is not available for a fit with a ",
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
    cb <- osa_cens_domain(av0, object$frame[["y"]][[rspec$resp_name]],
                          identical(fam[["type"]], "discrete"))
    ordinal <- identical(fam[["type"]], "ordinal")
    method <- osa_method %||%
      if (!is.null(tb) || !is.null(cb) || ordinal) "oneStepGeneric"
      else if (identical(fam[["family"]], "gaussian")) "fullGaussian"
      else "oneStepGeneric"
    if (!is.null(cb) && !identical(method, "oneStepGeneric")) {
      # a censored row's contribution is a probability MASS, so its
      # observation drops out of the tape; every method that
      # differentiates the observation hits a singular system there
      frm_stop("residuals(type = \"osa\") on a cens() fit needs ",
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
          frm_stop("residuals(type = \"osa\") needs trunc() bounds that are ",
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
      frm_stop("Family '", fam[["family"]], "' has no variance function; ",
               "pearson residuals are unavailable", call. = FALSE,
               package = frm_family_package(fam))
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
    frm_stop("simulate(censored = TRUE) is defined for left- and ",
             "right-censored rows; an interval-censored observation is an ",
             "interval, not a value", call. = FALSE)
  }
  point <- function(idx, side) {
    p <- unique(yobs[idx])
    if (length(p) > 1L) {
      frm_stop("simulate(censored = TRUE) needs one ", side,
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
#' model refuses `trunc()` when the frame is assembled). At `newdata` a
#' structured FAMILY is refused, because its structure indexes the rows
#' the model was fitted on: `mixture(groups = )`, a hidden Markov
#' family, a learning family, and [mixture_mvn()], which
#' `predict(newdata = )` refuses as well. A residual correlation term
#' is rebuilt on the new rows instead: rows that share a group are
#' drawn jointly, and a time the fit never saw is refused. A lag is
#' counted in the FITTED time levels, as the likelihood counts it, so
#' rows at times 2 and 4 are two levels apart even with nothing at
#' time 3 in newdata. brms counts a lag by a row's position among its
#' group's newdata rows instead, which makes the correlation of two
#' rows depend on which other rows newdata holds; frmtmb departs from
#' it on purpose.
#'
#' @section Group-level terms:
#' `re_formula` chooses which group-level terms the draws condition on,
#' read as [predict.frmtmb_fit()] reads it: `NULL` keeps every term,
#' `NA`, `~0` and `~1` keep none, and a one-sided formula keeps the terms
#' it names (a term the fit does not have is an error). A term that is
#' kept enters at its estimated effects. A term that is not kept is
#' REDRAWN from its estimated distribution in every replicate, which is
#' lme4's unconditional simulation. That is where `simulate()` and
#' `predict()` differ: a prediction is for an average group, so a
#' dropped term contributes nothing there, and a simulated response
#' needs a group, so here it gets a new one. When a formula keeps some
#' columns of a term and drops others, as `~ (1 | g)` does on a
#' `(1 + x | g)` fit, the dropped columns are drawn given the kept ones
#' at their estimates.
#'
#' A population smooth, `gp()` or `hsgp()` curve is not a group-level
#' term and is never redrawn. A factor-smooth term is, with the other
#' group-level terms.
#'
#' @section New data:
#' With `newdata` the draws are for its rows. The response column is
#' not needed. At a grouping level the fit saw, a kept term enters at
#' that level's estimate, and a redrawn term shares one draw across the
#' rows of the level, so newdata must carry the grouping column. A level
#' the fit never saw is an error unless `allow_new_levels = TRUE`, which
#' draws its effect from the term's estimated distribution, as
#' `predict()` does. Under `re_formula = NA` (or `~0`, `~1`) every term
#' is redrawn anyway, so an unseen level is one more fresh level and
#' needs nothing, except on a fit with a factor-smooth term: there an
#' unseen level takes the population curve under
#' `allow_new_levels = TRUE`, as in `predict()`, and is not redrawn.
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
#' is refused too. At `newdata` the same fitted window applies.
#'
#' @param object A `frmtmb_fit`.
#' @param nsim Number of simulated response vectors.
#' @param seed Optional RNG seed. Follows the [stats::simulate()]
#'   contract: the global RNG state is restored afterwards, and the
#'   seed used is attached as the `"seed"` attribute.
#' @param re_formula Which group-level terms the draws condition on:
#'   `NULL` (default) all of them, `NA` (or `~0`, `~1`) none, and a
#'   one-sided formula the terms it names. A term not kept is redrawn
#'   from its estimated distribution in each replicate (see Group-level
#'   terms).
#' @param censored Apply the fitted `cens()` mechanism to the draws
#'   (see Censored responses). Ignored without `cens()`.
#' @param newdata Optional data frame to simulate the responses of,
#'   instead of the fitted rows (see New data).
#' @param allow_new_levels With `newdata`: draw the effect of a grouping
#'   level the fit never saw from its term's estimated distribution,
#'   rather than refuse it. brms's spelling.
#' @param ... Refused: an argument the method does not have is an
#'   error naming it, rather than silently changing nothing.
#' @return A data frame with `nsim` columns and a `"seed"` attribute,
#'   with one row per fitted row, or per row of `newdata`.
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
#' # re_formula = NA redraws the group effects, which is the right choice
#' # for a parametric bootstrap over new groups
#' sims_m <- simulate(fit, nsim = 5, re_formula = NA, seed = 42)
#' apply(sims_m, 2, var) > apply(sims, 2, var)
#'
#' # draws for rows the fit never saw, at a known group and a new one
#' nd <- data.frame(x = c(-1, 1), g = factor(c("3", "new")))
#' simulate(fit, nsim = 3, seed = 1, newdata = nd, allow_new_levels = TRUE)
#'
#' # a posterior-predictive check by hand: does the fit reproduce the
#' # share of zeros in the data?
#' mean(dd$y == 0)
#' colMeans(simulate(fit, nsim = 20, seed = 1) == 0)
#' @export
simulate.frmtmb_fit <- function(object, nsim = 1, seed = NULL,
                                re_formula = NULL, censored = FALSE,
                                newdata = NULL, allow_new_levels = FALSE,
                                ...) {
  frm_check_dots(...)
  sim_fit_draws(object, nsim, seed, re_formula, censored, newdata,
                allow_new_levels)
}

#' The body of `simulate.frmtmb_fit()`, with one setting the public
#' method does not offer: `redraw_smooths = TRUE` makes `re_formula = NA`
#' redraw the penalized coefficients of population smooths as well,
#' which is `frm_bootstrap()`'s whole-model bootstrap and 0.62.0's
#' `simulate(re_formula = NA)` (see `sim_re_plan()`).
#'
#' @noRd
sim_fit_draws <- function(object, nsim = 1, seed = NULL, re_formula = NULL,
                          censored = FALSE, newdata = NULL,
                          allow_new_levels = FALSE,
                          redraw_smooths = FALSE) {
  # nsim reaches vapply()/replicate() as a length, where a length-2 or
  # character value reports "invalid 'length' argument" and names
  # neither simulate() nor nsim
  check_count(nsim, "nsim", min = 1L)
  check_flag(censored, "censored")
  check_flag(allow_new_levels, "allow_new_levels")
  # a value that is neither NULL, NA nor a formula used to read as
  # "condition on everything" with nothing said
  check_re_form(re_formula)
  if (!is.null(newdata) && !is.data.frame(newdata)) {
    frm_stop("simulate(): `newdata` must be a data frame, or NULL to ",
             "simulate the fitted rows, not ", arg_desc(newdata),
             call. = FALSE)
  }
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
    if (!is.null(re_formula)) {
      structure_gate(st, "re_form",
                     structure_generic(fam, "simulate(re_formula =)"),
                     context = "simulate")
    }
    if (isTRUE(censored)) {
      structure_gate(st, "cens_trunc",
                     structure_generic(fam, "simulate(censored = TRUE)"),
                     context = "simulate")
    }
  }
  if (!sim_can(fam)) {
    frm_stop("simulate(): family '", fam[["family"]], "' has no simulator yet",
             sim_note(fam), call. = FALSE,
             package = frm_family_package(fam))
  }
  plan <- sim_re_plan(object, re_formula, smooths = redraw_smooths)
  av <- object$frame[["aterm_values"]][[rspec$resp_name]]
  cwin <- NULL
  if (isTRUE(censored)) {
    if (is.null(av[["cens"]])) {
      frm_stop("simulate(censored = TRUE) needs a cens() response",
               call. = FALSE)
    }
    cwin <- cens_window(av, object$frame[["y"]][[rspec$resp_name]])
  }
  if (!is.null(newdata)) {
    sim_newdata_refuse_structure(fam)
    sim_newdata_group_cols(object, rspec, newdata, plan, allow_new_levels)
    n <- nrow(newdata)
    # the terms the family reads (trials, se, the trunc bounds) have to
    # follow the new rows; cens() does not, the draw is the latent one
    av_nd <- if (length(rspec$aterms)) {
      aterms_for_newdata(rspec, newdata)
    } else {
      list()
    }
    # every term redrawn: an unseen level is one more fresh level
    nl_ok <- allow_new_levels || isTRUE(plan$every)
    dz <- sim_newdata_design(object, rspec, newdata, nl_ok)
    # a row whose parameters are not finite at the estimates (a missing
    # covariate) draws NA, as predict() gives it, rather than NaN and a
    # warning from the family's random number generator
    ok0 <- dpv_row_finite(sim_newdata_dpars(object, rspec, newdata, nl_ok,
                                            dz))
    if (!all(ok0)) sim_newdata_warn_rows(which(!ok0))
    ac <- object$frame[["autocor"]][[rspec$resp_name]]
    ac_nd <- if (!is.null(ac)) {
      autocor_for_newdata(object, ac, rspec,
                          newdata[ok0, , drop = FALSE])
    }
    nls <- predict_new_level_spec(object, rspec, newdata, NULL, nl_ok)
    n_lost <- 0L
  } else {
    n <- stats::nobs(object)
  }
  out <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    b_use <- sim_draw_b(object, plan)
    if (is.null(newdata)) {
      dp <- with_cs_offsets(object, rspec, eval_dpars(object, b = b_use))
      dp <- dp[[rspec$resp_name]]
      ctx <- sim_context(object, rspec, dp, aterms = av, n = n,
                         extra = fit_extras(object))
      out[[s]] <- sim_draw(ctx)
    } else {
      fb <- object
      fb$estimates[["b"]] <- b_use
      dp <- sim_newdata_dpars(fb, rspec, newdata, nl_ok, dz,
                              predict_new_level_draw(fb, nls))
      ok <- ok0 & dpv_row_finite(dp)
      if (all(ok)) {
        ctx <- sim_context(fb, rspec, dp, aterms = av_nd, n = n,
                           extra = fit_extras(fb))
        ctx[["autocor"]] <- ac_nd
        out[[s]] <- sim_draw(ctx)
      } else {
        n_lost <- n_lost + sum(ok0 & !ok)
        out[[s]] <- sim_newdata_draw_rows(fb, rspec, dp, av_nd, ac_nd, ok,
                                          ok0)
      }
    }
    if (!is.null(cwin)) out[[s]] <- apply_censoring(out[[s]], cwin)
  }
  if (!is.null(newdata) && n_lost > 0L) sim_newdata_warn_draws(n_lost, nsim)
  names(out) <- paste0("sim_", seq_len(nsim))
  out <- lapply(out, function(v) {
    sim_restore_type(object, rspec, v, pad = is.null(newdata))
  })
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
#' `residuals()` keep; draws at `newdata` (`pad = FALSE`) are already one
#' per new row. `[glmmTMB test-simulate.R; lme4#737]`
#'
#' @noRd
sim_restore_type <- function(fit, rspec, v, pad = TRUE) {
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
  if (pad) napred(fit, v) else v
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
#' `frm_linpred(se.fit = TRUE)` builds a matrix `A` with one row per
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
#' @param re_formula `NULL` keeps every random effect, `NA` drops them all,
#'   a one-sided formula keeps the ones it names.
#' @param allow_new_levels Whether a grouping level the fit never saw is
#'   allowed.
#' @return A list with
#' \describe{
#'   \item{`eta`}{the linear predictor, exactly
#'     `frm_linpred(type = "link")`.}
#'   \item{`A`}{`n x p`; `d eta / d coef`.}
#'   \item{`coef_pos`}{length `p`; the rows of `V` the columns of `A`
#'     belong to, in `V`'s own order.}
#'   \item{`V`}{`p x p`; [frm_joint_cov()] subset to `coef_pos`.}
#'   \item{`coef_names`}{length `p`; the labels of those rows.}
#'   \item{`extra_var`}{length `n`; variance that is NOT coefficient
#'     uncertainty, kept separate rather than folded into `A V A'`
#'     because an exact `gp()`'s kriging variance and a new grouping
#'     level's marginal variance are not.}
#'   \item{`nonest`}{length `n`; rows that load on a direction the
#'     rank-deficient design could not identify.}
#' }
#' `var(eta)` is `rowSums((A %*% V) * A) + extra_var`, and
#' `frm_linpred(se.fit = TRUE)` is written that way.
#' @section A nonlinear body:
#' For a nonlinear predictor `A` is a JACOBIAN rather than a design, and
#' it is computed by taping the body against the coefficients it reaches
#' through. That includes a [ps()] block, whose coefficients enter the
#' body through a spline evaluated at an argument the parameters move.
#' `frm_linpred(se.fit = TRUE)` stays refused for a nonlinear predictor;
#' this is the route.
#'
#' The Jacobian is exact, and the delta method built on it is still a
#' first-order approximation, which for a warped curve is a stronger
#' assumption than it is for a linear one. `allow_new_levels = TRUE` is
#' refused there, and so is a contributing exact `gp()`, because neither
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
#' lb <- frm_lp_basis(fit, newdata = nd, re_formula = NA)
#' str(lb$A)
#' lb$coef_names
#'
#' # the covariance of the WHOLE grid, which frm_linpred() reduces to
#' # diagonal
#' Sigma <- lb$A %*% lb$V %*% t(lb$A)
#' all.equal(sqrt(diag(Sigma)),
#'           frm_linpred(fit, newdata = nd, re_formula = NA,
#'                       se.fit = TRUE)$se.fit)
#' @export
frm_lp_basis <- function(object, newdata = NULL, dpar = NULL, resp = NULL,
                         re_formula = NULL, allow_new_levels = FALSE) {
  require_frmtmb_fit(object, "frm_lp_basis()")
  require_fitted(object, "frm_lp_basis()")
  check_flag(allow_new_levels, "allow_new_levels")
  if (!is.null(newdata) && !is.data.frame(newdata)) {
    frm_stop("`newdata` must be a data frame, or NULL to use the training ",
             "data, not ", arg_desc(newdata), call. = FALSE)
  }
  check_re_form(re_formula)
  rr <- re_resolve(object, re_formula, "frm_lp_basis()")
  object <- rr$fit
  use_re <- re_form_keeps(rr$re_formula)
  resp <- resp %||% names(object$spec$responses)[1]
  rspec <- object$spec$responses[[resp]]
  if (is.null(rspec)) {
    stop_unknown_response(object, resp)
  }
  dpar <- dpar %||% if ("mu" %in% names(rspec$dpars)) "mu" else
    rspec$primary_dpars[1]
  lp <- object$frame[["linpreds"]][[linpred_key(resp, dpar)]]
  if (is.null(lp)) {
    frm_stop("frm_lp_basis(): unknown dpar '", dpar, "' for response '",
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
    frm_stop("frm_lp_basis(): allow_new_levels = TRUE is refused for a ",
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
      frm_stop("frm_lp_basis(): the nonlinear body names '", nm,
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
      frm_stop("frm_lp_basis(): the nonlinear body reaches '", nm,
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
    frm_stop("frm_lp_basis(): the nonlinear predictor '", lp[["dpar"]],
             "' reaches no estimated coefficient, so it has no design",
             call. = FALSE)
  }
  chat <- vapply(pos, function(k) est[[pm$comp[k]]][pm$idx[k]], 0)

  dl <- lp_basis_nl_data(object, newdata)
  # the same statement predict(newdata = ) makes: this door evaluates
  # the curve somewhere the fit never saw, and past a ps() knot span the
  # basis stops being a partition of unity. Collected rather than warned
  # from inside the closure because this function enters that closure
  # TWICE per term: once for the off-tape pass below and once for the
  # tape build. Warning in place would double-fire on every call. (RTMB
  # enters once; the second entry is ours.)
  #
  # NULL newdata is the in-sample case, where the fit-end coverage
  # report has already spoken, and a body with no ps() term has nothing
  # to say, so neither pays for the collector or for the extra pass.
  has_ps <- FALSE
  for (lpk in c(list(lp), lapply(parts, `[[`, "lp"))) {
    if (length(lpk[["ps_terms"]] %||% list())) has_ps <- TRUE
  }
  span <- if (is.null(newdata) || !has_ps) {
    FALSE
  } else {
    new.env(parent = emptyenv())
  }
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
                    dl[[nm]], ps_env(lpk, pars, bvec, check = span)),
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
           dl[[lp[["dpar"]]]], ps_env(lp, pars, bvec, check = span)),
         ad_overload_env(lp[["nl_env"]], lp[["nl_body"]]))
  }
  # One evaluation off the tape before the tape is built. A ps() term
  # whose argument names a nonlinear parameter arrives at the span check
  # as an advector during taping, where a comparison raises rather than
  # answers, so the taped pass alone would leave exactly the models
  # predict() warns about unchecked. Errors and incidental warnings are
  # swallowed: a diagnostic that fails must not take the basis with it.
  if (is.environment(span)) {
    tryCatch(suppressWarnings(eta_fun(chat)), error = function(e) NULL)
  }
  tp <- RTMB::MakeTape(eta_fun, chat)
  eta <- as.numeric(tp(chat))
  A <- tp$jacobian(chat)
  ps_span_flush(span)
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
