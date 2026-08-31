# frmtmb_spec + data -> frmtmb_frame: numeric design matrices, random-effect
# blocks, smooth bases, and the flat parameter template with index
# bookkeeping.
#
# Parameter layout (glmmTMB convention):
#   beta   - fixed effects of every primary (location) linear predictor of
#            every response (REML integrates these), including smooth
#            fixed (null-space) columns
#   betad  - fixed effects of all other dpar linear predictors
#   b      - all random-effect modes; each block level-major; |ID|-merged
#            blocks interleave their component coefficients within a level
#   theta  - all covariance / smoothing parameters
#   thetar - residual-correlation parameters (rescor only)
#
# Every linear predictor's Z spans the FULL b vector (sparse), so Z's
# column indices are b indices; |ID| merging is only a matter of column
# placement. Constant dpars keep their intercept column in betad but are
# excluded from estimation through the MakeADFun `map`.

# Sum of offset() terms of one linear predictor, evaluated against the
# combined model frame (whose columns are named by deparsed expressions).
extract_offset <- function(tt, mf, env) {
  oi <- attr(tt, "offset")
  if (is.null(oi) || !length(oi)) return(NULL)
  vars <- attr(tt, "variables")
  off <- 0
  for (i in oi) {
    expr <- vars[[i + 1L]]
    v <- mf[[deparse1(expr)]]
    if (is.null(v)) v <- eval(expr, mf, env)
    off <- off + as.numeric(v)
  }
  off
}

# The variables one dpar needs in the combined model frame: parametric
# terms, bar variables, and raw smooth variables (never the s() calls
# themselves, which model.frame cannot evaluate).
dpar_frame_rhs <- function(dp) {
  if (!is.null(dp$nl_body)) {
    parts <- lapply(dp$datavars, as.name)
    if (!length(parts)) return(1)
    out <- NULL
    for (p in parts) out <- if (is.null(out)) p else call("+", out, p)
    return(out)
  }
  parts <- list(reformulas::RHSForm(dp$fixed))
  for (rt in dp$re %||% list()) {
    parts <- c(parts, list(rt$bar[[2]], rt$bar[[3]]))
  }
  for (sspec in dp$smooth %||% list()) {
    for (tm in sspec$term) parts <- c(parts, list(as.name(tm)))
    if (!is.null(sspec$by) && sspec$by != "NA") {
      parts <- c(parts, list(as.name(sspec$by)))
    }
  }
  out <- NULL
  for (p in parts) out <- if (is.null(out)) p else call("+", out, p)
  out
}

extract_y <- function(resp, mf) {
  y <- mf[[deparse1(resp$resp_expr)]]
  if (is.null(y)) {
    y <- eval(resp$resp_expr, mf, resp$formula_env)
  }
  if (is.factor(y) && identical(resp$family$type, "ordinal")) {
    y <- as.numeric(y)   # category codes 1..K in level order
  } else if (is.factor(y)) {
    if (!identical(resp$family$family, "binomial") || nlevels(y) != 2L) {
      stop("Factor responses are only supported for binomial families ",
           "with 2 levels (ordinal families accept ordered factors)",
           call. = FALSE)
    }
    y <- as.numeric(y) - 1
  }
  if (is.matrix(y) && ncol(y) == 1L) {
    y <- drop(y)   # scale() and friends return n x 1 matrices (glmmTMB#937)
  }
  if (is.matrix(y)) {
    storage.mode(y) <- "double"
  } else {
    y <- as.numeric(as.vector(y))
  }
  if (any(!is.finite(y) & !is.na(y))) {
    stop("Non-finite (Inf/NaN) values in the response are not allowed",
         call. = FALSE)
  }
  y
}

# Data-dependent bases (poly, ns, scale) must be frozen at fit time: the
# combined model frame's terms carry predvars for every variable; patch
# them onto a sub-formula's terms by deparsed-variable match before any
# newdata evaluation (glmmTMB#402 and the largest bug class in
# lme4/glmmTMB prediction history).
patch_predvars <- function(tt, map) {
  if (is.null(map) || !length(map)) return(tt)
  vars <- attr(tt, "variables")
  pv <- vars
  changed <- FALSE
  for (i in seq_along(vars)[-1]) {
    key <- deparse1(vars[[i]])
    if (!is.null(map[[key]])) {
      pv[[i]] <- map[[key]]
      changed <- TRUE
    }
  }
  if (changed) attr(tt, "predvars") <- pv
  tt
}

assemble_frame <- function(spec, data, na.action = stats::na.omit) {
  # One combined model frame holds every response, every variable of every
  # dpar of every response, and the aterm variables, so na.omit keeps rows
  # aligned. Responses go on the RHS of the frame formula (a multivariate
  # LHS is not a valid model.frame response) and are extracted by name.
  rhs_comb <- NULL
  add_part <- function(part) {
    rhs_comb <<- if (is.null(rhs_comb)) part else call("+", rhs_comb, part)
  }
  for (resp in spec$responses) {
    add_part(resp$resp_expr)
    for (dp in resp$dpars) add_part(dpar_frame_rhs(dp))
    for (a in resp$aterms) {
      # literal constants (e.g. trials(10)) are not frame variables
      if (!is.numeric(a) && !is.logical(a)) add_part(a)
    }
  }
  env <- spec$responses[[1]]$formula_env
  fr_formula <- stats::as.formula(call("~", rhs_comb), env = env)
  mf <- stats::model.frame(fr_formula, data = data,
                           drop.unused.levels = TRUE,
                           na.action = na.action)
  n <- nrow(mf)
  if (n == 0L) stop("No complete observations after removing NAs",
                    call. = FALSE)
  if (anyNA(mf)) {
    stop("NA values remain in the model variables after applying ",
         "na.action; use na.omit (default) or na.exclude", call. = FALSE)
  }
  # freeze data-dependent bases: map deparsed variable -> predvar call
  predvar_map <- local({
    tt_all <- attr(mf, "terms")
    vars <- as.list(attr(tt_all, "variables"))[-1]
    pv <- as.list(attr(tt_all, "predvars") %||%
                    attr(tt_all, "variables"))[-1]
    stats::setNames(pv, vapply(vars, deparse1, ""))
  })

  y <- list()
  aterm_values <- list()
  extras <- list()
  for (resp in spec$responses) {
    y[[resp$resp_name]] <- extract_y(resp, mf)
    av <- lapply(resp$aterms, function(a) {
      v <- mf[[deparse1(a)]]
      if (is.null(v)) v <- eval(a, mf, resp$formula_env)
      as.numeric(v)
    })
    if (!is.null(av$weights)) {
      if (any(av$weights < 0)) {
        stop("weights() must be non-negative", call. = FALSE)
      }
      if (spec$rescor) {
        stop("weights() cannot be combined with rescor = TRUE",
             call. = FALSE)
      }
    }
    if (!is.null(av$cens) || !is.null(av$trunc_lb) ||
        !is.null(av$trunc_ub)) {
      if (is.null(resp$family$lcdf)) {
        stop("cens()/trunc() need a family with a CDF (currently: ",
             "gaussian, lognormal)", call. = FALSE)
      }
      if (!is.null(av$cens) && !all(av$cens %in% c(-1, 0, 1))) {
        stop("cens() codes must be -1 (left), 0 (observed), or 1 (right)",
             call. = FALSE)
      }
    }
    if (!is.null(resp$family$valid_y)) {
      resp$family$valid_y(y[[resp$resp_name]], av)
    }
    if (!is.null(resp$family$extra_pars)) {
      if (length(spec$responses) > 1) {
        stop("Families with extra parameters ('",
             resp$family$family, "') are not supported in multivariate ",
             "fits yet", call. = FALSE)
      }
      extras <- resp$family$extra_pars(y[[resp$resp_name]], av)
    }
    aterm_values[[resp$resp_name]] <- av
  }

  ## Phase 1: per-linpred design matrices and random-effect components.
  linpreds <- list()
  components <- list()
  beta_names <- character(0)
  betad_names <- character(0)
  betad_fixed_idx <- integer(0)

  for (resp in spec$responses) {
    for (dp in resp$dpars) {
      lp_key <- linpred_key(resp$resp_name, dp$name)
      is_primary <- dp$name %in% resp$primary_dpars
      par_name <- if (is_primary) "beta" else "betad"
      dp_prefix <- if (identical(dp$name, "mu")) "" else
        paste0(dp$name, ": ")
      if (length(spec$responses) > 1) {
        dp_prefix <- paste0(resp$resp_name, " ", dp_prefix)
      }

      if (!is.null(dp$nl_body)) {
        # nonlinear mu: no design of its own; evaluated from the nlpar
        # values and raw data columns inside the objective
        data_list <- lapply(stats::setNames(dp$datavars, dp$datavars),
                            function(v) {
                              val <- mf[[v]]
                              if (is.null(val)) {
                                stop("Variable '", v, "' from the ",
                                     "nonlinear formula not found",
                                     call. = FALSE)
                              }
                              val
                            })
        linpreds[[lp_key]] <- list(
          resp = resp$resp_name, dpar = dp$name, X = NULL,
          n_param_cols = 0L, Z = NULL, par = "beta", idx = integer(0),
          offset = NULL, link = dp$link, terms = NULL, xlevels = NULL,
          contrasts = NULL, smooths = list(), comp_ids = integer(0),
          constant = NULL, nl_body = dp$nl_body, data_list = data_list,
          nl_env = dp$nl_env
        )
        next
      }

      tt <- stats::terms(dp$fixed)
      X <- stats::model.matrix(tt, mf)
      contr <- attr(X, "contrasts")   # subsetting X below drops the attr
      if (isTRUE(resp$family$drop_intercept) && is_primary) {
        # thresholds replace the intercept; a threshold-only model
        # (y ~ 1) leaves X with zero columns, which is fine
        X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
      }
      # rank-deficient designs: drop aliased columns like lm() (lme4#144)
      if (ncol(X) > 1L) {
        qrX <- qr(X)
        if (qrX$rank < ncol(X)) {
          dropped <- colnames(X)[qrX$pivot[(qrX$rank + 1L):ncol(X)]]
          message("Fixed-effect design of '", lp_key,
                  "' is rank deficient; dropping column(s): ",
                  paste(dropped, collapse = ", "))
          X <- X[, setdiff(colnames(X), dropped), drop = FALSE]
        }
      }
      param_colnames <- colnames(X)
      n_param_cols <- ncol(X)
      off <- extract_offset(tt, mf, resp$formula_env)
      xlev <- stats::.getXlevels(stats::terms(mf), mf)
      comp_ids <- integer(0)

      if (length(dp$re)) {
        bars <- lapply(dp$re, `[[`, "bar")
        rt <- reformulas::mkReTrms(bars, fr = mf, reorder.terms = FALSE)
        fassign <- attr(rt$flist, "assign")
        for (k in seq_along(bars)) {
          d_k <- length(rt$cnms[[k]])
          len_k <- rt$Gp[k + 1L] - rt$Gp[k]
          cs_name <- dp$re[[k]]$covstruct
          if (cs_name %in% c("ar1", "cs") && d_k < 2L) {
            stop(cs_name, "() needs at least 2 terms per level",
                 call. = FALSE)
          }
          if (cs_name == "ar1" && "(Intercept)" %in% rt$cnms[[k]]) {
            stop("ar1() requires a factor without intercept on the left ",
                 "of the bar, e.g. ar1(times + 0 | g)", call. = FALSE)
          }
          fac <- rt$flist[[fassign[k]]]
          Zk <- Matrix::t(rt$Zt[rt$Gp[k] + seq_len(len_k), , drop = FALSE])
          components[[length(components) + 1L]] <- list(
            lp_key = lp_key, dpar = dp$name, resp = resp$resp_name,
            covstruct = cs_name, id = dp$re[[k]]$id,
            dim = d_k, n_levels = len_k %/% d_k,
            levels = levels(fac), cnms = rt$cnms[[k]],
            bar = bars[[k]], Zlocal = Zk,
            group_name = names(rt$flist)[fassign[k]],
            label = paste0(dp_prefix, deparse1(bars[[k]]))
          )
          comp_ids <- c(comp_ids, length(components))
        }
      }

      # Smooths: fixed (null-space) part into X, wiggly part as an
      # iid-Gaussian component whose variance is the inverse smoothing
      # parameter.
      sm_info <- list()
      for (sspec in dp$smooth %||% list()) {
        scl <- mgcv::smoothCon(sspec, data = mf, absorb.cons = TRUE)
        for (sm in scl) {
          re2 <- mgcv::smooth2random(sm, names(mf), type = 2)
          if (isTRUE(re2$fixed)) {
            stop("Fixed (fx = TRUE) smooths are not supported: ", sm$label,
                 call. = FALSE)
          }
          nr <- integer(0)
          sm_comp_ids <- integer(0)
          for (r in seq_along(re2$rand)) {
            Xr <- re2$rand[[r]]
            k_r <- ncol(Xr)
            components[[length(components) + 1L]] <- list(
              lp_key = lp_key, dpar = dp$name, resp = resp$resp_name,
              covstruct = "smooth", id = NULL,
              dim = k_r, n_levels = 1L,
              levels = NULL, cnms = paste0(sm$label, ".", seq_len(k_r)),
              bar = NULL, Zlocal = methods::as(Xr, "CsparseMatrix"),
              group_name = sm$label,
              label = paste0(dp_prefix, sm$label)
            )
            sm_comp_ids <- c(sm_comp_ids, length(components))
            nr <- c(nr, k_r)
          }
          comp_ids <- c(comp_ids, sm_comp_ids)
          nf <- if (is.null(re2$Xf)) 0L else ncol(re2$Xf)
          xf_idx <- integer(0)
          if (nf > 0L) {
            Xf <- re2$Xf
            colnames(Xf) <- paste0(sm$label, ".fx", seq_len(nf))
            xf_idx <- ncol(X) + seq_len(nf)
            X <- cbind(X, Xf)
          }
          sm_info[[length(sm_info) + 1L]] <- list(
            sm = sm, U = re2$trans.U, D = re2$trans.D,
            nr = nr, nf = nf, xf_idx = xf_idx,
            comp_ids = sm_comp_ids, block_ids = NULL,
            label = sm$label
          )
        }
      }

      cn <- colnames(X)
      if (!identical(dp$name, "mu")) cn <- paste(dp$name, cn, sep = "_")
      if (length(spec$responses) > 1) {
        cn <- paste(resp$resp_name, cn, sep = "_")
      }
      if (par_name == "beta") {
        idx <- length(beta_names) + seq_len(ncol(X))
        beta_names <- c(beta_names, cn)
      } else {
        idx <- length(betad_names) + seq_len(ncol(X))
        betad_names <- c(betad_names, cn)
        if (!is.null(dp$constant)) {
          betad_fixed_idx <- c(betad_fixed_idx, idx)
        }
      }

      linpreds[[lp_key]] <- list(
        resp = resp$resp_name,
        dpar = dp$name,
        X = X,
        n_param_cols = n_param_cols,
        param_colnames = param_colnames,
        Z = NULL,               # filled in phase 3
        par = par_name,
        idx = idx,
        offset = if (!is.null(off)) as.numeric(off),
        link = dp$link,
        terms = tt,
        xlevels = xlev,
        contrasts = contr,
        smooths = sm_info,
        comp_ids = comp_ids,
        constant = dp$constant
      )
    }
  }

  ## Phase 2: components -> blocks. Components sharing an |ID| key merge
  ## into one us block; everything else gets its own block.
  comp_group <- integer(length(components))
  id_keys <- vapply(components, function(cp) cp$id %||% "", "")
  group_defs <- list()
  for (ci in seq_along(components)) {
    key <- id_keys[ci]
    if (key != "" && key %in% names(group_defs)) {
      group_defs[[key]] <- c(group_defs[[key]], ci)
    } else if (key != "") {
      group_defs[[key]] <- ci
    } else {
      group_defs[[paste0(".solo", ci)]] <- ci
    }
  }

  re_blocks <- list()
  comp_block <- integer(length(components))   # component -> block index
  comp_offset <- integer(length(components))  # coef offset within level
  n_b <- 0L
  n_theta <- 0L

  for (gd in group_defs) {
    cps <- components[gd]
    if (length(cps) > 1L) {
      lv <- cps[[1]]$levels
      for (cp in cps) {
        if (!identical(cp$levels, lv)) {
          stop("|ID|-linked terms must share identical grouping-factor ",
               "levels (", cps[[1]]$label, " vs ", cp$label, ")",
               call. = FALSE)
        }
      }
      D <- sum(vapply(cps, `[[`, 0L, "dim"))
      cs_name <- "us"
      cnms <- unlist(lapply(cps, function(cp) {
        paste0(cp$lp_key, ":", cp$cnms)
      }))
      label <- paste0(paste(vapply(cps, `[[`, "", "label"),
                            collapse = " + "), " [ID]")
      n_levels <- cps[[1]]$n_levels
    } else {
      cp <- cps[[1]]
      D <- cp$dim
      cs_name <- cp$covstruct
      cnms <- cp$cnms
      label <- cp$label
      n_levels <- cp$n_levels
    }
    npar_k <- covstruct_registry[[cs_name]]$npar(D)
    blk_i <- length(re_blocks) + 1L
    ofs <- 0L
    for (k in seq_along(gd)) {
      comp_block[gd[k]] <- blk_i
      comp_offset[gd[k]] <- ofs
      ofs <- ofs + cps[[k]]$dim
    }
    re_blocks[[blk_i]] <- list(
      covstruct = cs_name,
      dim = D,
      n_levels = n_levels,
      b_idx = n_b + seq_len(D * n_levels),
      theta_idx = n_theta + seq_len(npar_k),
      levels = cps[[1]]$levels,
      cnms = cnms,
      group_name = cps[[1]]$group_name,
      term_label = label,
      dpar = cps[[1]]$dpar,
      components = lapply(seq_along(gd), function(k) {
        list(lp_key = cps[[k]]$lp_key, offset = comp_offset[gd[k]],
             dim = cps[[k]]$dim, bar = cps[[k]]$bar,
             cnms = cps[[k]]$cnms, label = cps[[k]]$label)
      })
    )
    n_b <- n_b + D * n_levels
    n_theta <- n_theta + npar_k
  }

  ## Phase 3: per-linpred Z over the full b vector; resolve smooth block
  ## ids; assemble the parameter template.
  for (key in names(linpreds)) {
    lp <- linpreds[[key]]
    if (length(lp$comp_ids)) {
      ii <- integer(0); jj <- integer(0); xx <- numeric(0)
      for (ci in lp$comp_ids) {
        cp <- components[[ci]]
        bk <- re_blocks[[comp_block[ci]]]
        Tk <- methods::as(cp$Zlocal, "TsparseMatrix")
        loc_col <- Tk@j + 1L
        lev <- (loc_col - 1L) %/% cp$dim + 1L
        cf <- (loc_col - 1L) %% cp$dim + 1L
        glob <- bk$b_idx[(lev - 1L) * bk$dim + comp_offset[ci] + cf]
        ii <- c(ii, Tk@i + 1L)
        jj <- c(jj, glob)
        xx <- c(xx, Tk@x)
      }
      lp$Z <- Matrix::sparseMatrix(i = ii, j = jj, x = xx,
                                   dims = c(n, n_b))
    }
    if (length(lp$smooths)) {
      lp$smooths <- lapply(lp$smooths, function(si) {
        si$block_ids <- comp_block[si$comp_ids]
        si
      })
    }
    lp$comp_ids <- NULL
    linpreds[[key]] <- lp
  }

  par_template <- list(beta = stats::setNames(numeric(length(beta_names)),
                                              beta_names))
  if (length(betad_names)) {
    par_template$betad <- stats::setNames(numeric(length(betad_names)),
                                          betad_names)
  }
  if (n_b) par_template$b <- numeric(n_b)
  if (n_theta) par_template$theta <- numeric(n_theta)
  if (spec$rescor) {
    K <- length(spec$responses)
    par_template$thetar <- numeric(K * (K - 1L) / 2L)
  }
  for (nm in names(extras)) {
    if (nm %in% names(par_template)) {
      stop("Extra-parameter name collides with the template: ", nm,
           call. = FALSE)
    }
    par_template[[nm]] <- extras[[nm]]
  }

  # Constant dpars: fix their betad entries at link(constant) via map.
  map <- list()
  if (length(betad_fixed_idx)) {
    mp <- seq_along(par_template$betad)
    mp[betad_fixed_idx] <- NA
    map$betad <- factor(mp)
    for (lp in linpreds) {
      if (!is.null(lp$constant)) {
        par_template$betad[lp$idx] <- lp$link$linkfun(lp$constant)
      }
    }
  }

  structure(
    list(spec = spec, n_obs = n, y = y, aterm_values = aterm_values,
         linpreds = linpreds, re_blocks = re_blocks,
         par_template = par_template, map = map,
         betad_fixed_idx = betad_fixed_idx,
         extra_names = names(extras),
         predvar_map = predvar_map,
         data_frame = mf,
         na_action = attr(mf, "na.action")),
    class = "frmtmb_frame"
  )
}

#' @export
print.frmtmb_frame <- function(x, ...) {
  cat("<frmtmb frame> ", x$n_obs, " observations, ",
      length(x$spec$responses), " response(s)\n", sep = "")
  for (lp in x$linpreds) {
    cat("  ", lp$resp, ".", lp$dpar, ": X[", nrow(lp$X), " x ",
        ncol(lp$X), "]", sep = "")
    if (!is.null(lp$Z)) cat(", Z[", nrow(lp$Z), " x ", ncol(lp$Z), "]",
                            sep = "")
    cat(" -> ", lp$par, "[", min(lp$idx), ":", max(lp$idx), "]", sep = "")
    if (!is.null(lp$constant)) cat("  (fixed at ", lp$constant, ")", sep = "")
    cat("\n")
  }
  for (bk in x$re_blocks) {
    cat("  RE block: ", bk$term_label, " [", bk$covstruct, "] dim=", bk$dim,
        " levels=", bk$n_levels, "\n", sep = "")
  }
  cat("  parameters:",
      paste0(names(x$par_template), "(",
             lengths(x$par_template), ")", collapse = ", "), "\n")
  invisible(x)
}
