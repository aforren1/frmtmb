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
  for (mexpr in dp$mo %||% list()) {
    for (v in all.vars(mexpr)) parts <- c(parts, list(as.name(v)))
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
    for (nm_at in names(resp$aterms)) {
      # literal constants (e.g. trials(10)) are not frame variables, and
      # interval bounds (cens_y2) may be NA on non-interval rows, so they
      # stay out of the na.omit frame (brms#1070)
      a <- resp$aterms[[nm_at]]
      if (nm_at != "cens_y2" && !is.numeric(a) && !is.logical(a)) {
        add_part(a)
      }
    }
  }
  env <- spec$responses[[1]]$formula_env
  fr_formula <- stats::as.formula(call("~", rhs_comb), env = env)
  # x | mi() responses may carry NAs (they become latent parameters);
  # rows are dropped only for NAs in every OTHER variable
  mi_cols <- vapply(
    Filter(function(r) isTRUE(r$aterms$mi), spec$responses),
    function(r) deparse1(r$resp_expr), ""
  )
  if (length(mi_cols)) {
    mf <- stats::model.frame(fr_formula, data = data,
                             drop.unused.levels = TRUE,
                             na.action = stats::na.pass)
    bad <- Reduce(`|`, lapply(setdiff(names(mf), mi_cols), function(cn) {
      v <- mf[[cn]]
      if (is.matrix(v)) rowSums(is.na(v)) > 0 else is.na(v)
    }), rep(FALSE, nrow(mf)))
    if (any(bad)) {
      dropped <- which(bad)
      mf <- mf[!bad, , drop = FALSE]
      attr(mf, "na.action") <- structure(
        dropped,
        class = if (identical(na.action, stats::na.exclude)) "exclude"
                else "omit"
      )
    }
  } else {
    mf <- stats::model.frame(fr_formula, data = data,
                             drop.unused.levels = TRUE,
                             na.action = na.action)
  }
  n <- nrow(mf)
  if (n == 0L) stop("No complete observations after removing NAs",
                    call. = FALSE)
  if (anyNA(mf[setdiff(names(mf), mi_cols)])) {
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
  mi_map <- list()   # per mi() response: missing rows + miss indices
  n_miss <- 0L
  miss_init <- numeric(0)
  for (resp in spec$responses) {
    y[[resp$resp_name]] <- extract_y(resp, mf)
    av <- lapply(resp$aterms[setdiff(names(resp$aterms),
                                     c("cens_y2", "se_sigma", "mi"))],
                 function(a) {
      v <- mf[[deparse1(a)]]
      if (is.null(v)) v <- eval(a, mf, resp$formula_env)
      as.numeric(v)
    })
    if (!is.null(resp$aterms$se_sigma)) {
      av$se_sigma <- resp$aterms$se_sigma   # logical flag, not data
    }
    if (isTRUE(resp$aterms$mi)) {
      if (!resp$family$family %in% c("gaussian", "student")) {
        stop("mi() responses need a gaussian or student model",
             call. = FALSE)
      }
      if (any(c("cens", "trunc_lb", "trunc_ub", "se") %in%
                names(resp$aterms))) {
        stop("mi() cannot be combined with cens(), trunc(), or se() ",
             "on the same response", call. = FALSE)
      }
      if (spec$rescor) {
        stop("mi() cannot be combined with rescor = TRUE", call. = FALSE)
      }
      yv <- y[[resp$resp_name]]
      if (is.matrix(yv)) {
        stop("mi() responses must be numeric vectors", call. = FALSE)
      }
      if (!is.null(av$mi_sd)) {
        # measurement error (brms me()): every true value is latent;
        # observed values get a N(latent, sd) term in the objective
        if (any(av$mi_sd <= 0)) {
          stop("mi(sd): measurement SDs must be positive", call. = FALSE)
        }
        obs <- which(!is.na(yv))
        if (!length(obs)) {
          stop("mi(sd): the response has no observed values",
               call. = FALSE)
        }
        rows <- seq_along(yv)
        mi_map[[resp$resp_name]] <- list(
          rows = rows, idx = n_miss + rows,
          se = av$mi_sd, obs = obs
        )
        n_miss <- n_miss + length(rows)
        miss_init <- c(miss_init,
                       ifelse(is.na(yv), mean(yv[obs]), yv))
        yv[is.na(yv)] <- 0
        y[[resp$resp_name]] <- yv
      } else {
        rows <- which(is.na(yv))
        if (length(rows)) {
          mi_map[[resp$resp_name]] <- list(rows = rows,
                                           idx = n_miss + seq_along(rows))
          n_miss <- n_miss + length(rows)
          miss_init <- c(miss_init,
                         rep(mean(yv[-rows]), length(rows)))
          yv[rows] <- 0
          y[[resp$resp_name]] <- yv
        }
      }
    }
    if (!is.null(resp$aterms$cens_y2)) {
      v <- as.numeric(eval(resp$aterms$cens_y2, data, resp$formula_env))
      if (!is.null(attr(mf, "na.action"))) {
        v <- v[-attr(mf, "na.action")]
      }
      av$cens_y2 <- v
    }
    for (vn in grep("^vint", names(av), value = TRUE)) {
      if (any(av[[vn]] != round(av[[vn]]))) {
        stop(vn, " values must be integers (use vreal() for reals)",
             call. = FALSE)
      }
    }
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
             "gaussian, lognormal, poisson)", call. = FALSE)
      }
      if (!is.null(av$cens) && identical(resp$family$type, "discrete")) {
        stop("cens() is not supported for discrete families yet ",
             "(truncation is)", call. = FALSE)
      }
      if (!is.null(av$cens) && !all(av$cens %in% c(-1, 0, 1, 2))) {
        stop("cens() codes must be -1 (left), 0 (observed), 1 (right), ",
             "or 2 (interval)", call. = FALSE)
      }
      if (!is.null(av$cens) && any(av$cens == 2)) {
        i2 <- av$cens == 2
        if (is.null(av$cens_y2)) {
          stop("Interval censoring (code 2) needs upper bounds: ",
               "cens(c, y2)", call. = FALSE)
        }
        if (anyNA(av$cens_y2[i2])) {
          stop("cens() upper bounds must not be NA on interval-censored ",
               "rows", call. = FALSE)
        }
        yv <- y[[resp$resp_name]]
        if (any(av$cens_y2[i2] <= yv[i2])) {
          stop("Interval upper bounds (y2) must exceed the lower bounds ",
               "(the response)", call. = FALSE)
        }
        # NA bounds on non-interval rows are legal and unused; make them
        # harmless for the taped CDF evaluation
        av$cens_y2[!i2 | is.na(av$cens_y2)] <-
          yv[!i2 | is.na(av$cens_y2)]
      }
    }
    if (!is.null(av$se)) {
      if (!resp$family$family %in% c("gaussian", "student")) {
        stop("se() is supported for gaussian and student families only",
             call. = FALSE)
      }
      if (any(av$se <= 0)) {
        stop("se() values must be positive", call. = FALSE)
      }
    }
    # glm/glmer compatibility: a proportion response with trials()
    # becomes integer counts before validation
    if (resp$family$family %in% c("binomial", "beta_binomial") &&
        !is.null(av$trials)) {
      yv <- y[[resp$resp_name]]
      if (is.numeric(yv) && !is.matrix(yv) && any(yv != round(yv))) {
        yc <- yv * av$trials
        if (max(abs(yc - round(yc))) > 1e-6) {
          stop(resp$family$family, ": a proportion response times ",
               "trials() must give integer counts", call. = FALSE)
        }
        y[[resp$resp_name]] <- round(yc)
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
          dist_cs <- c("ou", "exp", "gau", "mat")
          if (cs_name %in% c("ar1", "hetar1", "cs", "homcs", "toep",
                             "homtoep", "rr", dist_cs) && d_k < 2L) {
            stop(cs_name, "() needs at least 2 terms per level",
                 call. = FALSE)
          }
          if (cs_name %in% c("ar1", "hetar1", dist_cs) &&
              "(Intercept)" %in% rt$cnms[[k]]) {
            stop(cs_name, "() requires a factor without intercept on ",
                 "the left of the bar, e.g. ", cs_name,
                 "(times + 0 | g)", call. = FALSE)
          }
          fac <- rt$flist[[fassign[k]]]
          aux_A <- NULL
          aux_D <- NULL
          aux_kron <- NULL
          if (cs_name %in% dist_cs) {
            v <- all.vars(bars[[k]][[2]])
            if (length(v) != 1L || !is.factor(mf[[v]])) {
              stop(cs_name, "() needs a single factor built with ",
                   "num_factor(): ", cs_name, "(pos + 0 | g)",
                   call. = FALSE)
            }
            coords <- parse_num_levels(levels(mf[[v]]))
            if (is.matrix(coords)) {
              aux_D <- as.matrix(stats::dist(coords))
            } else {
              aux_D <- abs(outer(coords, coords, "-"))
            }
          }
          aux_Q <- NULL
          if (cs_name == "gr_prec") {
            if (d_k != 1L) {
              stop("gr(prec=) supports intercept-only terms:",
                   " (1 | gr(g, prec = Q))", call. = FALSE)
            }
            Q <- eval(dp$re[[k]]$cov_expr, data, resp$formula_env)
            if (is.null(rownames(Q)) ||
                !all(levels(fac) %in% rownames(Q))) {
              stop("gr(prec=): prec needs dimnames covering all ",
                   "grouping levels", call. = FALSE)
            }
            lv <- levels(fac)
            aux_Q <- methods::as(Matrix::Matrix(Q[lv, lv], sparse = TRUE),
                                 "generalMatrix")
          }
          if (cs_name == "gr_cov") {
            A <- eval(dp$re[[k]]$cov_expr, data, resp$formula_env)
            if (!is.matrix(A) || nrow(A) != ncol(A)) {
              stop("gr(cov=): cov must be a square matrix", call. = FALSE)
            }
            lv <- levels(fac)
            if (is.null(rownames(A)) || !all(lv %in% rownames(A))) {
              stop("gr(cov=): cov needs dimnames covering all grouping ",
                   "levels", call. = FALSE)
            }
            aux_A <- unname(A[lv, lv])
            if (d_k > 1L) {
              nl_k0 <- len_k %/% d_k
              r <- seq_len(len_k)
              l1 <- (r - 1L) %/% d_k + 1L
              c1 <- (r - 1L) %% d_k + 1L
              aux_kron <- list(
                ia = as.vector(outer(l1, l1,
                                     function(a, b) (b - 1L) * nl_k0 + a)),
                is = as.vector(outer(c1, c1,
                                     function(a, b) (b - 1L) * d_k + a))
              )
            }
          }
          if (cs_name == "equalto") {
            V <- eval(dp$re[[k]]$cov_expr, data, resp$formula_env)
            if (!is.matrix(V) || nrow(V) != d_k || ncol(V) != d_k) {
              stop("equalto(): V must be a ", d_k, " x ", d_k,
                   " matrix", call. = FALSE)
            }
            aux_A <- unname(V)
          }
          Zk <- Matrix::t(rt$Zt[rt$Gp[k] + seq_len(len_k), , drop = FALSE])
          components[[length(components) + 1L]] <- list(
            lp_key = lp_key, dpar = dp$name, resp = resp$resp_name,
            covstruct = cs_name, id = dp$re[[k]]$id,
            rank = dp$re[[k]]$rank,
            dim = d_k, n_levels = len_k %/% d_k,
            levels = levels(fac), cnms = rt$cnms[[k]],
            bar = bars[[k]], Zlocal = Zk, aux_A = aux_A,
            aux_D = aux_D, aux_kron = aux_kron, aux_Q = aux_Q,
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

      # Monotonic terms: one scale coefficient in beta (a zero column in
      # the stored X keeps the bookkeeping - names, idx, vcov - while
      # the objective and the numeric prediction paths supply the
      # simplex-weighted values); the simplex parameters join `extras`.
      mo_info <- list()
      for (mexpr in dp$mo %||% list()) {
        v <- eval(mexpr, mf, resp$formula_env)
        if (is.factor(v)) {
          if (!is.ordered(v)) {
            stop("mo(): factor variables must be ordered factors",
                 call. = FALSE)
          }
          codes <- as.integer(v) - 1L
          D_mo <- nlevels(v) - 1L
          mo_levels <- levels(v)
        } else {
          if (any(v < 0) || any(v != round(v))) {
            stop("mo(): variable must be an ordered factor or ",
                 "non-negative integers", call. = FALSE)
          }
          codes <- as.integer(v)
          D_mo <- max(codes)
          mo_levels <- NULL
        }
        if (D_mo < 2L) {
          stop("mo() needs at least 3 ordered categories", call. = FALSE)
        }
        lab <- paste0("mo", deparse1(mexpr))
        zname <- paste0("zeta", length(extras) + 1L)
        extras[[zname]] <- numeric(D_mo - 1L)
        X <- cbind(X, matrix(0, nrow(X), 1,
                             dimnames = list(NULL, lab)))
        mo_info[[length(mo_info) + 1L]] <- list(
          expr = mexpr, codes = codes, D = D_mo, levels = mo_levels,
          zeta = zname, col = ncol(X), label = lab
        )
      }

      # mi(x) predictor terms: one coefficient (zero placeholder column,
      # as with mo); the values are observed-or-latent, supplied by the
      # objective and the numeric prediction paths
      mi_info <- list()
      for (mexpr in dp$miterms %||% list()) {
        vn <- deparse1(mexpr)
        tgt <- spec$responses[[vn]]
        if (is.null(tgt) || !isTRUE(tgt$aterms$mi)) {
          stop("mi(", vn, ") needs a matching imputation model: ",
               "add bf(", vn, " | mi() ~ ...)", call. = FALSE)
        }
        if (identical(vn, resp$resp_name)) {
          stop("mi(", vn, ") cannot appear in its own model",
               call. = FALSE)
        }
        lab <- paste0("mi", vn)
        X <- cbind(X, matrix(0, nrow(X), 1, dimnames = list(NULL, lab)))
        mi_info[[length(mi_info) + 1L]] <- list(
          var = vn, col = ncol(X), label = lab
        )
      }

      # Category-specific ordinal effects cs(x): K-1 coefficients per
      # term (extras), entering the threshold-specific predictors.
      cs_info <- list()
      if (length(dp$csterms %||% list())) {
        if (!identical(resp$family$type, "ordinal") ||
            identical(resp$family$family, "cumulative")) {
          stop("cs() needs an sratio, cratio, or acat family",
               call. = FALSE)
        }
        K_cs <- max(y[[resp$resp_name]])
        for (cexpr in dp$csterms) {
          v <- as.numeric(eval(cexpr, mf, resp$formula_env))
          csname <- paste0("bcs", length(extras) + 1L)
          extras[[csname]] <- numeric(K_cs - 1L)
          cs_info[[length(cs_info) + 1L]] <- list(
            vals = v, par = csname,
            label = paste0("cs", deparse1(cexpr))
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
        mo = mo_info,
        mi = mi_info,
        cs = cs_info,
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
  n_b <- 0L      # parameter space (rr blocks hold factors here)
  n_c <- 0L      # coefficient space (the Z columns)
  n_theta <- 0L
  has_rr <- FALSE

  for (gd in group_defs) {
    cps <- components[gd]
    if (length(cps) > 1L) {
      if (any(vapply(cps, `[[`, "", "covstruct") == "rr")) {
        stop("rr() terms cannot share an |ID| key", call. = FALSE)
      }
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
      rank_k <- NULL
    } else {
      cp <- cps[[1]]
      D <- cp$dim
      cs_name <- cp$covstruct
      cnms <- cp$cnms
      label <- cp$label
      n_levels <- cp$n_levels
      rank_k <- cp$rank
    }
    if (cs_name == "rr") {
      if (is.null(rank_k) || rank_k > D) {
        stop("rr(): the rank d must not exceed the term dimension (",
             D, ")", call. = FALSE)
      }
      has_rr <- TRUE
      npar_k <- rr_npar(D, rank_k)
      nb_k <- rank_k * n_levels
    } else {
      npar_k <- covstruct_registry[[cs_name]]$npar(D)
      nb_k <- D * n_levels
    }
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
      rank = rank_k,
      n_levels = n_levels,
      b_idx = n_b + seq_len(nb_k),
      c_idx = n_c + seq_len(D * n_levels),
      theta_idx = n_theta + seq_len(npar_k),
      levels = cps[[1]]$levels,
      aux_A = cps[[1]]$aux_A,
      aux_D = cps[[1]]$aux_D,
      aux_kron = cps[[1]]$aux_kron,
      aux_Q = cps[[1]]$aux_Q,
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
    n_b <- n_b + nb_k
    n_c <- n_c + D * n_levels
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
        glob <- bk$c_idx[(lev - 1L) * bk$dim + comp_offset[ci] + cf]
        ii <- c(ii, Tk@i + 1L)
        jj <- c(jj, glob)
        xx <- c(xx, Tk@x)
      }
      lp$Z <- Matrix::sparseMatrix(i = ii, j = jj, x = xx,
                                   dims = c(n, n_c))
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
  if (n_theta) {
    th0 <- numeric(n_theta)
    for (bk in re_blocks) {
      th0[bk$theta_idx] <- if (bk$covstruct == "rr") {
        rr_start(bk$dim, bk$rank)
      } else {
        covstruct_registry[[bk$covstruct]]$start(bk$dim)
      }
    }
    par_template$theta <- th0
  }
  if (spec$rescor) {
    K <- length(spec$responses)
    par_template$thetar <- numeric(K * (K - 1L) / 2L)
  }
  if (n_miss) par_template$miss <- miss_init
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
         n_c = n_c, has_rr = has_rr, mi_map = mi_map,
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
