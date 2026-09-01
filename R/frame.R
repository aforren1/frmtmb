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

#' Sum of offset() terms of one linear predictor, evaluated against the
#' combined model frame (whose columns are named by deparsed expressions).
#'
#' @noRd
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

#' The variables one dpar needs in the combined model frame: parametric
#' terms, bar variables, and raw smooth variables (never the s() calls
#' themselves, which model.frame cannot evaluate).
#'
#' @noRd
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
  for (ent in dp$mo %||% list()) {
    for (v in all.vars(ent$expr)) parts <- c(parts, list(as.name(v)))
    if (!is.null(ent$mult)) parts <- c(parts, list(ent$mult))
  }
  for (ent in dp$miterms %||% list()) {
    if (!is.null(ent$mult)) parts <- c(parts, list(ent$mult))
  }
  for (cexpr in dp$csterms %||% list()) {
    for (v in all.vars(cexpr)) parts <- c(parts, list(as.name(v)))
  }
  for (ge in dp$gpterms %||% list()) {
    for (ex in ge$exprs) {
      for (v in all.vars(ex)) parts <- c(parts, list(as.name(v)))
    }
  }
  # car()/spde() need their grouping variable in the frame; the
  # adjacency and mesh matrices are external data, like gr(cov = A)'s.
  # The RAW variables go in, not the expression - a call-valued gr
  # (gr = factor(node)) is evaluated against the frame at assembly, the
  # way every other special resolves its arguments, and model.frame
  # cannot be asked to carry the call itself.
  for (ce in c(dp$carterms %||% list(), dp$spdeterms %||% list())) {
    for (v in all.vars(ce$gr_expr)) parts <- c(parts, list(as.name(v)))
  }
  out <- NULL
  for (p in parts) out <- if (is.null(out)) p else call("+", out, p)
  out
}

#' mo()/mi() interaction multipliers scale a single coefficient, so they
#' have to be one numeric column; a factor or character multiplier would
#' need contrast expansion, which the simplex/latent machinery has no
#' column for. Reject those up front: as.numeric() on a character vector
#' yields all-NA and the failure only surfaces as "NA/NaN gradient
#' evaluation" from the optimizer. `[brms#1828]`
#'
#' @noRd
check_special_mult <- function(mult, expr, fn) {
  if (is.logical(mult)) return(as.numeric(mult))
  if (!is.numeric(mult) || is.factor(mult)) {
    stop(fn, "() interactions support numeric multipliers only: ",
         deparse1(expr), " is ", class(mult)[1L],
         "; expand it to numeric indicator columns first", call. = FALSE)
  }
  as.numeric(mult)
}

#' ar1()/hetar1() correlate two levels by their POSITION in the ordering
#' factor, never by the distance between the labels: with times 1..6 and
#' 10 present, cor(t6, t10) is fitted as rho, not rho^4, and rho itself
#' comes out biased. glmmTMB reads the levels the same way and our
#' agreement tests pin that reading down, so the likelihood stays as it
#' is and the user gets told instead. Silence is right when the labels
#' carry no numbers: position is then the only meaning available.
#' `[glmmTMB#1278]`
#'
#' @noRd
warn_ar1_level_gaps <- function(bar, mf, cs_name) {
  v <- all.vars(bar[[2]])
  if (length(v) != 1L || !is.factor(mf[[v]])) return(invisible(NULL))
  lv <- levels(mf[[v]])
  pos <- suppressWarnings(as.numeric(lv))
  if (anyNA(pos) || any(pos != trunc(pos))) return(invisible(NULL))
  gap <- which(abs(diff(pos)) != 1)
  if (!length(gap)) return(invisible(NULL))
  i <- gap[1L]
  warning(cs_name, "(): the levels of '", v, "' are whole numbers but ",
          "not consecutive ('", lv[i], "' is followed by '", lv[i + 1L],
          "'), and ", cs_name, "() correlates levels by position, so ",
          "that gap counts as a single step. For irregularly spaced ",
          "positions use ou() over num_factor(): ou(num_factor(", v,
          ") + 0 | ...)", call. = FALSE)
  invisible(NULL)
}

#' Pull one response out of the combined model frame and coerce it to the
#' numeric form the objective needs. Ordinal factors become category
#' codes and keep their labels in a `y_levels` attribute, binomial
#' two-level factors become 0/1, and unsupported factor responses or
#' non-finite values are rejected here.
#'
#' @noRd
extract_y <- function(resp, mf) {
  y <- mf[[deparse1(resp$resp_expr)]]
  if (is.null(y)) {
    y <- eval(resp$resp_expr, mf, resp$formula_env)
  }
  lv <- NULL
  if (is.factor(y) && identical(resp$family$type, "ordinal")) {
    # brms accepts an unordered factor here and silently reads level
    # order as category order, which is alphabetical unless the user set
    # it. Match that behavior (refusing would break working brms code)
    # but say so, because the ordering is the whole model. mo() takes the
    # stricter route because there the variable is a predictor and the
    # user can always relabel it.
    if (!is.ordered(y)) {
      warning("Ordinal response '", deparse1(resp$resp_expr),
              "' is an unordered factor; its level order (",
              paste(levels(y), collapse = " < "),
              ") is taken as the category order. Use ",
              "factor(..., ordered = TRUE) to state it explicitly",
              call. = FALSE)
    }
    # the codes carry no meaning without the labels, and simulate() has
    # to hand draws back in the response's own type
    lv <- levels(y)
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
  # attached last: the numeric coercions above drop attributes
  if (!is.null(lv)) attr(y, "y_levels") <- lv
  y
}

#' Data-dependent bases (poly, ns, scale) must be frozen at fit time: the
#' combined model frame's terms carry predvars for every variable; patch
#' them onto a sub-formula's terms by deparsed-variable match before any
#' newdata evaluation (glmmTMB#402 and the largest bug class in
#' lme4/glmmTMB prediction history).
#'
#' @noRd
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

#' sparse.model.matrix names the columns of matrix-valued terms (poly,
#' ns, scale) by bare index instead of the term-prefixed dense names.
#' Names are load-bearing (frozen param_colnames, coefficient labels), so
#' take them from a one-row dense header over the same frame.
#'
#' This is the sparse counterpart of `stats::model.matrix()` used by
#' `assemble_frame()` when `sparse_x = TRUE`. It takes the terms of one
#' linear predictor, the combined model frame, and the contrasts, and
#' returns a sparse design matrix whose column names match the dense
#' path exactly. It errors if the two paths disagree on the column
#' count, because every later stage indexes coefficients by those names.
#'
#' @noRd
sparse_mm <- function(tt, mf, contrasts.arg = NULL) {
  X <- Matrix::sparse.model.matrix(tt, mf, contrasts.arg = contrasts.arg)
  # keeping the frame's terms attribute on the one-row slice routes
  # model.matrix through deparse-name column matching, never
  # re-evaluating data-dependent bases (poly, scale) on a single row
  mf1 <- mf[1L, , drop = FALSE]
  attr(mf1, "terms") <- attr(mf, "terms")
  hdr <- stats::model.matrix(tt, mf1, contrasts.arg = contrasts.arg)
  if (ncol(hdr) != ncol(X)) {
    stop("Internal error: sparse and dense fixed-effect designs disagree ",
         "on columns; refit without frmtmb_control(sparse_x = TRUE)",
         call. = FALSE)
  }
  colnames(X) <- colnames(hdr)
  X
}

#' Cheap rank screen for a sparse design: |diag(R)| of the sparse QR
#' collapses on aliased columns. The threshold is 100x looser than dense
#' qr()'s 1e-7 default, so anything the dense path would drop gets
#' flagged; a flagged design is re-checked densely, which decides the
#' exact drop set (identical to the dense path by construction).
#'
#' @noRd
sparse_maybe_deficient <- function(X) {
  d <- tryCatch(
    abs(Matrix::diag(Matrix::qrR(Matrix::qr(X), backPermute = FALSE))),
    error = function(e) NULL
  )
  is.null(d) || !all(is.finite(d)) || min(d) <= 0 ||
    min(d) < 1e-5 * max(d)
}

# Operators that reformulas expands structurally inside a grouping
# expression; everything else on the right of a bar is an ordinary call
# whose value has to exist as a single model-frame column.
grp_structural_ops <- c(":", "/", "*", "+", "%in%")

#' A grouping factor written as a call - `(1 | factor(x))`,
#' `(1 | interaction(a, b))` - is one variable of the combined model frame,
#' stored under its deparsed name. reformulas instead re-evaluates the
#' expression inside the frame, where the call's own arguments (x, a, b)
#' are not columns, and dies with an error raised several frames down.
#' Point the bar at the existing column by name instead; the original bar
#' is what the fit keeps, so prediction still evaluates the expression
#' against newdata. `[lme4#464, #156]`
#'
#' @noRd
resolve_group_calls <- function(bars, fr, env) {
  sub_call <- function(e) {
    if (!is.call(e)) return(e)
    if (as.character(e[[1L]])[1L] %in% grp_structural_ops) {
      for (i in seq_along(e)[-1]) e[[i]] <- sub_call(e[[i]])
      return(e)
    }
    key <- deparse1(e)
    if (!key %in% names(fr)) {
      # defensive: a call nested inside ':' need not be a frame variable
      fr[[key]] <<- eval(e, fr, env)
    }
    as.name(key)
  }
  bars <- lapply(bars, function(b) {
    b[[3L]] <- sub_call(b[[3L]])
    b
  })
  list(bars = bars, fr = fr)
}

# Internal censoring codes, shared with brms: -1 left, 0 observed,
# 1 right, 2 interval. The names are the accepted spellings.
cens_code_map <- c(none = 0, left = -1, right = 1, interval = 2)

#' brms lets cens() carry spelled-out codes as well as numbers, and a
#' plain as.numeric() would turn those into NAs with only a coercion
#' warning, so decode before the numeric conversion the other addition
#' terms use. Prefix matching mirrors brms prepare_cens(); unlike brms
#' it is case-insensitive, because "Right" reads as a typo, not garbage.
#'
#' @noRd
decode_cens <- function(v) {
  if (is.factor(v)) v <- as.character(v)
  if (!is.character(v)) return(as.numeric(v))
  key <- tolower(trimws(v))
  # the numeric codes the error message advertises can arrive as
  # strings (a character column, or a factor built from one); prefix
  # matching would reject every one of them, so read those directly.
  # Out-of-range numbers fall through to the -1/0/1/2 check, which
  # names the offending value
  num <- suppressWarnings(as.numeric(key))
  idx <- vapply(key, function(k) {
    hit <- if (nzchar(k)) which(startsWith(names(cens_code_map), k)) else
      integer(0)
    if (length(hit) == 1L) hit else NA_integer_
  }, integer(1L), USE.NAMES = FALSE)
  bad <- unique(v[is.na(idx) & is.na(num)])
  if (length(bad)) {
    stop("cens() cannot decode: ",
         paste0("\"", bad, "\"", collapse = ", "),
         "; use \"none\", \"left\", \"right\", or \"interval\" ",
         "(any unambiguous prefix), or the codes 0, -1, 1, 2",
         call. = FALSE)
  }
  out <- unname(cens_code_map[idx])
  ifelse(is.na(num), out, num)
}

# data2 must be a named list; anything else is a user mistake that would
# otherwise surface far away, inside one of the structural lookups.
validate_data2 <- function(data2) {
  if (is.null(data2)) return(list())
  nms <- names(data2)
  if (!is.list(data2) || (length(data2) && is.null(nms)) ||
      any(!nzchar(nms))) {
    stop("data2 must be a named list, e.g. data2 = list(W = W)",
         call. = FALSE)
  }
  if (anyDuplicated(nms)) {
    stop("data2 has duplicate names: ",
         paste(unique(nms[duplicated(nms)]), collapse = ", "),
         call. = FALSE)
  }
  data2
}

# Structural objects (adjacency matrices, precisions, covariance
# matrices, FEM triples) resolve from data2 before anything else, so a
# fit that names them there is self-contained: saveRDS() carries them on
# the fit and a later refit never reaches back into the calling
# environment that built them, which by then may be gone.
#
# brms accepts a bare name in data2 and nothing more. We keep that rule
# and add one: a compound expression is evaluated with data2 in front of
# the data mask, so car(solve(P)) finds P in data2 too. The historical
# data-then-formula-env evaluation stays as the fallback, so models
# written before data2 existed keep working.
lookup_structural <- function(expr, data2, data, env, what) {
  e2 <- NULL
  if (length(data2)) {
    if (is.symbol(expr)) {
      nm <- as.character(expr)
      if (nm %in% names(data2)) return(data2[[nm]])
    }
    # the wrapper list separates "evaluated to NULL" from "failed"
    val <- tryCatch(
      list(v = eval(expr, data2, list2env(as.list(data), parent = env))),
      error = function(e) {
        e2 <<- e
        NULL
      }
    )
    if (!is.null(val)) return(val$v)
  }
  tryCatch(eval(expr, data, env), error = function(e) {
    # the data2-mask attempt saw the widest scope, so when both paths
    # fail its error names the real cause (the fallback just repeats
    # "not found" for objects that only exist in data2)
    stop(structural_lookup_msg(expr, data2, what,
                               if (is.null(e2)) e else e2),
         call. = FALSE)
  })
}

structural_lookup_msg <- function(expr, data2, what, e) {
  txt <- deparse1(expr)
  held <- if (length(data2)) {
    paste0(" data2 holds: ", paste(names(data2), collapse = ", "), ".")
  } else {
    " data2 is empty."
  }
  if (is.symbol(expr)) {
    paste0(what, ": cannot find '", txt, "'. Pass structural objects ",
           "in data2, e.g. frm(..., data2 = list(", txt, " = ", txt,
           ")); a fit whose matrices come from data2 also survives ",
           "saveRDS() and refits in a new session.", held)
  } else {
    paste0(what, ": cannot evaluate '", txt, "' (", conditionMessage(e),
           "). Put the objects it needs in data2, e.g. ",
           "frm(..., data2 = list(...)).", held)
  }
}

#' Turn a parsed spec plus data into the numeric `frmtmb_frame` the
#' objective is built from. This is the second and last stage of the
#' formula-to-design-matrix pipeline: `parse_spec()` reads the formulas,
#' `assemble_frame()` reads the data.
#'
#' It builds one combined model frame over every response, every dpar
#' variable and every addition-term variable, so `na.action` drops rows
#' from all of them together. It then evaluates the responses and the
#' addition terms, builds the fixed-effect design of each linear
#' predictor (dropping aliased columns the way `lm()` does), and turns
#' every random-effect, smooth, Gaussian-process, spatial, monotonic,
#' missing-data and category-specific term into a component. Components
#' sharing an `|ID|` key merge into one block, blocks are laid out in
#' the flat `b` and `theta` vectors, and each linear predictor gets a
#' sparse `Z` over the whole `b` vector.
#'
#' The returned object carries the designs (`linpreds`), the
#' random-effect blocks (`re_blocks`), the response and addition-term
#' values, the `par_template` with the MakeADFun `map`, and everything
#' prediction needs to reproduce the designs on newdata: terms,
#' xlevels, contrasts, frozen predvars, and the model frame itself.
#' Structural matrices (adjacency, precision, covariance, FEM triples)
#' resolve through `lookup_structural()`: data2 first, then data, then
#' the formula environment.
#'
#' @noRd
assemble_frame <- function(spec, data, na.action = stats::na.omit,
                           sparse_x = FALSE, data2 = list()) {
  data2 <- validate_data2(data2)
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
    if (!is.null(resp$family$mix_groups)) {
      add_part(resp$family$mix_groups[[2L]])
    }
    for (nm_at in names(resp$aterms)) {
      # literal constants (e.g. trials(10)) are not frame variables, and
      # interval bounds (cens_y2) may be NA on non-interval rows, so they
      # stay out of the na.omit frame (brms#1070); a signed literal
      # (trunc(lb = -5)) parses as a unary call, not a numeric, and
      # would break the frame formula
      a <- resp$aterms[[nm_at]]
      signed_literal <- is.call(a) && length(a) == 2L &&
        as.character(a[[1]])[1] %in% c("-", "+") && is.numeric(a[[2]])
      if (nm_at != "cens_y2" && !is.numeric(a) && !is.logical(a) &&
          !signed_literal) {
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
  # Dropping rows changes the estimand and the n every later standard
  # error is built on, so the loss is reported instead of being inferred
  # from nobs(). One message per fit; suppressMessages() silences it.
  n_dropped <- length(attr(mf, "na.action"))
  if (n_dropped > 0L) {
    message(n_dropped, if (n_dropped == 1L) " row" else " rows",
            " removed because of missing values (na.action)")
  }
  n <- nrow(mf)
  if (n == 0L) {
    # zero rows in and zero rows left are different faults, and the
    # generic "after removing NAs" wording sends the second one hunting
    # for missing values that were never there
    stop(if (n_dropped > 0L) {
           "No complete observations after removing NAs"
         } else {
           "`data` has no rows; nothing to fit"
         }, call. = FALSE)
  }
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
  y_levels <- list()
  aterm_values <- list()
  extras <- list()
  mi_map <- list()   # per mi() response: missing rows + miss indices
  n_miss <- 0L
  miss_init <- numeric(0)
  mix_g <- list()    # per response: latent-class grouping structure
  for (resp in spec$responses) {
    # A name that is both a nonlinear parameter and a data column is
    # ambiguous, and the nonlinear body resolves it to the PARAMETER,
    # silently ignoring the column - the fit runs and reports numbers
    # for a model the user did not write. Refuse rather than document a
    # precedence nobody would remember. [brms#391, #734]
    if (length(resp$nlpars)) {
      clash <- intersect(resp$nlpars, names(data))
      if (length(clash)) {
        stop("Nonlinear parameter(s) ",
             paste0("'", clash, "'", collapse = ", "),
             " also name columns of the data. The nonlinear formula ",
             "would use the parameter and ignore the column; rename ",
             "one of them", call. = FALSE)
      }
    }
    yv0 <- extract_y(resp, mf)
    y_levels[[resp$resp_name]] <- attr(yv0, "y_levels")
    attr(yv0, "y_levels") <- NULL   # nothing on the tape carries labels
    y[[resp$resp_name]] <- yv0
    at_names <- setdiff(names(resp$aterms),
                        c("cens_y2", "se_sigma", "mi"))
    av <- stats::setNames(lapply(at_names, function(nm_at) {
      a <- resp$aterms[[nm_at]]
      v <- mf[[deparse1(a)]]
      if (is.null(v)) v <- eval(a, mf, resp$formula_env)
      if (nm_at == "cens") decode_cens(v) else as.numeric(v)
    }), at_names)
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
    # the joint-gaussian rescor likelihood has no censoring,
    # truncation, or known-se machinery; without this guard those
    # terms were silently dropped (wrong likelihood, no warning)
    if (spec$rescor &&
        (!is.null(av$cens) || !is.null(av$trunc_lb) ||
         !is.null(av$trunc_ub) || !is.null(av$se))) {
      stop("cens()/trunc()/se() cannot be combined with rescor = TRUE",
           call. = FALSE)
    }
    if (!is.null(av$cens) || !is.null(av$trunc_lb) ||
        !is.null(av$trunc_ub)) {
      if (is.null(resp$family$lcdf)) {
        stop("cens()/trunc() need a family with a CDF (currently: ",
             "gaussian, lognormal, poisson, exponential, weibull, ",
             "inverse.gaussian)", call. = FALSE)
      }
      if (!is.null(av$cens) && identical(resp$family$type, "discrete")) {
        stop("cens() is not supported for discrete families yet ",
             "(truncation is)", call. = FALSE)
      }
      if (!is.null(av$cens) && !all(av$cens %in% c(-1, 0, 1, 2))) {
        stop("cens() codes must be -1 (left), 0 (observed), 1 (right), ",
             "or 2 (interval), or the matching names \"left\", \"none\", ",
             "\"right\", \"interval\"; got: ",
             paste(unique(av$cens[!av$cens %in% c(-1, 0, 1, 2)]),
                   collapse = ", "), call. = FALSE)
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
    # cbind(successes, failures) reaches here already rewritten to
    # successes + trials(successes + failures); a fractional failure
    # count would otherwise surface as the generic "response must be
    # integer counts in [0, trials]" message, which names neither
    # column. [glmmTMB#1319]
    if (isTRUE(resp$cbind_resp)) {
      fails <- av$trials - y[[resp$resp_name]]
      if (any(fails < 0) || any(fails != round(fails))) {
        stop("cbind(successes, failures): the failure column must hold ",
             "non-negative integer counts", call. = FALSE)
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
    if (!is.null(resp$family$mix_groups)) {
      if (length(spec$responses) > 1L || spec$rescor) {
        stop("Group-level mixtures support univariate models",
             call. = FALSE)
      }
      if (any(c("cens", "trunc_lb", "trunc_ub") %in%
                names(resp$aterms))) {
        stop("Group-level mixtures cannot be combined with cens() or ",
             "trunc()", call. = FALSE)
      }
      if (isTRUE(resp$aterms$mi)) {
        stop("Group-level mixtures cannot be combined with mi() on ",
             "the same response", call. = FALSE)
      }
      gv <- factor(eval(resp$family$mix_groups[[2L]], mf,
                        resp$formula_env))
      G <- Matrix::sparseMatrix(i = seq_len(n), j = as.integer(gv),
                                x = 1, dims = c(n, nlevels(gv)))
      mix_g[[resp$resp_name]] <- list(G = G, Gt = Matrix::t(G),
                                      first = match(levels(gv), gv),
                                      gindex = as.integer(gv),
                                      levels = levels(gv))
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
      X <- if (sparse_x) sparse_mm(tt, mf) else stats::model.matrix(tt, mf)
      contr <- attr(X, "contrasts")   # subsetting X below drops the attr
      if (isTRUE(resp$family$drop_intercept) && is_primary) {
        # thresholds replace the intercept; a threshold-only model
        # (y ~ 1) leaves X with zero columns, which is fine
        X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
      }
      # rank-deficient designs: drop aliased columns like lm() (lme4#144).
      # Sparse X densifies a copy only when the cheap screen flags a
      # possible deficiency, so the dropped-column set never differs
      # from the dense path.
      alias_null <- NULL
      dropped_colnames <- NULL
      if (ncol(X) > 1L) {
        Xq <- if (!sparse_x) X
              else if (sparse_maybe_deficient(X)) as.matrix(X)
              else NULL
        if (!is.null(Xq)) {
          Xq <- as.matrix(Xq)
          qrX <- qr(Xq)
          if (qrX$rank < ncol(X)) {
            dropped <- colnames(X)[qrX$pivot[(qrX$rank + 1L):ncol(X)]]
            message("Fixed-effect design of '", lp_key,
                    "' is rank deficient; dropping column(s): ",
                    paste(dropped, collapse = ", "))
            # Directions the data could not identify: null(X) is the
            # orthogonal complement of the row space, so the trailing
            # columns of the complete Q of t(X) span it. Frozen here so
            # prediction can refuse rows that load on them. [lme4#303]
            alias_null <- qr.Q(qr(t(Xq)), complete = TRUE)[
              , (qrX$rank + 1L):ncol(Xq), drop = FALSE]
            rownames(alias_null) <- colnames(X)
            dropped_colnames <- dropped
            X <- X[, setdiff(colnames(X), dropped), drop = FALSE]
          }
        }
      }
      param_colnames <- colnames(X)
      n_param_cols <- ncol(X)
      off <- extract_offset(tt, mf, resp$formula_env)
      xlev <- stats::.getXlevels(stats::terms(mf), mf)
      comp_ids <- integer(0)

      if (length(dp$re)) {
        bars <- lapply(dp$re, `[[`, "bar")
        # bars keeps the user's expressions (labels, prediction); only
        # the copy handed to reformulas is name-resolved
        grp <- resolve_group_calls(bars, mf, resp$formula_env)
        rt <- reformulas::mkReTrms(grp$bars, fr = grp$fr,
                                   reorder.terms = FALSE)
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
          if (cs_name %in% c("ar1", "hetar1")) {
            warn_ar1_level_gaps(bars[[k]], mf, cs_name)
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
          aux_Qk <- NULL
          if (cs_name == "gr_prec") {
            Q <- lookup_structural(dp$re[[k]]$cov_expr, data2, data,
                                   resp$formula_env, "gr(prec = )")
            if (is.null(rownames(Q)) ||
                !all(levels(fac) %in% rownames(Q))) {
              stop("gr(prec=): prec needs dimnames covering all ",
                   "grouping levels", call. = FALSE)
            }
            lv <- levels(fac)
            aux_Q <- methods::as(Matrix::Matrix(Q[lv, lv], sparse = TRUE),
                                 "generalMatrix")
            if (d_k > 1L) aux_Qk <- kron_prec_parts(aux_Q, d_k)
          }
          if (cs_name == "gr_cov") {
            A <- lookup_structural(dp$re[[k]]$cov_expr, data2, data,
                                   resp$formula_env, "gr(cov = )")
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
            V <- lookup_structural(dp$re[[k]]$cov_expr, data2, data,
                                   resp$formula_env, "equalto()")
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
            aux_Qk = aux_Qk,
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

      # Gaussian-process terms: gp(x1, ...) is exact (a dense SE-kernel
      # block over the unique coordinate rows); gp(..., k=) is the
      # Hilbert-space approximation (tensor-product sine basis in Z,
      # spectral-density prior SDs). D up to 3; iso shares one
      # lengthscale, the default is one per dimension (brms).
      gp_info <- list()
      for (ge in dp$gpterms %||% list()) {
        Xc <- do.call(cbind, lapply(ge$exprs, function(ex) {
          as.numeric(eval(ex, mf, resp$formula_env))
        }))
        Dg <- ncol(Xc)
        iso <- isTRUE(ge$iso) || Dg == 1L
        vnames <- vapply(ge$exprs, deparse1, "")
        lab0 <- paste0("gp(", paste(vnames, collapse = ", "),
                       if (!is.null(ge$k)) paste0(", k = ", ge$k) else "",
                       ")")
        if (is.null(ge$k)) {
          posdf <- unique(as.data.frame(Xc))
          posdf <- posdf[do.call(order, posdf), , drop = FALSE]
          pos <- unname(as.matrix(posdf))
          npos <- nrow(pos)
          if (npos > 500L) {
            stop("gp() without k= builds a dense ", npos,
                 "-point covariance; use k= for the Hilbert-space ",
                 "approximation", call. = FALSE)
          }
          Zg <- Matrix::sparseMatrix(i = seq_len(nrow(Xc)),
                                     j = match(pos_rowkey(Xc),
                                               pos_rowkey(pos)),
                                     x = 1, dims = c(nrow(Xc), npos))
          components[[length(components) + 1L]] <- list(
            lp_key = lp_key, dpar = dp$name, resp = resp$resp_name,
            covstruct = "gp", id = NULL,
            dim = npos, n_levels = 1L,
            levels = NULL, cnms = paste0(lab0, ".", seq_len(npos)),
            bar = NULL, Zlocal = methods::as(Zg, "CsparseMatrix"),
            aux_D2 = lapply(seq_len(Dg), function(j) {
              outer(pos[, j], pos[, j], "-")^2
            }),
            gp_D = Dg, gp_iso = iso, gp_vars = vnames,
            group_name = lab0,
            label = paste0(dp_prefix, lab0)
          )
          gp_info[[length(gp_info) + 1L]] <- list(
            exprs = ge$exprs, type = "exact", positions = pos,
            comp_id = length(components), block_id = NULL, label = lab0
          )
        } else {
          # brms's input convention (brms:::.data_gp): rescale by the
          # largest pairwise distance over the DISTINCT coordinate rows,
          # center on that scale, then take a shared boundary
          # L_j = c_j * max(1, range of the whole centered matrix). The
          # same gp(x, k, c) call is then the same approximation here and
          # in brms. Distinct rows because brms's gr = TRUE default
          # collapses duplicate positions before computing the scale, so
          # ties would otherwise shift the center.
          cvec <- ge$c
          if (length(cvec) == 1L) cvec <- rep(cvec, Dg)
          if (length(cvec) != Dg) {
            stop("gp(): c = must be length 1 or the number of ",
                 "variables (", Dg, ")", call. = FALSE)
          }
          uq <- Xc[!duplicated(pos_rowkey(Xc)), , drop = FALSE]
          dmax <- gp_max_dist(uq)
          if (!isTRUE(dmax > 0)) {
            # a single scale over all coordinates, so it vanishes only
            # when every coordinate row is identical
            stop("gp(", paste(vnames, collapse = ", "),
                 "): the coordinates have no spread", call. = FALSE)
          }
          ctr <- colMeans(uq / dmax)
          Lb <- gp_choose_L(sweep(uq / dmax, 2, ctr), cvec)
          xc <- sweep(Xc / dmax, 2, ctr)
          m <- ge$k
          if (m^Dg > 1000) {
            stop("gp(): k = ", m, " over ", Dg, " dimensions gives ",
                 m^Dg, " basis columns (cap 1000); lower k=",
                 call. = FALSE)
          }
          idx <- as.matrix(do.call(expand.grid,
                                   rep(list(seq_len(m)), Dg)))
          omega <- sweep(idx * pi, 2, 2 * Lb, "/")
          Phi <- hsgp_basis(xc, omega, Lb)
          M_b <- nrow(omega)
          components[[length(components) + 1L]] <- list(
            lp_key = lp_key, dpar = dp$name, resp = resp$resp_name,
            covstruct = "hsgp", id = NULL,
            dim = M_b, n_levels = 1L,
            levels = NULL, cnms = paste0(lab0, ".", seq_len(M_b)),
            bar = NULL, Zlocal = methods::as(Phi, "CsparseMatrix"),
            aux_omega = omega,
            gp_D = Dg, gp_iso = iso, gp_vars = vnames,
            gp_dmax = dmax,
            group_name = lab0,
            label = paste0(dp_prefix, lab0)
          )
          gp_info[[length(gp_info) + 1L]] <- list(
            exprs = ge$exprs, type = "hsgp", center = ctr, L = Lb,
            dmax = dmax, omega = omega, comp_id = length(components),
            block_id = NULL, label = lab0
          )
        }
        comp_ids <- c(comp_ids, length(components))
      }

      # Spatial GMRF terms: car(M, gr = g) over an adjacency matrix and
      # spde(fem, gr = node) over a mesh. Both are one intercept per
      # location, so the block looks like (1 | g) with a structured
      # precision; a synthetic bar keeps the prediction, ranef() and
      # VarCorr() paths unchanged.
      for (ce in c(dp$carterms %||% list(), dp$spdeterms %||% list())) {
        is_car <- !is.null(ce$M_expr)
        fn <- if (is_car) "car" else "spde"
        gv <- eval(ce$gr_expr, mf, resp$formula_env)
        if (anyNA(gv)) {
          stop(fn, "(): the grouping variable '", deparse1(ce$gr_expr),
               "' has missing values", call. = FALSE)
        }
        aux_car <- NULL
        aux_spde <- NULL
        if (is_car) {
          # the adjacency matrix names its locations, so the block's
          # level order is the data's and the matrix is permuted to it
          locs <- if (is.factor(gv)) levels(droplevels(gv)) else {
            sort(unique(as.character(gv)))
          }
          j_loc <- match(as.character(gv), locs)
          M <- lookup_structural(ce$M_expr, data2, data,
                                 resp$formula_env, "car()")
          W <- car_adjacency(M, locs)
          aux_car <- car_aux(W, ce$type, ce$con_sd)
        } else {
          # the mesh names nothing, so the block's levels ARE the mesh
          # rows and the data is permuted to them: every node gets a
          # column, observed or not (an unobserved one keeps its prior)
          fem <- lookup_structural(ce$fem_expr, data2, data,
                                   resp$formula_env, "spde()")
          aux_spde <- spde_matrices(fem)
          n_node <- nrow(aux_spde[["M0"]])
          j_loc <- spde_node_index(gv, n_node, ce$gr_expr)
          locs <- as.character(seq_len(n_node))
        }
        Zc <- Matrix::sparseMatrix(i = seq_along(j_loc), j = j_loc, x = 1,
                                   dims = c(length(j_loc), length(locs)))
        components[[length(components) + 1L]] <- list(
          lp_key = lp_key, dpar = dp$name, resp = resp$resp_name,
          covstruct = fn, id = NULL,
          dim = 1L, n_levels = length(locs),
          levels = locs, cnms = "(Intercept)",
          bar = call("|", 1, ce$gr_expr),
          Zlocal = methods::as(Zc, "CsparseMatrix"),
          aux_car = aux_car, aux_spde = aux_spde,
          car_type = ce$type,
          group_name = deparse1(ce$gr_expr),
          label = paste0(dp_prefix, ce$label)
        )
        comp_ids <- c(comp_ids, length(components))
      }

      # Monotonic terms: one scale coefficient in beta (a zero column in
      # the stored X keeps the bookkeeping - names, idx, vcov - while
      # the objective and the numeric prediction paths supply the
      # simplex-weighted values); the simplex parameters join `extras`.
      mo_info <- list()
      mo_zetas <- list()   # simplexes are shared per mo() variable
      for (ent in dp$mo %||% list()) {
        mexpr <- ent$expr
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
        vkey <- deparse1(mexpr)
        zname <- mo_zetas[[vkey]]
        if (is.null(zname)) {
          zname <- paste0("zeta", length(extras) + 1L)
          extras[[zname]] <- numeric(D_mo - 1L)
          mo_zetas[[vkey]] <- zname
        }
        mult <- NULL
        if (!is.null(ent$mult)) {
          mult <- check_special_mult(eval(ent$mult, mf, resp$formula_env),
                                     ent$mult, "mo")
        }
        lab <- paste0("mo", vkey,
                      if (!is.null(ent$mult)) {
                        paste0(":", deparse1(ent$mult))
                      } else "")
        X <- cbind(X, matrix(0, nrow(X), 1,
                             dimnames = list(NULL, lab)))
        mo_info[[length(mo_info) + 1L]] <- list(
          expr = mexpr, codes = codes, D = D_mo, levels = mo_levels,
          zeta = zname, col = ncol(X), label = lab,
          mult = mult, mult_expr = ent$mult
        )
      }

      # mi(x) predictor terms: one coefficient (zero placeholder column,
      # as with mo); the values are observed-or-latent, supplied by the
      # objective and the numeric prediction paths
      mi_info <- list()
      for (ent in dp$miterms %||% list()) {
        vn <- deparse1(ent$expr)
        tgt <- spec$responses[[vn]]
        if (is.null(tgt) || !isTRUE(tgt$aterms$mi)) {
          stop("mi(", vn, ") needs a matching imputation model: ",
               "add bf(", vn, " | mi() ~ ...)", call. = FALSE)
        }
        if (identical(vn, resp$resp_name)) {
          stop("mi(", vn, ") cannot appear in its own model",
               call. = FALSE)
        }
        mult <- NULL
        if (!is.null(ent$mult)) {
          mult <- check_special_mult(eval(ent$mult, mf, resp$formula_env),
                                     ent$mult, "mi")
        }
        lab <- paste0("mi", vn,
                      if (!is.null(ent$mult)) {
                        paste0(":", deparse1(ent$mult))
                      } else "")
        X <- cbind(X, matrix(0, nrow(X), 1, dimnames = list(NULL, lab)))
        mi_info[[length(mi_info) + 1L]] <- list(
          var = vn, col = ncol(X), label = lab,
          mult = mult, mult_expr = ent$mult
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
        alias_null = alias_null,
        dropped_colnames = dropped_colnames,
        Z = NULL,               # filled in phase 3
        par = par_name,
        idx = idx,
        offset = if (!is.null(off)) as.numeric(off),
        link = dp$link,
        terms = tt,
        xlevels = xlev,
        contrasts = contr,
        smooths = sm_info,
        gps = gp_info,
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
    } else if (cs_name %in% c("gp", "hsgp")) {
      # parameter count depends on the gp dimension count, not the
      # block dimension (positions / basis size)
      npar_k <- gp_npar(cps[[1]]$gp_D, cps[[1]]$gp_iso)
      nb_k <- D * n_levels
    } else if (cs_name == "car") {
      # the CAR type decides whether there is a mixing parameter
      npar_k <- car_npar(cps[[1]]$car_type)
      nb_k <- D * n_levels
    } else if (cs_name == "spde") {
      npar_k <- spde_npar()
      nb_k <- D * n_levels
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
      aux_D2 = cps[[1]]$aux_D2,
      aux_kron = cps[[1]]$aux_kron,
      aux_Q = cps[[1]]$aux_Q,
      aux_Qk = cps[[1]]$aux_Qk,
      aux_car = cps[[1]]$aux_car,
      aux_spde = cps[[1]]$aux_spde,
      car_type = cps[[1]]$car_type,
      aux_omega = cps[[1]]$aux_omega,
      gp_D = cps[[1]]$gp_D,
      gp_iso = cps[[1]]$gp_iso,
      gp_vars = cps[[1]]$gp_vars,
      gp_dmax = cps[[1]]$gp_dmax,
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
    if (length(lp$gps)) {
      lp$gps <- lapply(lp$gps, function(gi) {
        gi$block_id <- comp_block[gi$comp_id]
        gi
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
      } else if (bk$covstruct == "hsgp") {
        hsgp_start(bk$gp_D, bk$gp_iso)
      } else if (bk$covstruct == "gp") {
        gp_start(bk$gp_D, bk$gp_iso)
      } else if (bk$covstruct == "car") {
        car_start(bk$car_type)
      } else if (bk$covstruct == "spde") {
        spde_start()
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
    list(spec = spec, n_obs = n, y = y, y_levels = y_levels,
         aterm_values = aterm_values,
         linpreds = linpreds, re_blocks = re_blocks,
         n_c = n_c, has_rr = has_rr, mi_map = mi_map, mix_g = mix_g,
         par_template = par_template, map = map,
         betad_fixed_idx = betad_fixed_idx,
         extra_names = names(extras),
         predvar_map = predvar_map,
         sparse_x = isTRUE(sparse_x),
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
