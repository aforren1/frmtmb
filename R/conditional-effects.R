# Conditional-effects displays, diagnostic plot method, pp_check.

#' Addition-term values for the conditional-effects grid. A grid row is
#' an artificial observation, so an aterm that changes the predictive
#' distribution (trials, se, truncation bounds) must not be taken at a
#' reference value: the mean number of trials is rarely a whole number,
#' and a mean truncation bound is nobody's bound. Those terms are read
#' only from variables the user pinned in `conditions`; literal bounds
#' apply as written. Everything else (`vint`/`vreal` payloads a custom
#' family needs) is evaluated against the grid when it can be.
#'
#' @noRd
ce_aterms <- function(rspec, nd, cset, n) {
  skip <- c("cens", "cens_y2", "se_sigma", "mi", "mi_sd", "weights")
  strict <- c("trials", "se", "trunc_lb", "trunc_ub")
  av <- list()
  for (nm in setdiff(names(rspec$aterms), skip)) {
    ex <- rspec$aterms[[nm]]
    vars <- all.vars(ex)
    pinned <- !length(vars) || all(vars %in% names(cset))
    if (nm %in% strict && !pinned) {
      stop("conditional_effects() cannot evaluate ",
           aterm_label(nm, ex), " on the effect grid: its value would ",
           "be a reference value, not a real one. Pin ",
           paste(setdiff(vars, names(cset)), collapse = ", "),
           " in conditions = list(...).", call. = FALSE)
    }
    v <- tryCatch(as.numeric(eval(ex, nd, rspec$formula_env)),
                  error = function(e) NULL)
    if (!is.null(v) && !length(v) %in% c(1L, n)) v <- NULL
    if (is.null(v) && nm %in% strict) {
      stop("conditional_effects() could not evaluate ",
           aterm_label(nm, ex), " on the effect grid",
           call. = FALSE)
    }
    if (!is.null(v)) av[[nm]] <- v
  }
  if (!is.null(rspec$aterms[["se_sigma"]]))
    av[["se_sigma"]] <- rspec$aterms[["se_sigma"]]
  av
}

#' The `trials()` variables a grid holds at 1 unless the user pins
#' them. A grid row is ONE artificial observation, so brms sets the
#' number of trials to 1 and says so; the mean number of trials is
#' rarely a whole number, and an expected count over 7.4 trials is not
#' a quantity anybody asked for.
#'
#' @noRd
ce_trial_vars <- function(rspec, base) {
  ex <- rspec$aterms[["trials"]]
  if (is.null(ex)) return(character(0))
  intersect(all.vars(ex), names(base))
}

#' Reference value a predictor is held at when it is not varied.
#'
#' Missing values are dropped rather than propagated: an `mi()` fit
#' models a predictor that HAS gaps, and a reference value of `NA` is a
#' grid of nothing. The complete cases are what it summarizes.
#'
#' @noRd
ce_ref_value <- function(col) {
  if (is.matrix(col)) {
    matrix(colMeans(col, na.rm = TRUE), 1, ncol(col))
  } else if (is.factor(col)) {
    factor(levels(col)[1L], levels = levels(col))
  } else if (is.numeric(col)) {
    mean(col, na.rm = TRUE)
  } else if (is.logical(col)) {
    FALSE
  } else {
    sort(unique(col))[1L]
  }
}

#' Grid of values for the varied (first) predictor.
#'
#' `int_cond` is the user's `int_conditions` entry for this variable and
#' replaces the grid outright. `stepwise` marks a variable the model
#' only defines at whole steps - a monotonic `mo()` predictor, whose
#' simplex assigns one increment per LEVEL and nothing in between - and
#' takes brms's rule for it, one grid point per step over the observed
#' range.
#'
#' @noRd
ce_grid_values <- function(col, resolution, nm = "the predictor",
                           int_cond = NULL, stepwise = FALSE) {
  if (!is.null(int_cond)) return(ce_int_cond(int_cond, col))
  if (is.factor(col)) {
    factor(levels(col), levels = levels(col))
  } else if (is.numeric(col)) {
    if (!any(is.finite(col))) {
      stop("Variable '", nm, "' has no finite values to build an effect ",
           "grid from", call. = FALSE)
    }
    lo <- min(col, na.rm = TRUE)
    hi <- max(col, na.rm = TRUE)
    if (stepwise) seq(lo, hi, by = 1) else seq(lo, hi, length.out = resolution)
  } else if (is.logical(col)) {
    c(FALSE, TRUE)
  } else {
    sort(unique(col))
  }
}

#' Values for the second predictor of an `"x:z"` effect.
#'
#' The value is EXACT: `mean +/- sd` as it comes, not `signif(, 3)`.
#' Rounding belongs in the label `effect2__` carries, where brms puts
#' it; rounding the value itself evaluates the model at a covariate the
#' user did not ask for, and how wrong that is depends on the
#' coefficient rather than on anything visible in the display.
#'
#' @noRd
ce_second_values <- function(col, int_cond = NULL) {
  if (!is.null(int_cond)) return(ce_int_cond(int_cond, col))
  if (is.numeric(col) && !is.matrix(col)) {
    mean(col, na.rm = TRUE) + c(-1, 0, 1) * stats::sd(col, na.rm = TRUE)
  } else {
    ce_grid_values(col, resolution = 0)
  }
}

#' One `int_conditions` entry, resolved against the column it replaces.
#'
#' brms takes either the values themselves or a function of the observed
#' column (`quantile`, say), and sorts numeric values so the grid is
#' ordered whatever the user passed. Names are kept: they become the
#' `effect2__` labels.
#'
#' @noRd
ce_int_cond <- function(int_cond, col) {
  v <- if (is.function(int_cond)) int_cond(col) else int_cond
  if (!length(v)) {
    stop("int_conditions must give at least one value per variable",
         call. = FALSE)
  }
  if (is.numeric(v) && !is.matrix(v)) {
    v <- sort(v)
  } else if (is.character(v) && is.factor(col)) {
    v <- factor(v, levels = levels(col))
  }
  v
}

#' Plottable variables of ONE linear predictor: its fixed-effect terms,
#' its smooth terms, its monotonic terms and its `mi()` terms. An `mo()`
#' or `mi()` column is a placeholder in the design matrix and its
#' variable never reaches `terms`, so what is stored with the term is
#' the only place the name survives.
#'
#' @noRd
ce_lp_vars <- function(lp) {
  v <- if (is.null(lp[["terms"]])) {
    character(0)
  } else {
    all.vars(stats::delete.response(lp[["terms"]]))
  }
  # a factor-smooth's own grouping factor is not a predictor to display:
  # the curve is drawn at the population level, which drops that term,
  # so varying its levels would draw the same line several times. It
  # stays plottable when the model also has it as a fixed effect, since
  # the terms object above contributes it then.
  for (si in lp[["smooths"]] %||% list()) {
    v <- c(v, setdiff(si$sm$term, si$group_var %||% character(0)))
  }
  for (m in lp[["mo"]] %||% list()) {
    v <- c(v, all.vars(m$expr), all.vars(m$mult_expr))
  }
  for (m in lp[["mi"]] %||% list()) {
    v <- c(v, m$var, all.vars(m$mult_expr))
  }
  unique(v)
}

#' Every variable `conditional_effects()` can vary for one display.
#'
#' A nonlinear predictor has no fixed-effect terms of its own: its
#' covariates are read straight out of the data by the nl body, and the
#' rest of the model lives in the nonlinear parameters' own predictors.
#' Enumerating the `mu` terms alone therefore found nothing to plot on
#' exactly the fits whose display is most wanted.
#'
#' The walk follows the parameters a body names rather than the whole
#' `nlpars` set, and recurses, so a chain of `nlf()` formulas reaches
#' the covariates at its far end and a nonlinear `sigma` does not pull
#' in the parameters of an unrelated nonlinear `mu`.
#'
#' @noRd
ce_plot_vars <- function(x, rspec, lp, resp, seen = character(0)) {
  v <- ce_lp_vars(lp)
  if (!is.null(lp[["nl_body"]])) {
    v <- c(v, names(lp[["data_list"]]),
           setdiff(all.vars(lp[["nl_body"]]), rspec$nlpars))
    reach <- c(lp[["nl_pars"]] %||% rspec$nlpars, lp[["nl_dpar_refs"]])
    for (np in setdiff(reach, seen)) {
      lpn <- x$frame[["linpreds"]][[linpred_key(resp, np)]]
      if (!is.null(lpn)) {
        v <- c(v, ce_plot_vars(x, rspec, lpn, resp, c(seen, np)))
      }
    }
  }
  unique(v)
}

#' The variables a monotonic term reads, walked the way ce_plot_vars()
#' walks the plottable ones so a nonlinear chain reaches the `mo()`
#' terms at its far end. These are the variables whose grid is stepwise:
#' the model is defined at their levels and nowhere between.
#'
#' @noRd
ce_step_vars <- function(x, rspec, lp, resp, seen = character(0)) {
  v <- unlist(lapply(lp[["mo"]] %||% list(), function(m) all.vars(m$expr)))
  if (!is.null(lp[["nl_body"]])) {
    reach <- c(lp[["nl_pars"]] %||% rspec$nlpars, lp[["nl_dpar_refs"]])
    for (np in setdiff(reach, seen)) {
      lpn <- x$frame[["linpreds"]][[linpred_key(resp, np)]]
      if (!is.null(lpn)) {
        v <- c(v, ce_step_vars(x, rspec, lpn, resp, c(seen, np)))
      }
    }
  }
  unique(v)
}

#' Every plottable variable of every linear predictor of one response.
#'
#' The fallback for a model whose SELECTED predictor has nothing to
#' plot: `bf(y ~ 1, theta1 ~ x) + mixture(...)` has a covariate, it is
#' just not on `mu1`, and refusing to draw it names the one predictor
#' the search looked at rather than the model the user fitted. Only
#' reached when the selected predictor is empty, so a model whose
#' `dpar` does have terms keeps enumerating that dpar's terms alone
#' (which is the deliberate difference from brms recorded as finding 6).
#'
#' @noRd
ce_plot_vars_any <- function(x, rspec, resp) {
  v <- character(0)
  for (dp in names(rspec$dpars)) {
    lpn <- x$frame[["linpreds"]][[linpred_key(resp, dp)]]
    if (!is.null(lpn)) v <- c(v, ce_plot_vars(x, rspec, lpn, resp))
  }
  unique(v)
}

#' Fitted interactions of ONE linear predictor as `"a:b"` effect
#' pairs. brms plots interaction displays by default alongside the main
#' effects, and a display that hides a fitted interaction invites
#' reading the main-effect curves as the whole story. The display
#' itself takes two variables, so a term of order three or more
#' contributes the pair of its leading two variables (the rest sit at
#' reference values, and `conditions =` pins them elsewhere). Term
#' labels are parsed back to variables so `log(x):f` still yields the
#' `x:f` pair; a component that is not one variable (a matrix column)
#' drops its term from the default.
#'
#' @noRd
ce_lp_pairs <- function(lp) {
  tt <- lp[["terms"]]
  if (is.null(tt)) return(character(0))
  ord <- attr(tt, "order")
  fac <- attr(tt, "factors")
  if (is.null(ord) || !any(ord >= 2L)) return(character(0))
  out <- character(0)
  for (k in which(ord >= 2L)) {
    lbls <- rownames(fac)[fac[, k] > 0]
    vs <- lapply(lbls, function(l) all.vars(str2lang(l)))
    if (!all(lengths(vs) == 1L)) next
    vv <- unique(unlist(vs))
    if (length(vv) >= 2L) {
      out <- c(out, paste(vv[1L], vv[2L], sep = ":"))
    }
  }
  unique(out)
}

#' The default interaction displays, walking nonlinear parameters the
#' way ce_plot_vars() does so a chain of nlf() formulas contributes the
#' interactions of the predictors at its far end.
#'
#' @noRd
ce_plot_pairs <- function(x, rspec, lp, resp, seen = character(0)) {
  p <- ce_lp_pairs(lp)
  if (!is.null(lp[["nl_body"]])) {
    reach <- c(lp[["nl_pars"]] %||% rspec$nlpars, lp[["nl_dpar_refs"]])
    for (np in setdiff(reach, seen)) {
      lpn <- x$frame[["linpreds"]][[linpred_key(resp, np)]]
      if (!is.null(lpn)) {
        p <- c(p, ce_plot_pairs(x, rspec, lpn, resp, c(seen, np)))
      }
    }
  }
  unique(p)
}

#' The columns brms's `conditional_effects()` frame carries, in brms's
#' order: the varied predictor(s) first, then every other model variable
#' at the value it is held at, then `cond__`, then the `effect1__` /
#' `effect2__` copies its `plot()` reads. The band columns are appended
#' by the caller.
#'
#' Carrying them is not decoration. brms's own `plot()` facets on
#' `cond__`, so a ported faceting call has nothing to facet on without
#' it, and the held values are the only record of WHERE the other
#' covariates sat while this one moved.
#'
#' @noRd
ce_frame <- function(nd, ev, v2 = NULL, cond = NULL, cats = NULL) {
  d <- nd[c(ev, setdiff(names(nd), ev))]
  if (!is.null(cond)) d[["cond__"]] <- cond
  if (!is.null(cats)) d[["cats__"]] <- cats
  d[["effect1__"]] <- nd[[ev[1L]]]
  e2 <- if (!is.null(cats)) {
    cats
  } else if (length(ev) == 2L) {
    nd[[ev[2L]]]
  }
  if (!is.null(e2)) d[["effect2__"]] <- ce_effect2(e2, v2)
  d
}

#' The moderator's DISPLAY label. brms rounds a numeric moderator to two
#' decimals here and only here, and orders the levels descending so a
#' legend reads down the plot. Names on an `int_conditions` entry
#' replace the numbers outright, which is how a moderator gets labels
#' like "low"/"high".
#'
#' The factor's LEVELS come from the distinct values, not from the
#' rounded ones, and the rounding is only the label text. Building the
#' levels from `round(v, 2)` gave two curves one level whenever two
#' values rounded together. `plot()` groups on this column, so it drew
#' them as one series, and a user's own names were dropped because the
#' name-count guard no longer matched. `int_conditions` makes that
#' easy to reach (`c(lo = 0.001, hi = 0.002)`); the default
#' `mean +/- sd` rarely collides. Two decimals are kept unless they
#' would collide, because assigning a duplicated label to `levels<-`
#' MERGES the levels and would undo the fix.
#'
#' @noRd
ce_effect2 <- function(v, v2 = NULL) {
  if (!is.numeric(v) || is.matrix(v)) {
    return(if (is.factor(v)) v else factor(v))
  }
  uv <- sort(unique(v), decreasing = TRUE)
  out <- factor(match(v, uv), levels = seq_along(uv))
  nm <- names(v2)
  labs <- if (!is.null(nm) && length(nm) == length(uv) && all(nzchar(nm))) {
    nm[order(v2, decreasing = TRUE)]
  } else {
    d <- 2L
    while (d < 15L && anyDuplicated(round(uv, d))) d <- d + 1L
    as.character(round(uv, d))
  }
  if (anyDuplicated(labs)) labs <- make.unique(labs)
  levels(out) <- labs
  out
}

#' The scale a Wald band for the EXPECTED RESPONSE is symmetric on.
#'
#' It has to reduce exactly to the link-scale band the ordinary branch
#' draws, since that is what every family whose mean is the inverse link
#' of mu has always got: with `m = linkinv(eta)` and
#' `se_m = |mu_eta| se_eta`, `linkinv(link(m) +/- z se_m / |mu_eta|)` IS
#' `linkinv(eta +/- z se_eta)`. So the mu link is tried first and kept
#' only when it can hold the mean - a binomial mean is a COUNT under
#' `trials()`, which the logit link's domain does not contain. The
#' fallbacks are the log scale, which cannot cross zero, and the
#' response scale itself for a mean that can be negative (a truncated
#' gaussian).
#'
#' @noRd
ce_band_scale <- function(link, m) {
  fin <- is.finite(m)
  if (any(fin)) {
    t <- suppressWarnings(link$linkfun(m))
    back <- suppressWarnings(link$linkinv(t))
    d <- suppressWarnings(link$mu_eta(t))
    ok <- all(is.finite(t[fin])) && all(is.finite(d[fin])) &&
      all(abs(d[fin]) > 0) &&
      max(abs(back[fin] - m[fin])) <= 1e-8 * max(1, max(abs(m[fin])))
    if (isTRUE(ok)) return(list(t = t, inv = link$linkinv, dm = d))
    if (all(m[fin] > 0)) {
      return(list(t = log(m), inv = exp, dm = m))
    }
  }
  list(t = m, inv = identity, dm = rep(1, length(m)))
}

#' brms's `method =` vocabulary, accepted alongside frmtmb's own.
#'
#' The values agree wherever both spellings resolve, so a ported call
#' failing at `match.arg()` was a stumble and nothing else.
#'
#' @noRd
ce_method <- function(method) {
  if (length(method) == 1L && is.character(method)) {
    method <- switch(method,
      posterior_epred = "epred",
      posterior_predict = "predict",
      posterior_linpred = stop(
        "conditional_effects(method = \"posterior_linpred\") has no ",
        "frmtmb spelling: ask for the linear predictor's own display ",
        "with dpar = instead", call. = FALSE),
      method)
  }
  match.arg(method, c("epred", "predict"))
}

#' The dots this function itself accepts, with the rest reported against
#' THIS function rather than against the `predict()` it used to forward
#' them to - a function the user did not call, and only on the branches
#' that forward, so an ordinal display swallowed them in silence.
#'
#' @noRd
ce_dots <- function(dots) {
  anl <- FALSE
  for (nm in c("allow_new_levels", "allow.new.levels")) {
    if (nm %in% names(dots)) {
      anl <- isTRUE(dots[[nm]])
      dots[[nm]] <- NULL
    }
  }
  if (length(dots)) {
    warning("conditional_effects() is ignoring unknown argument(s): ",
            paste(names(dots), collapse = ", "), call. = FALSE)
  }
  anl
}

#' The grouping variables a curve conditions on. Smooth, gp() and hsgp()
#' blocks are not groups: their "levels" are basis functions, and
#' blanking their variable would remove the curve rather than the group.
#'
#' @noRd
ce_group_vars <- function(x) {
  bks <- Filter(function(bk) {
    !bk[["covstruct"]] %in% c("smooth", "gp", "hsgp")
  }, x$frame[["re_blocks"]] %||% list())
  g <- vapply(bks, function(bk) bk[["group_name"]] %||% "", "")
  g <- unique(unlist(strsplit(g, ":", fixed = TRUE)))
  g[nzchar(g)]
}

#' Which display one call asks for: `"cats"` the per-category
#' probabilities (brms's `categorical = TRUE`), `"cats_mean"` the
#' expected category number (brms's ordinal default), `"linpred"` one
#' linear predictor's own curve.
#'
#' frmtmb keeps the per-category display as its DEFAULT for a polytomous
#' family, which is the layout brms's own message asks the user to
#' switch to. What was a defect is that `categorical =` was accepted and
#' did nothing, so the other layout could not be asked for at all.
#'
#' @noRd
ce_display_kind <- function(rspec, dpar, categorical) {
  poly <- isTRUE(rspec$family[["type"]] %in% c("ordinal", "categorical"))
  if (!is.null(categorical)) check_flag(categorical, "categorical")
  if (!is.null(dpar) || !poly) {
    if (isTRUE(categorical)) {
      stop("conditional_effects(categorical = TRUE) needs an ordinal or ",
           "categorical family and no dpar = : it draws one curve per ",
           "response category, and ",
           if (!is.null(dpar)) {
             paste0("dpar = \"", dpar, "\" asks for that linear ",
                    "predictor instead")
           } else {
             paste0("family '", rspec$family[["family"]],
                    "' has no response categories")
           }, call. = FALSE)
    }
    return("linpred")
  }
  if (is.null(categorical) || isTRUE(categorical)) return("cats")
  if (!identical(rspec$family[["type"]], "ordinal")) {
    stop("conditional_effects(categorical = FALSE) has nothing to draw ",
         "for family '", rspec$family[["family"]], "': the expected ",
         "category number needs ORDERED categories, and a nominal ",
         "family's are not ordered. Use categorical = TRUE (the ",
         "default here), or dpar = for one category's predictor",
         call. = FALSE)
  }
  "cats_mean"
}

#' Whether the display is per response CATEGORY rather than one curve.
#'
#' The contract, not the family name: a family whose `type` says the
#' modelled response is a set of categories predicts an `n x K`
#' probability matrix on the response scale, and that matrix is what the
#' display shows - one curve per column, brms's `categorical = TRUE`.
#' `"ordinal"` covers `cumulative()`, `sratio()`, `cratio()`, `acat()`;
#' `"categorical"` is the nominal family, which arrives at merge and
#' takes the same display through the same predict() call. Naming a
#' `dpar` opts back into the ordinary linear-predictor display.
#'
#' @noRd
ce_cats_display <- function(rspec, dpar) {
  is.null(dpar) &&
    isTRUE(rspec$family[["type"]] %in% c("ordinal", "categorical"))
}

#' The dpar `conditional_effects()` actually PREDICTS, which is not
#' always the dpar it labels the display with.
#'
#' `NULL` means the expected response, which is what
#' `predict(type = "response")` returns with no `dpar =` and what
#' `method = "epred"` has always been documented to draw. It is the
#' answer whenever no dpar was named and the family's mean is not the
#' inverse link of `mu` - a zero-inflated or hurdle family, a mixture -
#' or the response carries `trunc()` bounds, because in all of those the
#' mu predictor alone is a different quantity from the mean. Taking the
#' estimate from mu regardless plotted `(1 - zi)` times too little on a
#' zero-inflated fit, one component's mean on a mixture, and a NEGATIVE
#' mean on a response truncated below at zero.
#'
#' Exported rather than left to each caller to re-derive: `mean_is_mu()`
#' is a structural test of a family's `mean_fn` body and has no business
#' crossing a package boundary, but the DECISION it feeds does. A
#' sampling extension that draws the same curves per posterior draw has
#' to make the same choice, and there is exactly one right answer per
#' model.
#'
#' `dpar` is the resolved (defaulted) label, `dpar_given` whether the
#' user named one; the two category displays never take this path
#' because their quantity is a probability, not a dpar.
#'
#' @noRd
ce_pred_dpar <- function(rspec, dpar, dpar_given = FALSE,
                         categorical = FALSE, cats_mean = FALSE) {
  mean_display <- !categorical && !cats_mean && !dpar_given &&
    (!mean_is_mu(rspec$family) || has_trunc(rspec))
  if (mean_display) NULL else dpar
}

#' One conditional-effects grid: every predictor at its reference value
#' or at its `conditions` override, with the varied predictor(s) replaced
#' by their grid values. Split out because `band = "boot"` needs every
#' grid of the call BEFORE any refit happens, so that one bootstrap can
#' serve all of them.
#'
#' @noRd
ce_build_nd <- function(base, ev, v1, v2, cset, n, n2,
                        na_vars = character(0)) {
  nd <- data.frame(.ce_row = seq_len(n))
  for (nm in names(base)) {
    val <- if (nm %in% names(cset)) {
      cnd <- cset[[nm]]
      if (is.factor(base[[nm]])) {
        factor(cnd, levels = levels(base[[nm]]))
      } else {
        cnd
      }
    } else if (nm %in% na_vars) {
      # a NEW group: the level is unobserved, so the column says so
      # rather than naming an arbitrary observed one
      if (is.factor(base[[nm]])) {
        factor(NA, levels = levels(base[[nm]]))
      } else {
        base[[nm]][NA_integer_][1L]
      }
    } else {
      ce_ref_value(base[[nm]])
    }
    nd[[nm]] <- if (is.matrix(val)) {
      matrix(val, n, ncol(val), byrow = TRUE)
    } else {
      rep(val, length.out = n)
    }
  }
  # brms's row order: the FIRST effect varies slowest, so the moderator
  # moves within a block of one x value. Both orders hold the same
  # points; matching brms means a script that indexes rows positionally
  # ports, and it costs nothing.
  nd[[ev[1L]]] <- rep(v1, each = n2)
  if (length(ev) == 2L) nd[[ev[2L]]] <- rep(v2, times = length(v1))
  nd$.ce_row <- NULL
  nd
}

#' Population-level display values of one grid under one (re)fit: the
#' same numbers `estimate__` carries, flattened. The ordinal display is
#' per category, so the K-column probability matrix flattens
#' column-major, which is the category-major row order the ordinal data
#' frame is built in.
#'
#' @noRd
ce_boot_one <- function(fit, nd, categorical, resp, dpar,
                        re_form = NA, allow_new_levels = FALSE) {
  p <- if (categorical) {
    predict(fit, newdata = nd, type = "response", resp = resp,
            re.form = re_form, allow_new_levels = allow_new_levels)
  } else {
    predict(fit, newdata = nd, type = "response", dpar = dpar,
            resp = resp, re.form = re_form,
            allow_new_levels = allow_new_levels)
  }
  as.vector(p)
}

#' The population switch, in brms's spelling. conditional_effects() IS
#' the brms function, so it takes brms's `re_formula`; frmtmb's fit
#' surface spells the same setting `re.form` after lme4, and a user
#' who reaches for that spelling here is told which one this function
#' takes instead of hitting a matched-by-multiple-arguments error from
#' the internal predict() calls.
#'
#' @noRd
ce_re_formula <- function(re_formula, dots) {
  if ("re.form" %in% names(dots)) {
    stop("conditional_effects() spells this argument `re_formula` ",
         "(brms's spelling; the fit surface's predict() and ",
         "simulate() spell it `re.form` after lme4). Pass ",
         "re_formula = ", call. = FALSE)
  }
  if (!is.null(re_formula) && !inherits(re_formula, "formula") &&
        !(length(re_formula) == 1L && is.na(re_formula))) {
    stop("`re_formula` must be NA to draw the population-level curve ",
         "(the default), NULL to condition on the grid reference ",
         "group levels, or a one-sided formula naming the terms to ",
         "keep, not ", arg_desc(re_formula), call. = FALSE)
  }
  re_formula
}

#' Identity of the grids one bootstrap was run over.
#'
#' The reuse check has to compare the grid CONTENT, not its shape. A
#' bootstrap taken under `conditions = list(f = "a")` holds percentiles
#' for a different curve than the one `conditions = list(f = "b")`
#' draws, yet the two calls agree on every structural feature - same
#' effect names, same row counts, same column layout - so a key built
#' from those accepted the wrong draws silently and produced a band that
#' need not contain its own estimate. The grid values and the condition
#' sets therefore go into the key themselves.
#'
#' `serialize()` rather than a hash: no new dependency, and an exact
#' comparison has no collision to reason about. The data frames are
#' stripped to their columns first, so the key does not turn on row
#' names or on the class attribute. `prob` is deliberately absent - the
#' same draws answer any coverage.
#'
#' @noRd
ce_boot_key <- function(grids, categorical, resp, dpar, lens,
                        nspec = list()) {
  serialize(list(
    resp = resp, dpar = dpar, categorical = categorical, lens = lens,
    new_level = ce_new_level_key(nspec),
    grids = lapply(grids, function(g) {
      list(eff = g$eff, ci = g$ci, cset = g$cset, nd = as.list(g$nd))
    })
  ), NULL, xdr = FALSE)
}

#' WHICH GROUP a bootstrap's draws belong to, as a comparable value.
#'
#' The grid alone cannot say. `ce_ref_value()` holds an unvaried factor
#' at `levels(col)[1L]` and `ce_new_level_spec()` uses
#' `bk[["levels"]][1L]` as its placeholder, so on sleepstudy both are
#' `"308"` and a population grid and a new-group grid are BYTE
#' IDENTICAL. Without this in the key, `boot =` reuse accepted a
#' population bootstrap for a `re_formula = NULL` call and handed back
#' the population band, which is the exact defect the per-replicate
#' draw was added to fix, reached through the path the help page
#' recommends; the converse handed a population call a band four times
#' too wide.
#'
#' The block object itself is deliberately not in here: what identifies
#' the draws is which blocks are redrawn, where they are written, and
#' whether they are drawn or zeroed.
#'
#' @noRd
ce_new_level_key <- function(nspec) {
  lapply(nspec, function(sp) {
    list(vars = sp$vars, parts = sp$parts, idx = sp$idx,
         dim = sp$dim, rr = sp$rr, draw = sp$draw)
  })
}

#' What a NEW group's effects are, per random-effect block, so that
#' `band = "boot"` under `re_formula = NULL` means what the wald band
#' means.
#'
#' The wald band adds an unseen level's MARGINAL variance to the linear
#' predictor (`lp_extra_var()`). A bootstrap has no variance to add: it
#' has replicates, so the analog is to DRAW that level's effects once
#' per replicate. Without the draw every refit predicted the new
#' level's ZERO modes and the interval was the POPULATION interval
#' under another name, bit-identically so on an ordinal fit, where
#' `band = "boot"` is the only band this function allows with
#' `re_formula`.
#'
#' The mechanism is a placeholder level: the boot grid carries an
#' OBSERVED level so the design maps it, and each replicate overwrites
#' that level's coefficients, which is exactly `z_i' u` through the
#' ordinary `Z`. The covariance drawn from, and the blocks that get no
#' draw at all, are `lp_extra_var()`'s own choices rather than a second
#' opinion: a block whose levels ARE the structure (`gr_cov`,
#' `gr_prec`, `car`, `spde`) has no marginal covariance for an unseen
#' level, so its placeholder entries are ZEROED, which is what the
#' wald band assumes for it too.
#'
#' @noRd
ce_new_level_spec <- function(fit, na_vars, base) {
  out <- list()
  for (bk in fit$frame[["re_blocks"]] %||% list()) {
    if (bk[["covstruct"]] %in% c("smooth", "gp", "hsgp")) next
    gv <- unlist(strsplit(bk[["group_name"]] %||% "", ":", fixed = TRUE))
    gv <- gv[nzchar(gv)]
    if (!length(gv) || !all(gv %in% na_vars)) next
    lvl <- bk[["levels"]][1L]
    parts <- strsplit(lvl, ":", fixed = TRUE)[[1L]]
    if (length(parts) != length(gv)) next
    if (!all(vapply(gv, function(v) !is.null(base[[v]]), TRUE))) next
    # rr blocks carry standard-normal FACTORS, rank per level; every
    # other block carries the coefficients themselves, dim per level
    rr <- identical(bk[["covstruct"]], "rr")
    d <- if (rr) bk[["rank"]] else bk[["dim"]]
    out[[length(out) + 1L]] <- list(
      bk = bk, rr = rr, dim = d,
      idx = bk[["b_idx"]][seq_len(d)],   # level one occupies the first d
      vars = gv, parts = parts,
      draw = rr || !bk[["covstruct"]] %in%
        c("gr_cov", "gr_prec", "car", "spde")
    )
  }
  out
}

#' The boot grids: the same points, with the placeholder level in place
#' of the `NA` the returned frame reports.
#'
#' @noRd
ce_boot_grids <- function(grids, nspec, base) {
  for (gi in seq_along(grids)) {
    for (sp in nspec) {
      for (k in seq_along(sp$vars)) {
        v <- sp$vars[[k]]
        col <- base[[v]]
        val <- if (is.factor(col)) {
          factor(sp$parts[[k]], levels = levels(col))
        } else if (is.numeric(col)) {
          as.numeric(sp$parts[[k]])
        } else {
          sp$parts[[k]]
        }
        grids[[gi]]$nd[[v]] <- rep(val, length.out = nrow(grids[[gi]]$nd))
      }
    }
  }
  grids
}

#' One replicate's new-group draw, written into the placeholder level.
#'
#' @noRd
ce_draw_new_levels <- function(f, nspec) {
  b <- f$estimates[["b"]]
  th <- f$estimates[["theta"]]
  # a refit that came back without the blocks this spec was built from
  # gets no draw rather than a length error inside the replicate
  need <- max(unlist(lapply(nspec, function(sp) sp$idx)), 0L)
  if (!length(b) || length(b) < need) return(f)
  for (sp in nspec) {
    if (!isTRUE(sp$draw)) {
      b[sp$idx] <- 0
      next
    }
    if (sp$rr) {
      b[sp$idx] <- stats::rnorm(sp$dim)
      next
    }
    bk <- sp$bk
    S <- covstruct_registry[[bk[["covstruct"]]]]$vcov(th[bk[["theta_idx"]]],
                                                      bk)
    if (is_student_block(bk)) S <- S * student_var_factor(bk[["dist_nu"]])
    S <- as.matrix(S)
    L <- tryCatch(t(chol(S)), error = function(e) {
      ev <- eigen(S, symmetric = TRUE)
      ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow(S))
    })
    b[sp$idx] <- as.numeric(L %*% stats::rnorm(nrow(S)))
  }
  f$estimates[["b"]] <- b
  f
}

#' ONE parametric bootstrap for every grid of the call.
#'
#' The draws are grid predictions, not coefficients: `frm_bootstrap()`
#' already takes an arbitrary `FUN` of the refit, so the per-draw work
#' rides on the documented mechanism and nothing is stored for the calls
#' that do not ask for it. It also keeps the refit's own conditional
#' modes in play, which a saved coefficient vector could not do, and it
#' covers the ordinal per-category display for free.
#'
#' A refit whose predictions fail (a degenerate draw, a level that
#' disappeared) contributes an NA row rather than aborting the call:
#' `frm_bootstrap()` guards the refit but not `FUN`.
#'
#' @noRd
ce_boot_draws <- function(x, grids, categorical, resp, dpar, boot,
                          seed, re_form = NA, anl = FALSE,
                          nspec = list()) {
  lens <- vapply(grids, function(g) {
    length(ce_boot_one(x, g$nd, categorical, resp, dpar, re_form, anl))
  }, 1L)
  tot <- sum(lens)
  nkey <- ce_new_level_key(nspec)
  key <- ce_boot_key(grids, categorical, resp, dpar, lens, nspec)
  if (inherits(boot, "frmtmb_boot")) {
    if (!identical(boot$ce_key, key)) {
      stop("boot = was not produced by a conditional_effects(band = ",
           "\"boot\") call on this grid: its draws are ",
           if (is.null(boot$ce_key)) {
             "coefficients or another quantity"
           } else if (!identical(boot$ce_new %||% list(), nkey)) {
             # the grids can be byte identical here, so this reason has
             # to be checked before the grid one or it would never be
             # the one reported
             if (length(nkey)) {
               paste0("predictions for a different group: this call ",
                      "conditions on a NEW group (re_formula = NULL), ",
                      "and those draws do not carry that group's ",
                      "effects, so their percentiles are the ",
                      "population band")
             } else {
               paste0("predictions for a different group: those draws ",
                      "carry a NEW group's effects (they came from a ",
                      "re_formula = NULL call) and this call is the ",
                      "population curve, so their percentiles are too ",
                      "wide for it")
             }
           } else {
             paste0("predictions over a different grid (the effects, ",
                    "the resolution, the conditions or the data are ",
                    "not the same)")
           },
           ", so its percentiles would not be a band for this curve. ",
           "Pass a number of draws instead, or reuse attr(ce, \"boot\") ",
           "from an otherwise identical call", call. = FALSE)
    }
    return(list(bs = boot, lens = lens,
                offsets = cumsum(c(0L, lens))))
  }
  if (!is.null(boot) &&
      !(is.numeric(boot) && length(boot) == 1L && is.finite(boot) &&
        boot >= 2)) {
    stop("boot = takes a single number of bootstrap draws (>= 2) or a ",
         "frmtmb_boot object; got ", class(boot)[1L], call. = FALSE)
  }
  nsim <- if (is.null(boot)) 200L else as.integer(boot)
  FUN <- function(f) {
    # the new group is drawn ONCE per replicate, before the grids are
    # evaluated, so every panel of the call shares that group
    if (length(nspec)) f <- ce_draw_new_levels(f, nspec)
    v <- tryCatch(
      unlist(lapply(grids, function(g) {
        ce_boot_one(f, g$nd, categorical, resp, dpar, re_form, anl)
      }), use.names = FALSE),
      error = function(e) NULL
    )
    if (!is.numeric(v) || length(v) != tot) rep(NA_real_, tot) else v
  }
  if (is.null(boot)) {
    message("conditional_effects(band = \"boot\"): refitting the model ",
            nsim, " times (one bootstrap shared by all ",
            length(grids), " grid(s)). Pass boot = <draws> for a ",
            "cheaper run, or boot = attr(ce, \"boot\") to reuse this one.")
  }
  bs <- frm_bootstrap(x, FUN = FUN, nsim = nsim, seed = seed)
  bs$ce_key <- key
  # kept beside the key so a mismatch can say WHICH of the two it is
  bs$ce_new <- nkey
  list(bs = bs, lens = lens, offsets = cumsum(c(0L, lens)))
}

#' Component label per position of the outer (optimized) parameter
#' vector, in `obj$par` order. `outer_par_names()` gives the names; the
#' profile band needs to know which block each position belongs to so a
#' design row can be written as a `lincomb` over `obj$par`.
#'
#' @noRd
ce_outer_comp <- function(fit) {
  tpl <- fit$frame[["par_template"]]
  random <- c("b", "miss")
  if (fit$REML || isTRUE(fit$control$profile)) random <- c(random, "beta")
  comp <- character(0)
  for (cp in names(tpl)) {
    if (cp %in% random) next
    len <- length(tpl[[cp]])
    if (cp == "betad" && length(fit$frame[["betad_fixed_idx"]])) {
      len <- len - length(fit$frame[["betad_fixed_idx"]])
    }
    comp <- c(comp, rep(cp, len))
  }
  comp
}

#' Refusals for `band = "profile"`, all raised before any grid is built.
#'
#' The likelihood-root search inverts the LR for ONE linear combination
#' of the outer parameters. That covers the linear predictor of a single
#' distributional parameter, and (the link being monotone) anything the
#' link maps it to. It does not cover a quantity that mixes several
#' predictors, that runs through inner parameters, or that is not a
#' linear combination at all - so those are refused rather than
#' approximated.
#'
#' @noRd
ce_profile_check <- function(x, rspec, lp, dpar_given, categorical) {
  if (categorical) {
    stop("band = \"profile\" cannot cover an ordinal category ",
         "probability: it is not a linear combination of the ",
         "parameters (the thresholds and any cs() coefficients enter ",
         "every category), so there is no single likelihood root to ",
         "invert. Use band = \"boot\", or dpar = \"mu\" for the latent ",
         "predictor", call. = FALSE)
  }
  if (x$REML) {
    stop("band = \"profile\" requires an ML fit: REML integrates the ",
         "fixed effects out of the outer problem, so the effect grid is ",
         "not a function of the parameters the likelihood is profiled ",
         "over. Use band = \"boot\", or refit with REML = FALSE",
         call. = FALSE)
  }
  if (isTRUE(x$control$profile)) {
    stop("band = \"profile\" needs a fit without ",
         "frmtmb_control(profile = TRUE): the profiled coefficients are ",
         "not outer parameters there. Use band = \"boot\"",
         call. = FALSE)
  }
  if (!is.null(lp[["nl_body"]])) {
    stop("band = \"profile\" is not available for a nonlinear ",
         "predictor: a grid value is not a linear combination of the ",
         "nonlinear parameters. Use band = \"boot\"", call. = FALSE)
  }
  if (!dpar_given && (!mean_is_mu(rspec$family) || has_trunc(rspec))) {
    stop("band = \"profile\" cannot cover the expected response of ",
         "family '", rspec$family[["family"]], "': it runs through more ",
         "than one distributional parameter (zero inflation, a hurdle, ",
         "a dispersion) or through truncation bounds, so it is not a ",
         "single linear combination of the parameters. Use ",
         "band = \"boot\", or name one predictor with dpar =",
         call. = FALSE)
  }
  key <- linpred_key(lp[["resp"]], lp[["dpar"]])
  sm <- Filter(function(bk) {
    bk[["covstruct"]] %in% c("smooth", "gp", "hsgp") &&
      any(vapply(bk[["components"]], function(cp) cp$lp_key == key, TRUE))
  }, x$frame[["re_blocks"]])
  if (length(sm)) {
    stop("band = \"profile\" cannot cover a predictor carrying a ",
         "smooth, gp() or hsgp() term: the basis coefficients are ",
         "inner (random) parameters, which a likelihood-root search ",
         "over the outer parameter vector does not move. Use ",
         "band = \"boot\"", call. = FALSE)
  }
  invisible(NULL)
}

#' Profile-likelihood band for one grid, on the LINK scale.
#'
#' Row `i` of the grid has `eta_i = a_i' par + c_i` with `a_i` the design
#' row over this predictor's coefficients and `c_i` whatever does not
#' move with them (the offset). [TMB::tmbroot()] inverts the likelihood
#' ratio along `lincomb = a_i`, and `c_i` shifts the interval back. The
#' link is monotone, so the caller maps the two endpoints through it.
#'
#' A numeric grid is profiled at `profile_points` points and the
#' endpoints are interpolated linearly between them: a root search is
#' two constrained optimizations, and 100 of them per effect is not a
#' default anyone would wait for. Interpolation happens on the link
#' scale, where the endpoints are smooth in the predictor.
#'
#' @noRd
ce_profile_eta_ci <- function(x, lp, nd, v1, n1, n2, prob,
                              profile_points) {
  ed <- lp_eta_design(x, lp, nd, FALSE, FALSE)
  comp <- ce_outer_comp(x)
  par <- x$opt$par
  X <- as.matrix(ed[["X"]])
  if (lp[["par"]] == "beta") {
    pos <- which(comp == "beta")[lp[["idx"]]]
    keep <- rep(TRUE, ncol(X))
  } else {
    tpl_len <- length(x$frame[["par_template"]][["betad"]])
    est_rank <- match(lp[["idx"]], setdiff(seq_len(tpl_len),
                                      x$frame[["betad_fixed_idx"]]))
    keep <- !is.na(est_rank)
    pos <- which(comp == "betad")[est_rank[keep]]
  }
  X <- X[, keep, drop = FALSE]
  if (length(comp) != length(par) || anyNA(pos) || length(pos) != ncol(X)) {
    stop("band = \"profile\" cannot line this fit's coefficients up ",
         "with its outer parameter vector (a mapped, fixed or profiled ",
         "coefficient block). Use band = \"boot\"", call. = FALSE)
  }
  eta <- ed[["eta"]]
  n <- length(eta)
  target <- 0.5 * stats::qchisq(prob, df = 1)
  sel <- if (is.numeric(v1) && !is.matrix(v1) && n1 > profile_points) {
    unique(round(seq(1, n1, length.out = max(2L, profile_points))))
  } else {
    seq_len(n1)
  }
  # the grid varies the first effect slowest (ce_build_nd), so the rows
  # of grid point `i` are (i - 1) * n2 + 1 ... i * n2
  rows <- as.vector(outer((sel - 1L) * n2, seq_len(n2), "+"))
  lo <- up <- rep(NA_real_, n)
  fails <- 0L
  for (r in rows) {
    if (!is.finite(eta[r])) next
    a <- X[r, ]
    if (all(a == 0)) {
      # nothing this fit estimates moves the row: the offset IS the
      # linear predictor, and it carries no uncertainty
      lo[r] <- up[r] <- eta[r]
      next
    }
    v <- numeric(length(par))
    v[pos] <- a
    const <- eta[r] - sum(a * par[pos])
    ci <- tryCatch(
      suppressWarnings(TMB::tmbroot(x$obj, lincomb = v, target = target)),
      error = function(e) c(NA_real_, NA_real_)
    )
    if (length(ci) != 2L || !all(is.finite(ci))) {
      fails <- fails + 1L
      next
    }
    lo[r] <- min(ci) + const
    up[r] <- max(ci) + const
  }
  if (length(sel) < n1) {
    xs <- as.numeric(v1)
    for (j in seq_len(n2)) {
      idx <- (seq_len(n1) - 1L) * n2 + j
      fill <- function(y) {
        if (sum(is.finite(y)) < 2L) return(rep(NA_real_, n1))
        # na.rm = FALSE: an interval touching a failed point stays NA
        # rather than being bridged over silently
        stats::approx(xs[sel], y, xout = xs, na.rm = FALSE)$y
      }
      lo[idx] <- fill(lo[idx[sel]])
      up[idx] <- fill(up[idx[sel]])
    }
  }
  list(lower = lo, upper = up, fails = fails, tried = length(rows))
}

#' Conditional effects of predictors
#'
#' For each requested effect, predicts over a grid of that predictor
#' with every other predictor held at a reference value (numeric: mean;
#' factor: first level; matrix covariate: column means) and random
#' effects excluded (`re.form = NA`). Confidence bands are Wald
#' intervals computed on the link scale and back-transformed. Smooth
#' terms are included, so this also covers what brms calls
#' `conditional_smooths()`.
#'
#' @param x A `frmtmb_fit`.
#' @param effects Character vector of variable names, or `"x:z"` pairs;
#'   for a pair, the first variable is varied over its range while the
#'   second is held at its levels (factors) or at mean and mean plus or
#'   minus one SD (numeric). Default: every fixed-effect and smooth
#'   variable of the selected linear predictor, plus one `"a:b"` pair
#'   per fitted interaction (brms's default); a term of order three or
#'   more contributes its leading pair.
#' @param resp,dpar Response and distributional parameter, as in
#'   [predict.frmtmb_fit()].
#' @param resolution Number of grid points for a varied numeric
#'   predictor.
#' @param prob Coverage of the confidence bands (brms spelling).
#' @param re_formula The population switch, in brms's spelling: `NA`
#'   (the default) draws the population-level curve, `NULL` conditions
#'   on a NEW, unobserved group, and a one-sided formula keeps the named
#'   terms for one. A new group's conditional modes are zero, so its
#'   curve IS the population curve and what the group costs is spread:
#'   the band carries the random-effect variance on top of the
#'   coefficient uncertainty, and the grouping column of the returned
#'   frame is `NA` to say which group it is. `band = "boot"` carries it
#'   too, by drawing that group's effects once per bootstrap replicate
#'   rather than by adding a variance. To condition on an OBSERVED
#'   group, name it in `conditions` (`conditions = list(g = "3")`).
#'   brms draws a new group's random effects afresh from the fitted
#'   covariance in every posterior draw, so its curve is stochastic
#'   around this one; a maximum-likelihood fit has the mode and the
#'   variance instead of draws. The fit surface's
#'   [predict.frmtmb_fit()] spells the same
#'   setting `re.form` after lme4; `conditional_effects()` takes brms's
#'   name because it is brms's function, and says so if handed the
#'   other spelling. `band = "profile"` exists only for the
#'   population-level curve.
#' @param band How the confidence band is built: `"wald"` (default,
#'   the delta method on the scale the band is symmetric on),
#'   `"profile"` (likelihood-root
#'   inversion per grid point) or `"boot"` (parametric-bootstrap
#'   percentiles). See the band section. Only for `method = "epred"`.
#' @param boot For `band = "boot"`: `NULL` (default) runs one
#'   [frm_bootstrap()] of 200 refits and says so, a number runs that
#'   many, and a `frmtmb_boot` object from an earlier identical call
#'   (`attr(ce, "boot")`) is reused without refitting anything.
#' @param profile_points For `band = "profile"`: how many points of a
#'   numeric grid are profiled, the band being interpolated between them
#'   on the link scale.
#' @param seed Seed for `band = "boot"`, passed to [frm_bootstrap()].
#' @param method `"epred"` (default): Wald bands for the expected
#'   response, which for a family whose mean is not the inverse link of
#'   its `mu` predictor (zero-inflated, hurdle, `trials()`, truncated)
#'   is the MEAN and not that predictor. brms's own spellings
#'   `"posterior_epred"` and `"posterior_predict"` are accepted as
#'   aliases. `"predict"`: prediction intervals - quantile bands from
#'   `ndraws` responses simulated from the family at each grid point
#'   (observation noise; random effects stay excluded, as in brms with
#'   `re_formula = NA`), around the expected response on the same
#'   scale as the draws (a count under `trials()`, the truncated mean
#'   under `trunc()`). The
#'   draws respect the response's addition terms: literal `trunc()`
#'   bounds apply, and `trials()`, `se()` or variable `trunc()` bounds
#'   must be pinned in `conditions` (a grid row is an artificial
#'   observation, so a reference value for those is meaningless and is
#'   an error rather than a silent default).
#' @param ndraws Simulated responses per grid point for
#'   `method = "predict"`. For the draws method: how many evenly spaced
#'   posterior draws the curves are computed over (default all).
#' @param conditions Named list overriding reference values, e.g.
#'   `list(x2 = 1, g = "b")`; or a data frame whose rows define
#'   multiple condition sets (brms style), labeled by a `cond__`
#'   column from its row names.
#' @param surface Accepted for brms compatibility. `TRUE` (a fitted
#'   surface over two predictors) is refused: ask for the two-variable
#'   effect `"x1:x2"` instead, which varies the first predictor at three
#'   values of the second.
#' @param data The original model data. Only needed when the model frame
#'   does not store a raw variable (e.g. a variable used only inside
#'   `poly()`).
#' @param int_conditions Named list giving the values one or both
#'   variables of an effect are evaluated at, in place of the defaults
#'   (the range of a numeric grid, `mean +/- sd` for a numeric
#'   moderator). An element is either the values themselves or a
#'   function of the observed column, e.g.
#'   `int_conditions = list(z = c(-1, 0, 1))` or
#'   `list(z = function(v) quantile(v, c(0.1, 0.9)))`. Names on a
#'   numeric vector become the moderator's `effect2__` labels. brms's
#'   argument, with brms's meaning.
#' @param categorical Per-category display for a polytomous family:
#'   `TRUE` (the default there) draws one curve per response category
#'   and keys the effect `"x:cats__"`, as brms's `categorical = TRUE`
#'   does; `FALSE` draws the expected CATEGORY NUMBER,
#'   `sum(k * p_k)`, which is brms's default. Only for an ordinal or
#'   categorical family with no `dpar`; a nominal family has no ordered
#'   categories to average and refuses `FALSE`.
#' @param ... `allow_new_levels`, passed to [predict.frmtmb_fit()].
#'   Anything else is reported as unknown, by name, against
#'   `conditional_effects()`.
#' @return A named list of data frames (one per effect), in brms's
#'   column layout: the varied variable(s), then every other model
#'   variable at the value it is held at, then `cond__` (the condition
#'   label, always present), `effect1__` and, for a two-variable effect,
#'   `effect2__` (the moderator as a display label: rounded to two
#'   decimals, levels descending), then `estimate__`, `se__` (on the
#'   scale the band is symmetric on), `lower__` and `upper__`. Printing
#'   it draws the plots. A polytomous fit adds a `cats__` column, one
#'   block of rows per response category, and keys the effect
#'   `"x:cats__"`.
#'   `plot(ce, points = TRUE)` overlays the raw observations (the brms
#'   argument), each panel showing only the observations that belong to
#'   its own condition; see the faceting section. No points are drawn
#'   for a per-category ordinal display, a non-mean `dpar`, or a matrix
#'   response (a message says so).
#' @section Several conditions become one faceted page:
#' A `conditions` data frame of several rows gives the effect one panel
#' per row, laid out as small multiples on a SINGLE page with a shared
#' scale, the way brms's `facet_wrap("cond__")` does.
#' `plot(ce, ncol = )` sets the number of columns; the default lays the
#' panels out roughly square, as brms's `ncol = NULL` does.
#'
#' The panels are drawn with \pkg{tinyplot} when it is installed, which
#' supplies the shared axes and a single outer legend. Without it the
#' fallback is a grid of ordinary base-graphics panels, still one page
#' and still honoring `ncol`, labeled by condition on the y axis. The
#' per-category ordinal display always takes the fallback grid: its
#' grouping slot already carries the response category.
#'
#' `points = TRUE` draws only the observations belonging to each
#' condition, matching brms's `make_point_frame()`. A condition claims
#' the rows of the data that MATCH it on the variables it sets, so with
#' one condition per level of a factor each observation appears once, in
#' its own panel. Two rules bound that, both as in brms:
#'
#' * A NUMERIC condition variable is dropped from the match, because a
#'   reference value such as `list(x2 = 0.37)` is a point on a
#'   continuum that names no observation; matching on it would empty
#'   every panel. A grouping factor stored as a number is exempt, being
#'   a label rather than a continuum. A condition left with nothing to
#'   match on therefore claims EVERY observation, and its panel differs
#'   from the others only in its curve.
#' * A condition variable that names no column of the model data cannot
#'   select anything either, so the points are left unsplit and every
#'   panel draws all of them.
#'
#' Unlike brms, a variable the `conditions` data frame does not set is
#' not silently pinned to its first level for this purpose, so an
#' unmentioned factor does not drop observations from every panel.
#' @section Ordinal responses:
#' `cumulative()`, `sratio()`, `cratio()` and `acat()` have no mean, so
#' the display is per CATEGORY, as brms's `categorical = TRUE` is: each
#' effect data frame gains a `cats__` factor of the response's own
#' levels and carries the fitted category probability in `estimate__`,
#' with one curve per category in the plot (a second predictor gets a
#' panel of its own). A one-variable effect is keyed `"x:cats__"` there,
#' brms's key for that layout.
#'
#' That is the DEFAULT here and it is brms's non-default: brms's own
#' default summarizes the categories into an expected category number
#' and warns that it is treating an ordered factor as continuous.
#' `categorical = FALSE` asks for that summary anyway, `sum(k * p_k)`
#' with its own delta-method band (the category weights go on the
#' gradient, so the covariances between the category probabilities are
#' kept), keyed `"x"` and directly comparable with brms's default
#' curve.
#'
#' `se__` is then on the probability scale, and the band is a Wald
#' interval on the logit of the probability so it cannot leave `[0, 1]`.
#' The standard errors are the delta method over the joint covariance of
#' the coefficients, the thresholds AND the `cs()` coefficients: a
#' category probability depends on all of them, and holding the
#' thresholds fixed would understate every band.
#'
#' `method = "predict"` is refused there (the category probabilities are
#' already the whole predictive distribution). Naming a distributional
#' parameter, `dpar = "mu"`, opts back into the ordinary display of the
#' latent linear predictor.
#' @section Confidence bands:
#' `band` picks how `lower__` and `upper__` are found. The estimate is
#' the same curve in all three cases; only the band changes.
#'
#' - `"wald"` (default, and free): the delta-method standard error,
#'   back-transformed from the scale the band is symmetric on. For the
#'   ordinary display that scale is the link scale. For the expected
#'   response of a family whose mean runs through several distributional
#'   parameters (zero inflation, a hurdle, `trials()`, truncation) the
#'   standard error is the delta method over EVERY predictor's
#'   coefficients jointly, so the cross-parameter covariances are in the
#'   band; the band's scale is then the `mu` link's if the mean lives on
#'   it, the log scale if the mean is positive (so the band cannot cross
#'   zero), and the response scale otherwise. The two rules agree
#'   exactly wherever the mean IS the inverse link of `mu`.
#' - `"profile"`: one likelihood-root search ([TMB::tmbroot()]) per grid
#'   point. A grid point's linear predictor is a linear combination of
#'   the coefficients, so the search inverts the likelihood ratio along
#'   that combination; the link is monotone, so the two endpoints map
#'   through it. The band is not symmetric and does not assume a
#'   quadratic log-likelihood, which is what makes it worth its cost
#'   near a boundary or at a small sample size.
#' - `"boot"`: percentiles of the grid predictions across the refits of
#'   ONE [frm_bootstrap()]. One bootstrap serves every effect, condition
#'   set and ordinal category of the call; `attr(ce, "boot")` returns
#'   it, and passing it back as `boot =` costs no refits at all. Draws
#'   taken under `re_formula = NULL` carry a new group's effects and
#'   draws taken without it do not, so the two are not
#'   interchangeable and `boot =` refuses the swap by name rather than
#'   returning the wrong band.
#'
#' A Student-t random-effect block enters either band as a GAUSSIAN
#' with the t variance, `nu / (nu - 2)` times the scale matrix: the
#' delta method has no other shape to offer, and the bootstrap draws
#' the same way so that the two bands stay comparable. It is the right
#' variance around a heavier-tailed truth, not the right quantile.
#'
#' What `se__` means follows the band: the Wald standard error (link
#' scale, or the probability scale on the ordinal display) for `"wald"`
#' and `"profile"` - the profile changes the endpoints, not the standard
#' error - and the standard deviation of the bootstrap draws, on the
#' displayed scale, for `"boot"`.
#'
#' The three bands answer the same question and need not give the same
#' interval. Wald and bootstrap agree in the middle of a grid, and the
#' remaining difference there is the bootstrap's own Monte Carlo error:
#' on a zero-inflated fit the widths differ by 13% at 200 refits and 4%
#' at 800. At the ENDS of a grid where the estimate is poorly
#' determined they need not converge at all: on a steep zero-inflated
#' shape the Wald band stays about 28% wider at 3000 refits, because it
#' is symmetric on its own scale and the percentile band is not. The
#' point estimate is the same curve either way.
#'
#' Cost, and how it is capped. A root search is two constrained
#' optimizations, so `band = "profile"` profiles at most
#' `profile_points` (default 25) of a numeric grid, spread over its
#' range and including both ends, and interpolates the endpoints
#' linearly between them on the link scale. `resolution` still governs
#' the estimate curve. A factor grid is profiled at every level. Points
#' whose search does not converge become `NA` and one warning names how
#' many. `band = "boot"` costs its refits once, however many effects are
#' asked for.
#'
#' Refusals. `band` other than `"wald"` needs `method = "epred"`: a
#' prediction interval is already a simulation quantile. `band =
#' "profile"` additionally needs the displayed quantity to BE a linear
#' combination of the parameters, so it is refused (naming `"boot"`) for
#' an ordinal category probability, for an expected response that runs
#' through several distributional parameters or through truncation
#' bounds (`dpar =` names one predictor and opts back in), for a
#' nonlinear predictor, for a predictor carrying `s()`, `gp()` or
#' `hsgp()` (basis coefficients are inner parameters), and for a REML
#' fit or `frmtmb_control(profile = TRUE)` (the coefficients are not
#' outer parameters there).
#'
#' A nonlinear predictor (`nl = TRUE`) has no delta-method standard
#' error at all, so `band = "boot"` is its only band and the other two
#' are refused. The same holds for the per-category display of a nominal
#' family, whose probabilities have no threshold Jacobian to
#' differentiate.
#' @section Draws objects:
#' On a `frmtmb_draws` object from `frmtmb.sample::frm_sample()` the
#' same grids are evaluated once per posterior draw, and `estimate__`,
#' `lower__`,
#' `upper__` and `se__` are the pointwise mean, quantiles and standard
#' deviation of the drawn curves. There is no `band =` or `method =` to
#' choose (the band IS the posterior quantile band), a nonlinear
#' predictor and a nominal per-category display work without a delta
#' method, and the method runs on formula-route draws that have no
#' maximum-likelihood fit behind them. `ndraws` thins the draws
#' evenly for a cheaper curve. Draws from a sampling run with
#' `laplace = TRUE`
#' are refused: the sampled vector no longer aligns with the model's
#' parameter template.
#' @section Which predictors are plotted by default:
#' Every variable of the selected linear predictor that the display can
#' vary: its fixed-effect terms, its smooth terms and its `mo()` terms
#' (whose design columns are placeholders, so the variable is read from
#' the term itself). On a nonlinear fit they live one level down - the
#' covariates the nonlinear formula reads, plus the terms of each
#' nonlinear parameter's own predictor - and all of those are collected
#' too. Matrix-valued columns are excluded, a grid over a matrix
#' covariate not being a curve. Every fitted interaction whose
#' components are single plottable variables adds an `"a:b"` display,
#' as in brms: a fitted interaction hidden by default invites reading
#' the main-effect curves as the whole story. The display takes two
#' variables, so a three-way or deeper term contributes its leading
#' pair, with the remaining variables at their reference values until
#' `conditions =` pins them. Naming `effects =` overrides the search.
#' @examples
#' set.seed(5)
#' dd <- data.frame(x = rnorm(120), f = factor(rep(c("a", "b"), 60)))
#' dd$y <- rnorm(120, 1 + 0.5 * dd$x + (dd$f == "b"), 1)
#' fit <- frm(bf(y ~ x * f), family = gaussian(), data = dd)
#' ce <- conditional_effects(fit, effects = c("x", "x:f"))
#' plot(ce, ask = FALSE)
#' # prediction intervals instead of epred bands
#' ce_p <- conditional_effects(fit, effects = "x", method = "predict")
#' \donttest{
#' # a likelihood-profile band: asymmetric, and no quadratic assumption
#' ce_pr <- conditional_effects(fit, effects = "x", band = "profile",
#'                              resolution = 20, profile_points = 5)
#'
#' # a bootstrap band, reused for a second effect without refitting
#' ce_b <- conditional_effects(fit, effects = "x", band = "boot",
#'                             boot = 25, seed = 1)
#' ce_b2 <- conditional_effects(fit, effects = "x", band = "boot",
#'                              boot = attr(ce_b, "boot"))
#' }
#' @export
conditional_effects <- function(x, ...) UseMethod("conditional_effects")

#' A structured family's refusal, shared by the fit and draws methods so
#' the two speak with one message.
#'
#' The expected response of a family whose likelihood does not factorize
#' over rows can depend on the observed responses of a whole sequence or
#' group, which the synthetic grid this function builds does not carry.
#' A structured family whose per-row mean IS rowwise (a group-level
#' mixture) declares `conditional_effects` and passes straight through.
#'
#' @noRd
ce_structure_check <- function(rspec) {
  structure_gate(fam_structure(rspec$family), "conditional_effects",
                 structure_generic(rspec$family, "conditional_effects()"))
  invisible(NULL)
}

#' Effect grids for one conditional_effects() call: the plottable
#' variables, the deduplicated effects, the condition sets and one grid
#' (with its newdata) per effect x condition set. Shared by the fit and
#' draws methods so the two build identical curves.
#'
#' @noRd
ce_grids_build <- function(x, rspec, lp, effects, resp, dpar, resolution,
                           conditions, data, int_conditions = list(),
                           na_vars = character(0)) {
  base <- data %||% x$frame[["data_frame"]]

  vars <- ce_plot_vars(x, rspec, lp, resp)
  # a model whose SELECTED predictor has no terms still has covariates
  # somewhere (bf(y ~ 1, theta1 ~ x)): naming the one predictor the
  # search looked at was a refusal to draw a model that has something
  # to draw
  if (is.null(effects) && !length(intersect(vars, names(base)))) {
    vars <- ce_plot_vars_any(x, rspec, resp)
  }
  vars <- vars[vars %in% names(base)]
  step_vars <- ce_step_vars(x, rspec, lp, resp)
  # a matrix column is a whole function per row, so there is no single
  # value to hold it at while another predictor varies. Remember which
  # ones were dropped: on a scalar-on-function fit they are the only
  # candidates there were, and the refusal has to say so.
  mat_vars <- vars[vapply(vars, function(v) is.matrix(base[[v]]), TRUE)]
  vars <- setdiff(vars, mat_vars)
  if (is.null(effects)) {
    # brms's default: the main effects AND the fitted two-way
    # interactions, each side of a pair having survived the same
    # plottability filter as the main effects
    prs <- ce_plot_pairs(x, rspec, lp, resp)
    prs <- prs[vapply(strsplit(prs, ":", fixed = TRUE),
                      function(v) all(v %in% vars), TRUE)]
    effects <- c(vars, prs)
    if (!length(effects) && length(mat_vars)) {
      terms_on <- vapply(lp[["smooths"]] %||% list(), function(si) {
        if (any(smooth_pred_vars(si$sm) %in% mat_vars)) si$label else ""
      }, "")
      terms_on <- terms_on[nzchar(terms_on)]
      stop("conditional_effects() has nothing to draw for dpar '", dpar,
           "': the only predictor(s) of that parameter are the matrix ",
           "column(s) ", paste0("`", mat_vars, "`", collapse = ", "),
           if (length(terms_on)) {
             paste0(" (carried by ", paste(terms_on, collapse = ", "), ")")
           },
           ", which the display excludes. A matrix column is a whole ",
           "function per row, so it has neither a one-dimensional axis ",
           "to vary along nor a single value to hold the other ",
           "predictors at. Draw the coefficient function with ",
           "predict(newdata = ) over a grid you build yourself: one row ",
           "per grid point, the matrix column holding the grid, and the ",
           "weight column an indicator of the point",
           call. = FALSE)
    }
    if (!length(effects)) {
      stop("No plottable predictors found for dpar '", dpar, "'",
           call. = FALSE)
    }
  }
  # one grid per effect: a repeated name would otherwise stack the same
  # grid twice inside its own data frame
  effects <- unique(effects)
  # a misspelled int_conditions name silently conditioned on nothing
  unknown <- setdiff(names(int_conditions), names(base))
  if (length(unknown)) {
    warning("int_conditions names no variable of the model data: ",
            paste(unknown, collapse = ", "), call. = FALSE)
  }

  # a data-frame `conditions` defines one condition set per row (brms
  # style); a named list is a single condition set
  cond_sets <- if (is.data.frame(conditions)) {
    stats::setNames(lapply(seq_len(nrow(conditions)), function(r) {
      as.list(conditions[r, , drop = FALSE])
    }), rownames(conditions))
  } else {
    list(conditions)
  }

  tv <- ce_trial_vars(rspec, base)
  if (length(tv)) {
    unpinned <- character(0)
    for (i in seq_along(cond_sets)) {
      miss <- setdiff(tv, names(cond_sets[[i]]))
      unpinned <- union(unpinned, miss)
      for (v in miss) cond_sets[[i]][[v]] <- 1
    }
    if (length(unpinned)) {
      message("conditional_effects(): holding the trials variable(s) ",
              paste(unpinned, collapse = ", "), " at 1, so the display ",
              "is a probability per trial (brms's default too). Pin ",
              "them in conditions = list(...) for a count.")
    }
  }
  grids <- list()
  for (eff in effects) {
    ev <- strsplit(eff, ":", fixed = TRUE)[[1L]]
    if (length(ev) > 2L) {
      stop("Effects support at most two variables: '", eff, "'",
           call. = FALSE)
    }
    missing_ev <- setdiff(ev, names(base))
    if (length(missing_ev)) {
      stop("Variable '", missing_ev[1L], "' is not stored in the model ",
           "frame; pass the original data via data =", call. = FALSE)
    }
    v1 <- ce_grid_values(base[[ev[1L]]], resolution, ev[1L],
                         int_conditions[[ev[1L]]],
                         stepwise = ev[1L] %in% step_vars)
    v2 <- if (length(ev) == 2L) {
      ce_second_values(base[[ev[2L]]], int_conditions[[ev[2L]]])
    }
    n1 <- length(v1)
    n2 <- max(1L, length(v2))
    for (ci in seq_along(cond_sets)) {
      grids[[length(grids) + 1L]] <- list(
        eff = eff, ev = ev, ci = ci, v1 = v1, v2 = v2, n1 = n1, n2 = n2,
        n = n1 * n2, cset = cond_sets[[ci]],
        nd = ce_build_nd(base, ev, v1, v2, cond_sets[[ci]], n1 * n2, n2,
                         na_vars)
      )
    }
  }
  # brms exempts a random-effects grouping factor from the rule that
  # drops numeric condition variables from the raw-point match: a group
  # id is a label an observation carries, even when it is stored as a
  # number, so it matches exactly
  groups <- vapply(x$frame[["re_blocks"]] %||% list(),
                   function(bk) bk[["group_name"]] %||% "", "")
  groups <- unique(unlist(strsplit(groups, ":", fixed = TRUE)))

  list(base = base, effects = effects, cond_sets = cond_sets,
       grids = grids, groups = groups[nzchar(groups)])
}

#' Split the raw observations across the condition sets the way brms's
#' make_point_frame() does, so each panel shows its own condition's data
#' rather than all of it.
#'
#' A condition row claims the observations that MATCH it on the
#' variables it sets. Two rules make that workable, both measured
#' against brms rather than invented here:
#'
#' A numeric condition variable is dropped from the match, because a
#' reference value like `z = 0.37` is a point on a continuum that names
#' no observation; filtering on it would empty every panel. (This is
#' brms's `select_points = 0` default. A grouping factor stored as a
#' number is exempt: it is a label, not a continuum.) A condition that
#' has nothing left to match on therefore claims every observation, and
#' the panels differ only in their curves.
#'
#' A condition variable that names no column of the model data cannot
#' select anything either, so the points are left unsplit and every
#' panel draws all of them, as brms's ggplot layer does when the point
#' frame carries no `cond__`.
#'
#' @noRd
ce_points_by_cond <- function(pts, base, cond_sets, ev, groups) {
  if (is.null(pts) || length(cond_sets) < 2L) return(pts)
  cvars <- intersect(names(cond_sets[[1L]]), setdiff(names(base), ev))
  if (!length(cvars)) return(pts)
  labs <- names(cond_sets) %||% as.character(seq_along(cond_sets))
  parts <- lapply(seq_along(cond_sets), function(i) {
    cs <- cond_sets[[i]][cvars]
    # brms treats NA and its "zero__" placeholder as "unset here"
    set <- vapply(cs, function(v) {
      length(v) == 1L && !is.na(v) && !identical(as.character(v), "zero__")
    }, TRUE)
    mv <- cvars[set]
    keep <- !vapply(mv, function(v) is.numeric(base[[v]]), TRUE) |
      mv %in% groups
    mv <- mv[keep]
    rows <- if (!length(mv)) {
      seq_len(nrow(pts))
    } else {
      # "\r" cannot occur in a factor level, so pasting the columns
      # together compares the whole condition in one shot
      lhs <- do.call(paste, c(lapply(mv, function(v) {
        as.character(base[[v]])
      }), sep = "\r"))
      rhs <- do.call(paste, c(lapply(cs[mv], as.character), sep = "\r"))
      which(lhs %in% rhs)
    }
    if (!length(rows)) return(NULL)
    d <- pts[rows, , drop = FALSE]
    d$cond__ <- labs[i]
    d
  })
  out <- do.call(rbind, parts)
  if (is.null(out)) {
    out <- pts[0L, , drop = FALSE]
    out$cond__ <- character(0)
  }
  out$cond__ <- factor(out$cond__, levels = labs)
  rownames(out) <- NULL
  out
}

#' Assemble the per-effect data frames into the classed result, with the
#' display attributes and the raw points for plot(points = TRUE).
#' Shared by the fit and draws methods.
#'
#' @noRd
ce_finalize <- function(dfs_by_eff, effects, rspec, resp, dpar, band,
                        base, categorical, cond_sets = list(),
                        groups = character(0), cats_key = FALSE) {
  out <- list()
  for (eff in effects) {
    ev <- strsplit(eff, ":", fixed = TRUE)[[1L]]
    df <- do.call(rbind, dfs_by_eff[[eff]])
    # brms keys the per-category layout "x:cats__", the category being
    # the second display dimension, and its own plot() reads the pair
    # back out of the effects attribute
    key <- eff
    if (isTRUE(cats_key) && length(ev) == 1L) {
      ev <- c(ev, "cats__")
      key <- paste(ev, collapse = ":")
    }
    attr(df, "effects") <- ev
    attr(df, "response") <- resp
    attr(df, "dpar") <- dpar
    attr(df, "band") <- band
    # raw observations for plot(..., points = TRUE): only meaningful on
    # the expected-response display, and only when the response is a
    # plain numeric column of the data (not cbind()/matrix responses)
    default_dpar <- if ("mu" %in% names(rspec$dpars)) "mu" else
      rspec$primary_dpars[1]
    if (!categorical && identical(dpar, default_dpar) &&
        resp %in% names(base) && is.numeric(base[[resp]]) &&
        is.null(dim(base[[resp]]))) {
      pdf_ <- data.frame(x = base[[ev[1L]]], y = base[[resp]])
      if (length(ev) == 2L) pdf_$grp <- base[[ev[2L]]]
      attr(df, "points_df") <- ce_points_by_cond(pdf_, base, cond_sets,
                                                 ev, groups)
    }
    out[[key]] <- df
  }
  structure(out, class = "frmtmb_conditional_effects")
}

#' @rdname conditional_effects
#' @exportS3Method brms::conditional_effects
#' @export
conditional_effects.frmtmb_fit <- function(x, effects = NULL, resp = NULL,
                                           dpar = NULL, resolution = 100,
                                           prob = 0.95,
                                           method = c("epred", "predict"),
                                           band = c("wald", "profile",
                                                    "boot"),
                                           re_formula = NA,
                                           ndraws = 400, boot = NULL,
                                           profile_points = 25,
                                           seed = NULL,
                                           conditions = list(),
                                           surface = FALSE,
                                           data = NULL,
                                           int_conditions = list(),
                                           categorical = NULL, ...) {
  method <- ce_method(method)
  band <- match.arg(band)
  # prob becomes a normal quantile that RECYCLES along the grid, so a
  # length-2 prob drew a band whose coverage alternated point by point.
  # A length-2 resolution silently used only its first element.
  check_probability(prob, "prob")
  check_count(resolution, "resolution", min = 1L)
  check_count(ndraws, "ndraws", min = 1L)
  check_count(profile_points, "profile_points", min = 1L)
  check_flag(surface, "surface")
  check_named_list(conditions, "conditions", "conditions = list(z = 0)")
  check_named_list(int_conditions, "int_conditions",
                   "int_conditions = list(z = c(-1, 0, 1))")
  dots <- list(...)
  re_formula <- ce_re_formula(re_formula, dots)
  # the dots used to go straight to predict(), so an argument this
  # function does not know was reported against a function the user
  # never called - and only on the branches that forward, which left
  # the ordinal display discarding arguments in silence
  allow_new_levels <- ce_dots(dots)
  pop_level <- !inherits(re_formula, "formula") &&
    length(re_formula) == 1L && is.na(re_formula)
  resp <- resp %||% names(x$spec$responses)[1L]
  rspec <- x$spec$responses[[resp]]
  ce_structure_check(rspec)
  if (isTRUE(surface)) {
    stop("conditional_effects(surface = TRUE) is not implemented: the ",
         "display draws curves with bands, not a fitted surface. Ask ",
         "for the two-variable effect instead, e.g. ",
         "effects = \"x1:x2\", which varies x1 over its range at three ",
         "values of x2 (or at its levels)", call. = FALSE)
  }
  # a per-category effect display is on the CATEGORIES, not the latent
  # scale; naming a dpar explicitly is the way back to the predictor
  kind <- ce_display_kind(rspec, dpar, categorical)
  categorical <- identical(kind, "cats")
  cats_mean <- identical(kind, "cats_mean")
  if ((categorical || cats_mean) && method == "predict") {
    stop("method = \"predict\" has no meaning on an ordinal family: the ",
         "category probabilities conditional_effects() draws ARE the ",
         "predictive distribution, so there is no further observation ",
         "noise to add. Use method = \"epred\" (the default), or ask ",
         "for the latent predictor with dpar = \"mu\"", call. = FALSE)
  }
  # the delta method for a category probability runs through the ordinal
  # THRESHOLDS (ord_prob_se); a nominal family has none, so its bands
  # come from refits until someone writes that Jacobian
  if (categorical && band != "boot" &&
      !identical(rspec$family[["type"]], "ordinal")) {
    stop("conditional_effects() has no analytic standard error for the ",
         "category probabilities of family '", rspec$family[["family"]],
         "': the delta method it uses is written for ordinal thresholds. ",
         "Use band = \"boot\"", call. = FALSE)
  }
  if (band != "wald" && method == "predict") {
    stop("band = \"", band, "\" does not apply to method = \"predict\": ",
         "a prediction interval is already a quantile of simulated ",
         "responses, not a Wald band, and adding parameter uncertainty ",
         "to it twice is not an interval for anything. Use ",
         "method = \"epred\"", call. = FALSE)
  }
  # properties of the FAMILY, so they are settled before any grid is
  # built: raised from inside the loop they arrived after the point
  # estimate had already been computed, and a family whose expected
  # response is not one number per row died there on a replacement
  # length instead
  if (method == "predict") {
    if (!sim_can(rspec$family)) {
      stop("method = 'predict' needs a family with a simulator",
           sim_note(rspec$family), call. = FALSE)
    }
    # the band is built from ROWWISE draws, so a family that only draws
    # whole is refused even though simulate() accepts it:
    # mixture_mvn()'s sim_ctx returns an n by D matrix, which the
    # per-row quantiles cannot consume
    if (is.null(rspec$family[["sim"]])) {
      stop("method = 'predict' needs a family that draws row by row; '",
           rspec$family[["family"]], "' draws the response whole",
           call. = FALSE)
    }
  }
  dpar_given <- !is.null(dpar)
  dpar <- dpar %||% if ("mu" %in% names(rspec$dpars)) "mu" else
    rspec$primary_dpars[1]
  lp <- find_linpred(x, resp, dpar)
  # the display quantity when no dpar is named and the family's mean is
  # not the inverse link of mu: the EXPECTED RESPONSE, which is what
  # method = "epred" has always been documented to draw and what
  # fitted() and predict(type = "response") return. Taking the estimate
  # from the mu predictor alone plotted (1 - zi) times too little on a
  # zero-inflated fit and a sign-changing error on a hurdle one.
  pred_dpar <- ce_pred_dpar(rspec, dpar, dpar_given, categorical,
                            cats_mean)
  # `dpar` is resolved by now, so NULL here means, and only means, that
  # the display quantity is the expected response
  mean_display <- is.null(pred_dpar)
  # a dpar whose RESPONSE scale is not its own link inverse (a
  # mixture's mixing weight, which is a softmax over the component
  # predictors) is displayed on that scale, not on the predictor the
  # band happens to be computed from
  hook <- if (!categorical && !cats_mean) {
    dpar_report_hook(rspec$family, dpar, rspec)
  }
  if (band == "profile") {
    ce_profile_check(x, rspec, lp, dpar_given, categorical || cats_mean)
    if (!pop_level) {
      stop("band = \"profile\" draws its band over the OUTER ",
           "parameters, so it exists only for the population-level ",
           "curve (re_formula = NA). Use band = \"wald\" or ",
           "\"boot\" to condition on random effects", call. = FALSE)
    }
  }
  if ((categorical || cats_mean) && !pop_level &&
      identical(rspec$family[["type"]], "ordinal") && band != "boot") {
    stop("the ordinal per-category delta method is written for the ",
         "population-level curve; with re_formula use band = ",
         "\"boot\", or ask for dpar = \"mu\"", call. = FALSE)
  }
  # a nonlinear predictor has no delta-method standard error, so it
  # needs a band that never asks for one: boot refits, predict simulates
  if (!is.null(lp[["nl_body"]]) && band != "boot" && method != "predict") {
    stop("conditional_effects() cannot put a ", band, " band on a ",
         "nonlinear predictor: predict() has no standard error for it. ",
         "Use band = \"boot\", which refits instead of differentiating, ",
         "method = \"predict\", whose band is a quantile of simulated ",
         "responses, or display one nonlinear parameter with dpar = \"",
         rspec$nlpars[1L], "\"", call. = FALSE)
  }
  # every grid of the call is built before any band is: one bootstrap
  # covers all of them, which is the whole point of doing it here rather
  # than per effect
  # re_formula = NULL keeps the random effects, and the frequentist
  # answer for "the group this curve belongs to" is a NEW one: its
  # modes are zero, so the curve is the population one and the band
  # carries the random-effect variance on top (lp_extra_var() adds it
  # for an unobserved level). brms draws a new group's effects per
  # posterior draw instead, which is the paradigm difference; what
  # frmtmb did before was to take the FIRST OBSERVED level silently.
  na_vars <- if (!pop_level) ce_group_vars(x) else character(0)
  gb <- ce_grids_build(x, rspec, lp, effects, resp, dpar, resolution,
                       conditions, data, int_conditions, na_vars)
  base <- gb$base
  effects <- gb$effects
  cond_sets <- gb$cond_sets
  grids <- gb$grids
  # a new level is not in the fit's factor levels, so the design builder
  # has to be told it may meet one
  anl <- allow_new_levels || length(setdiff(na_vars, names(conditions)))
  z <- stats::qnorm(1 - (1 - prob) / 2)
  bd <- if (band == "boot") {
    # a new group's effects are DRAWN per replicate rather than left at
    # their zero modes; without that the interval collapsed onto the
    # population one while the frame's NA grouping column claimed
    # otherwise
    nspec <- if (length(na_vars)) {
      ce_new_level_spec(x, na_vars, base)
    } else {
      list()
    }
    bgrids <- if (length(nspec)) ce_boot_grids(grids, nspec, base) else grids
    ce_boot_draws(x, bgrids, categorical, resp, pred_dpar, boot, seed,
                  re_formula, anl, nspec)
  }
  pfail <- c(0L, 0L)

  dfs_by_eff <- list()
  for (gi in seq_along(grids)) {
    g <- grids[[gi]]
    ev <- g$ev
    nd <- g$nd
    n <- g$n
    ci <- g$ci
    cset <- g$cset
    bcols <- if (band == "boot") {
      bd$bs$t[, bd$offsets[gi] + seq_len(bd$lens[gi]), drop = FALSE]
    }
    # brms's frame always carries cond__, with one level when there is
    # one condition set: its own plot() facets on the column, so a
    # ported faceting call needs it there rather than only when
    # conditions = was passed
    clev <- names(cond_sets) %||% as.character(seq_along(cond_sets))
    cond <- factor(clev[ci], levels = clev)
    if (categorical) {
      if (identical(rspec$family[["type"]], "ordinal")) {
        ed <- lp_eta_design(x, lp, nd, !pop_level, anl)
        ps <- ord_prob_se(x, rspec, lp, ed, nd, !pop_level)
      } else {
        # a nominal family has no thresholds, so the ordinal delta
        # method does not apply; under band = "boot" (the only band
        # allowed here) the draws supply the se and the bounds
        P <- predict(x, newdata = nd, type = "response", resp = resp,
                     re.form = re_formula, allow_new_levels = anl)
        ps <- list(P = P, se = matrix(NA_real_, nrow(P), ncol(P)))
      }
      cats <- colnames(ps$P)
      df <- do.call(rbind, lapply(seq_along(cats), function(k) {
        d <- ce_frame(nd, ev, g$v2, cond,
                      cats = factor(cats[k], levels = cats))
        pk <- ps$P[, k]
        sk <- ps$se[, k]
        d$estimate__ <- pk
        d$se__ <- sk
        # the band is a Wald interval on the LOGIT of the probability:
        # on the probability scale itself it would leave [0, 1] near a
        # category that is nearly certain or nearly impossible
        sl <- sk / pmax(pk * (1 - pk), .Machine$double.eps)
        d$lower__ <- stats::plogis(stats::qlogis(pk) - z * sl)
        d$upper__ <- stats::plogis(stats::qlogis(pk) + z * sl)
        d
      }))
    } else if (cats_mean) {
      # brms's ordinal default: the expected CATEGORY NUMBER,
      # sum_k k p_k. One quantity per row, so its delta-method standard
      # error is the per-category one with the same category weights
      # applied to the gradient before the quadratic form.
      ed <- lp_eta_design(x, lp, nd, !pop_level, anl)
      ps <- ord_prob_se(x, rspec, lp, ed, nd, !pop_level,
                        weights = seq_len(ordinal_ncat(x)))
      df <- ce_frame(nd, ev, g$v2, cond)
      df$estimate__ <- as.vector(ps$P)
      df$se__ <- as.vector(ps$se)
      df$lower__ <- df$estimate__ - z * df$se__
      df$upper__ <- df$estimate__ + z * df$se__
    } else {
      df <- ce_frame(nd, ev, g$v2, cond)
      if (band == "boot" || method == "predict") {
        # the draws (or simulated responses) are the band AND the
        # standard error here, so the delta method is not asked for:
        # that is what lets a nonlinear predictor, which has no
        # analytic se, reach a band at all
        df$estimate__ <- as.vector(predict(x, newdata = nd,
                                           type = "response",
                                           dpar = pred_dpar,
                                           resp = resp,
                                           re.form = re_formula,
                                           allow_new_levels = anl))
        df$se__ <- NA_real_
        df$lower__ <- NA_real_
        df$upper__ <- NA_real_
      } else if (mean_display || !is.null(hook)) {
        if (mean_display) {
          # the mean runs through the addition terms as well as the
          # dpars, so a term whose value on a grid row would be a
          # reference value rather than a real one is refused here, on
          # the same rule (and with the same message) method =
          # "predict" has always used
          ce_aterms(rspec, nd, cset, n)
        }
        # the expected response (or the reported dpar) and ITS standard
        # error: for the mean the delta method runs over every dpar's
        # linear predictor jointly (predict_mean_se()), so the
        # cross-dpar covariances are in the band rather than dropped
        p <- predict(x, newdata = nd, type = "response", dpar = pred_dpar,
                     resp = resp, re.form = re_formula, se.fit = TRUE,
                     allow_new_levels = anl)
        # a reported probability gets a logit band, which cannot leave
        # (0, 1); anything else keeps the predictor's own link
        bl <- if (!is.null(hook) && all(is.finite(p$fit)) &&
                    all(p$fit > 0 & p$fit < 1)) {
          get_link("logit")
        } else {
          lp[["link"]]
        }
        bs <- ce_band_scale(bl, p$fit)
        df$estimate__ <- p$fit
        # abs(): a decreasing link (1/mu) has a negative derivative,
        # and a standard error is not negative
        df$se__ <- p$se.fit / abs(bs$dm)
        df$lower__ <- bs$inv(bs$t - z * df$se__)
        df$upper__ <- bs$inv(bs$t + z * df$se__)
        lo <- pmin(df$lower__, df$upper__)
        df$upper__ <- pmax(df$lower__, df$upper__)
        df$lower__ <- lo
      } else {
        p <- predict(x, newdata = nd, type = "link", dpar = dpar,
                     resp = resp, re.form = re_formula, se.fit = TRUE,
                     allow_new_levels = anl)
        df$estimate__ <- lp[["link"]]$linkinv(p$fit)
        df$se__ <- p$se.fit
        df$lower__ <- lp[["link"]]$linkinv(p$fit - z * p$se.fit)
        df$upper__ <- lp[["link"]]$linkinv(p$fit + z * p$se.fit)
      }
      if (method == "predict") {
        fam <- rspec$family
        # the NATURAL scale, which is the one the density and the
        # simulator consume: a mixture's theta reports as a softmax
        # probability on the response scale and would reach log_pi()
        # already normalized
        dpv <- dpars_natural(x, rspec, nd, re_formula, anl)
        avc <- ce_aterms(rspec, nd, cset, n)
        # sim_response(), not fam$sim(): trunc() bounds are respected by
        # rejection, as everywhere else responses are drawn
        sims <- replicate(ndraws, sim_response(fam, dpv, avc, n,
                                               extra = fit_extras(x)))
        # the point estimate moves onto the response scale the bands
        # live on: a binomial band is a count, not a probability, and a
        # truncated band is centered on the truncated mean
        df$estimate__ <- response_mean(fam, dpv, avc)
        df$lower__ <- apply(sims, 1, stats::quantile, (1 - prob) / 2)
        df$upper__ <- apply(sims, 1, stats::quantile, 1 - (1 - prob) / 2)
        df$se__ <- apply(sims, 1, stats::sd)
      }
    }
    # the estimate is the fit's; only the band changes with `band`
    if (band == "boot") {
      df$lower__ <- ce_pctl(bcols, (1 - prob) / 2)
      df$upper__ <- ce_pctl(bcols, 1 - (1 - prob) / 2)
      df$se__ <- apply(bcols, 2, function(v) {
        if (sum(is.finite(v)) < 2L) NA_real_ else stats::sd(v, na.rm = TRUE)
      })
    } else if (band == "profile") {
      pci <- ce_profile_eta_ci(x, lp, nd, g$v1, g$n1, g$n2, prob,
                               profile_points)
      # the link is monotone but not necessarily increasing (inverse,
      # 1/mu), so the endpoints are sorted after the transform
      lo <- lp[["link"]]$linkinv(pci$lower)
      up <- lp[["link"]]$linkinv(pci$upper)
      df$lower__ <- pmin(lo, up)
      df$upper__ <- pmax(lo, up)
      pfail <- pfail + c(pci$fails, pci$tried)
    }
    dfs_by_eff[[g$eff]] <- c(dfs_by_eff[[g$eff]], list(df))
  }
  if (pfail[1L] > 0L) {
    warning("band = \"profile\": the likelihood-root search did not ",
            "converge at ", pfail[1L], " of ", pfail[2L],
            " profiled grid point(s); their bounds are NA",
            call. = FALSE)
  }

  out <- ce_finalize(dfs_by_eff, effects, rspec, resp, dpar, band, base,
                     categorical, cond_sets, gb$groups, categorical)
  # the bootstrap rides along so a second call can reuse it: refits are
  # the expensive part and nobody should pay for them twice
  if (!is.null(bd)) attr(out, "boot") <- bd$bs
  out
}

#' Pointwise percentile of a draws matrix (draws in rows), NA where a
#' column has no usable draw rather than an error out of `quantile()`.
#'
#' @noRd
ce_pctl <- function(m, p) {
  apply(m, 2, function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) NA_real_ else unname(stats::quantile(v, p))
  })
}

#' @export
print.frmtmb_conditional_effects <- function(x, ...) {
  plot(x, ...)
  invisible(x)
}

#' @export
plot.frmtmb_conditional_effects <- function(x, ask = NULL, points = FALSE,
                                            ncol = NULL, ...) {
  if (!is.null(ncol)) check_count(ncol, "ncol", min = 1L)
  # a condition set is a FACET, not a page: several conditions used to
  # draw several full pages that overwrote each other on a normal
  # device, so the page prompt counts effects, not conditions
  ask <- ask %||% (length(x) > 1L && grDevices::dev.interactive())
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask), add = TRUE)
  }
  for (nm in names(x)) {
    df <- x[[nm]]
    if (points && is.null(attr(df, "points_df"))) {
      message("points = TRUE: no observations to draw for effect '", nm,
              "' (the display is per-category, on a non-mean ",
              "distributional parameter, or the response is not a ",
              "plain numeric column)")
    }
    if (!is.null(df$cond__) && length(unique(df$cond__)) > 1L) {
      ce_plot_facets(df, points = points, ncol = ncol)
    } else {
      ce_plot_one(df, points = points)
    }
  }
  invisible(x)
}

#' Panel grid for a faceted display. brms hands `ncol` straight to
#' facet_wrap(), whose NULL default lays the panels out roughly square;
#' asking for more columns than there are panels only wastes the page,
#' so the request is capped.
#'
#' @noRd
ce_facet_layout <- function(n, ncol = NULL) {
  nc <- if (is.null(ncol)) ceiling(sqrt(n)) else min(as.integer(ncol), n)
  nc <- max(1L, as.integer(nc))
  c(nrow = as.integer(ceiling(n / nc)), ncol = nc)
}

#' The raw observations belonging to one condition. A points frame with
#' no `cond__` column was never split (no condition variable named a
#' column of the model data), so every panel draws all of it, which is
#' what brms's ggplot layer does when the point frame lacks the facet
#' variable.
#'
#' @noRd
ce_points_at <- function(pts, cv) {
  if (is.null(pts) || is.null(pts[["cond__"]])) return(pts)
  out <- pts[as.character(pts$cond__) == cv, , drop = FALSE]
  if (!nrow(out)) NULL else out
}

#' One page of small multiples, one panel per condition set. tinyplot
#' does the shared axes and the outer legend natively; without it the
#' fallback is a par(mfrow) grid of the same base panels, which is still
#' ONE page honoring `ncol`.
#'
#' The per-category ordinal display already spends its grouping slot on
#' the category and its panel slot on a second predictor, so it takes
#' the base grid rather than a second facet dimension.
#'
#' @noRd
ce_plot_facets <- function(df, points = FALSE, ncol = NULL) {
  lv <- unique(as.character(df$cond__))
  lay <- ce_facet_layout(length(lv), ncol)
  pts <- if (points) attr(df, "points_df")
  if (is.null(df[["cats__"]]) &&
      requireNamespace("tinyplot", quietly = TRUE)) {
    ce_facet_tinyplot(df, lv, lay, pts)
  } else {
    ce_facet_base(df, lv, lay, points)
  }
  invisible(NULL)
}

#' @noRd
ce_facet_base <- function(df, lv, lay, points) {
  pts <- if (points) attr(df, "points_df")
  # the panels share one scale: a per-panel range would make curves of
  # different heights look alike, which is the opposite of what small
  # multiples are for
  ylim <- range(df$lower__, df$upper__, df$estimate__,
                if (!is.null(pts)) pts$y, na.rm = TRUE)
  op <- graphics::par(mfrow = c(lay[["nrow"]], lay[["ncol"]]),
                      mar = c(4, 4, 2.5, 1))
  on.exit(graphics::par(op), add = TRUE)
  for (cv in lv) {
    sub <- df[as.character(df$cond__) == cv, , drop = FALSE]
    for (a in c("effects", "response", "dpar")) {
      attr(sub, a) <- attr(df, a)
    }
    attr(sub, "points_df") <- ce_points_at(pts, cv)
    ce_plot_one(sub, cond = cv, points = points, ylim = ylim)
  }
}

#' @noRd
ce_facet_tinyplot <- function(df, lv, lay, pts) {
  ev <- attr(df, "effects")
  xv <- df[[ev[1L]]]
  ylim <- range(df$lower__, df$upper__, df$estimate__,
                if (!is.null(pts)) pts$y, na.rm = TRUE)
  args <- list(
    x = if (is.numeric(xv)) xv else factor(xv),
    y = df$estimate__, ymin = df$lower__, ymax = df$upper__,
    facet = factor(as.character(df$cond__), levels = lv),
    facet.args = list(ncol = lay[["ncol"]]),
    type = if (is.numeric(xv)) "ribbon" else "pointrange",
    xlab = ev[1L],
    ylab = paste0(attr(df, "response"), " (", attr(df, "dpar"), ")"),
    ylim = ylim,
    # tinyplot leaves the facet layout in par() otherwise, and the next
    # effect of the same call may be an ordinary single base panel
    restore.par = TRUE
  )
  if (length(ev) == 2L) {
    args$by <- factor(df[[ev[2L]]])
    args$legend <- list(title = ev[2L])
  }
  do.call(tinyplot::tinyplot, args)
  if (!is.null(pts) && nrow(pts)) {
    if (is.null(pts[["cond__"]])) {
      np <- nrow(pts)
      pts <- pts[rep(seq_len(np), times = length(lv)), , drop = FALSE]
      pts$cond__ <- rep(lv, each = np)
    }
    # tinyplot_add() inherits the previous layer's aesthetics, so the
    # band has to be cleared explicitly: the points carry none, and
    # there are usually more of them than there are grid rows, which
    # ends the call on a length mismatch. Clearing `by` the same way is
    # right only when the curves HAD one; asking to clear a grouping
    # that was never set walks off the end of tinyplot's palette
    add <- list(
      x = if (is.numeric(xv)) pts$x else factor(pts$x, levels = levels(args$x)),
      y = pts$y, facet = factor(as.character(pts$cond__), levels = lv),
      ymin = NULL, ymax = NULL,
      # one uniform translucent color, as in the single-panel display
      type = "p", pch = 16, cex = 0.5,
      col = grDevices::adjustcolor("black", 0.25)
    )
    # single-bracket assignment: `add$by <- NULL` would DROP the element
    # rather than pass a NULL that clears the inherited grouping
    if (!is.null(args$by)) add["by"] <- list(NULL)
    do.call(tinyplot::tinyplot_add, add)
  }
}

#' Draw one conditional-effect panel: the estimate over the varied
#' predictor with its band, lines and a shaded band for a numeric
#' predictor, points and error bars for a discrete one, split by the
#' optional second predictor.
#'
#' @noRd
ce_plot_one <- function(df, cond = NULL, points = FALSE, ylim = NULL) {
  pts <- if (points) attr(df, "points_df")
  ev <- attr(df, "effects")
  if (!is.null(df[["cats__"]])) {
    # an ordinal display carries one curve per response category, so the
    # category takes the grouping slot; a second predictor then needs a
    # panel of its own rather than a second set of colors
    ylab <- paste0("P(", attr(df, "response"), ")")
    if (!is.null(cond)) ylab <- paste0(ylab, " | ", cond)
    ylim <- ylim %||% range(df$lower__, df$upper__, na.rm = TRUE)
    # "cats__" is the category dimension itself, which already has the
    # grouping slot; only a real second PREDICTOR needs its own panel
    if (length(ev) == 2L && !identical(ev[2L], "cats__")) {
      for (lv in unique(df[[ev[2L]]])) {
        sub <- df[df[[ev[2L]]] == lv, , drop = FALSE]
        ce_draw_panel(sub, ev[1L], factor(sub$cats__,
                                          levels = levels(df$cats__)),
                      "category",
                      paste0(ylab, " | ", ev[2L], " = ", lv), ylim)
      }
    } else {
      ce_draw_panel(df, ev[1L], df$cats__, "category", ylab, ylim)
    }
    return(invisible(NULL))
  }
  ylab <- paste0(attr(df, "response"), " (", attr(df, "dpar"), ")")
  if (!is.null(cond)) ylab <- paste0(ylab, " | ", cond)
  grp <- if (length(ev) == 2L) factor(df[[ev[2L]]])
  ylim <- ylim %||% range(df$lower__, df$upper__,
                          if (!is.null(pts)) pts$y, na.rm = TRUE)
  ce_draw_panel(df, ev[1L], grp, ev[2L], ylab, ylim, pts = pts)
}

#' Draw one panel: the estimate over the varied predictor `xv` with its
#' band, lines and a shaded band for a numeric predictor, points and
#' error bars for a discrete one, split by the optional grouping factor
#' `grp`.
#'
#' @noRd
ce_draw_panel <- function(df, xv, grp, grp_title, ylab, ylim,
                          pts = NULL) {
  v1 <- df[[xv]]
  ev <- c(xv, grp_title)
  pt_col <- grDevices::adjustcolor("black", 0.25)

  if (is.numeric(v1)) {
    graphics::plot(range(v1), ylim, type = "n", xlab = ev[1L],
                   ylab = ylab)
    if (!is.null(pts)) {
      graphics::points(pts$x, pts$y, pch = 16, cex = 0.5, col = pt_col)
    }
    if (is.null(grp)) {
      graphics::polygon(c(v1, rev(v1)), c(df$lower__, rev(df$upper__)),
                        col = grDevices::adjustcolor("black", 0.15),
                        border = NA)
      graphics::lines(v1, df$estimate__, lwd = 2)
    } else {
      for (k in seq_along(levels(grp))) {
        i <- grp == levels(grp)[k]
        graphics::polygon(c(v1[i], rev(v1[i])),
                          c(df$lower__[i], rev(df$upper__[i])),
                          col = grDevices::adjustcolor(k, 0.12),
                          border = NA)
        graphics::lines(v1[i], df$estimate__[i], col = k, lwd = 2)
      }
      graphics::legend("topleft", legend = levels(grp), col =
                         seq_along(levels(grp)), lwd = 2, title = ev[2L],
                       bty = "n")
    }
  } else {
    xi <- as.integer(factor(v1))
    if (!is.null(grp)) {
      xi <- xi + (as.integer(grp) - (nlevels(grp) + 1) / 2) * 0.15
    }
    graphics::plot(range(xi) + c(-0.5, 0.5), ylim, type = "n",
                   xaxt = "n", xlab = ev[1L], ylab = ylab)
    graphics::axis(1, at = seq_len(nlevels(factor(v1))),
                   labels = levels(factor(v1)))
    if (!is.null(pts)) {
      xp <- as.integer(factor(pts$x, levels = levels(factor(v1))))
      # deterministic spread, no RNG: replotting looks identical and the
      # user's random seed is left alone
      off <- ((seq_along(xp) * 7L) %% 17L - 8L) / 100
      graphics::points(xp + off, pts$y, pch = 16, cex = 0.5,
                       col = pt_col)
    }
    cols <- if (is.null(grp)) 1L else as.integer(grp)
    graphics::arrows(xi, df$lower__, xi, df$upper__, angle = 90,
                     code = 3, length = 0.05, col = cols)
    graphics::points(xi, df$estimate__, pch = 16, col = cols)
    if (!is.null(grp)) {
      graphics::legend("topleft", legend = levels(grp),
                       col = seq_along(levels(grp)), pch = 16,
                       title = ev[2L], bty = "n")
    }
  }
}

#' Diagnostic plots for a fit
#'
#' Panel 1: Pearson residuals against fitted values with a lowess
#' trend. Panel 2: normal QQ plot of the Pearson residuals. For
#' simulation-based residuals that are exact for discrete families, use
#' [dharma_residuals()] or `residuals(type = "osa")`.
#'
#' On an ordinal fit [fitted()] is a matrix of category probabilities,
#' so panel 1 uses the expected category index `sum_k k * P(y = k)` -
#' the same scalar the Pearson residual is taken against - and labels
#' the axis accordingly.
#'
#' @param x A `frmtmb_fit`.
#' @param which Subset of `1:2`.
#' @param ask Whether to prompt between plots; defaults to the usual
#'   interactive-device rule.
#' @param ... Unused.
#' @return `x`, invisibly. Called for the plots it draws.
#'
#' @srrstats {RE6.0} A `frmtmb_fit` has a default `plot()` method, so
#'   `plot(fit)` works without the user naming a function. It draws the
#'   two standard regression diagnostics: Pearson residuals against
#'   fitted values with a lowess trend, and a normal QQ plot of those
#'   residuals.
#' @srrstats {RE6.2} The first panel of the default `plot()` method plots
#'   the fitted values (on the horizontal axis, against the Pearson
#'   residuals), so the model's fitted response is visualized by default.
#'   Bands for the fitted response are available through
#'   [conditional_effects()], which draws Wald or predictive intervals.
#' @srrstats {RE6.1} The method is a real S3 method dispatched on the
#'   class of the returned object (`plot.frmtmb_fit`), registered in
#'   `NAMESPACE`, so the generic reaches it and no separate signposting
#'   is needed. [conditional_effects()] and [pp_check()] have their own
#'   plot methods for effect displays and posterior-predictive checks,
#'   and this page points at [dharma_residuals()] and
#'   `residuals(type = "osa")` for residuals that stay exact under
#'   discrete families.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#' # both panels, side by side
#' op <- par(mfrow = c(1, 2))
#' plot(fit, ask = FALSE)
#' par(op)
#'
#' # just the QQ panel
#' plot(fit, which = 2)
#' @export
plot.frmtmb_fit <- function(x, which = 1:2, ask = NULL, ...) {
  r <- residuals(x, type = "pearson")
  ask <- ask %||% (length(which) > 1L && grDevices::dev.interactive())
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask), add = TRUE)
  }
  if (1L %in% which) {
    # fitted() is a K-column probability matrix on an ordinal fit, and
    # the panel needs one number per row: the expected category index,
    # which is what the residual on the vertical axis was taken against
    rspec <- x$spec$responses[[1L]]
    ordinal <- identical(rspec$family[["type"]], "ordinal")
    ft <- if (ordinal) {
      napred(x, ord_cat_moments(x, rspec)$mean)
    } else {
      fitted(x)
    }
    graphics::plot(ft, r,
                   xlab = if (ordinal) "Expected category" else
                     "Fitted values",
                   ylab = "Pearson residuals")
    graphics::abline(h = 0, lty = 2)
    ok <- is.finite(ft) & is.finite(r)
    graphics::lines(stats::lowess(ft[ok], r[ok]), col = 2, lwd = 2)
  }
  if (2L %in% which) {
    stats::qqnorm(r, main = "Pearson residuals")
    stats::qqline(r, lty = 2)
  }
  invisible(x)
}

#' Predictive check against simulated responses
#'
#' The frequentist analog of brms's `pp_check()`: responses are
#' simulated from the fitted model (marginally over the random effects)
#' and handed to the corresponding bayesplot `ppc_*` function
#' (bayesplot must be installed, but not necessarily attached).
#'
#' @param object A `frmtmb_fit` for a univariate model.
#' @param ... Passed to the `ppc_*` function.
#' @return A ggplot object, as returned by the bayesplot `ppc_*`
#'   function that `type` selects.
#' @examples
#' if (requireNamespace("bayesplot", quietly = TRUE)) {
#'   set.seed(1)
#'   dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#'   dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
#'   fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
#'
#'   # the observed density against draws from the fit
#'   pp_check(fit, ndraws = 20)
#'
#'   # any bayesplot ppc_* check, named by its suffix. A statistic the
#'   # model was not fitted to is the informative one: here, the share
#'   # of zeros, which is how zero inflation shows up.
#'   pp_check(fit, type = "stat", stat = function(y) mean(y == 0),
#'            ndraws = 50)
#' }
#' @export
pp_check <- function(object, ...) {
  # frmtmb ships its own generic so pp_check(fit) works without
  # attaching bayesplot; the methods are ALSO registered on
  # bayesplot::pp_check, so whichever generic sits in front on the
  # search path dispatches to the same code
  UseMethod("pp_check")
}

#' @rdname pp_check
#' @param type The bayesplot check, i.e. the part after `ppc_`
#'   (`"dens_overlay"`, `"hist"`, `"stat"`, `"scatter_avg"`, ...).
#' @param ndraws Number of simulated response vectors.
#' @param re_formula The random-effect switch, in brms's spelling
#'   (`pp_check()` is a brms function). On a fit it is passed to
#'   [simulate()] and defaults to `NA`, which simulates new random
#'   effects; on draws it is passed to `posterior_predict()` and
#'   defaults to `NULL`, because a draw already carries its own.
#' @param re.form lme4's spelling of `re_formula`, accepted as an alias.
#'   Pass one or the other, not both; see the *Argument spellings*
#'   section of `frmtmb.sample::posterior_epred()`.
#' @exportS3Method bayesplot::pp_check
#' @export
pp_check.frmtmb_fit <- function(object, type = "dens_overlay",
                                ndraws = 10,
                                re_formula = arg_unset(),
                                re.form = arg_unset(), ...) {
  # NA, not NULL: a fit has ONE estimate of the random effects, so
  # conditioning on it would compare the data against draws that already
  # know each group's deviation. New levels per replicate are what makes
  # this the frequentist analog of the posterior predictive check.
  re_form <- re_form_arg(re_formula, re.form, "pp_check()", default = NA)
  rspec <- single_response(object, "pp_check()")
  y <- object$frame[["y"]][[1L]]
  if (is.matrix(y)) {
    stop("pp_check() on a fit supports vector responses", call. = FALSE)
  }
  sims <- na_unpad(object, simulate(object, nsim = ndraws,
                                    re.form = re_form))
  # ordinal draws come back as ordered factors carrying the response's
  # levels; bayesplot compares them with y, which is the 1..K codes
  yrep <- if (identical(rspec$family[["type"]], "ordinal")) {
    matrix(unlist(lapply(sims, as.integer), use.names = FALSE),
           nrow = nrow(sims))
  } else {
    as.matrix(sims)
  }
  yrep <- t(yrep)
  fun <- get(paste0("ppc_", type), envir = asNamespace("bayesplot"))
  fun(as.numeric(y), yrep, ...)
}
