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
      stop("conditional_effects(method = \"predict\") cannot evaluate ",
           aterm_label(nm, ex), " on the effect grid: its value would ",
           "be a reference value, not a real one. Pin ",
           paste(setdiff(vars, names(cset)), collapse = ", "),
           " in conditions = list(...).", call. = FALSE)
    }
    v <- tryCatch(as.numeric(eval(ex, nd, rspec$formula_env)),
                  error = function(e) NULL)
    if (!is.null(v) && !length(v) %in% c(1L, n)) v <- NULL
    if (is.null(v) && nm %in% strict) {
      stop("conditional_effects(method = \"predict\") could not ",
           "evaluate ", aterm_label(nm, ex), " on the effect grid",
           call. = FALSE)
    }
    if (!is.null(v)) av[[nm]] <- v
  }
  if (!is.null(rspec$aterms$se_sigma)) av$se_sigma <- rspec$aterms$se_sigma
  av
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
#' @noRd
ce_grid_values <- function(col, resolution, nm = "the predictor") {
  if (is.factor(col)) {
    factor(levels(col), levels = levels(col))
  } else if (is.numeric(col)) {
    if (!any(is.finite(col))) {
      stop("Variable '", nm, "' has no finite values to build an effect ",
           "grid from", call. = FALSE)
    }
    seq(min(col, na.rm = TRUE), max(col, na.rm = TRUE),
        length.out = resolution)
  } else if (is.logical(col)) {
    c(FALSE, TRUE)
  } else {
    sort(unique(col))
  }
}

#' Values for the second predictor of an `"x:z"` effect.
#'
#' @noRd
ce_second_values <- function(col) {
  if (is.numeric(col) && !is.matrix(col)) {
    signif(mean(col, na.rm = TRUE) +
             c(-1, 0, 1) * stats::sd(col, na.rm = TRUE), 3)
  } else {
    ce_grid_values(col, resolution = 0)
  }
}

#' Plottable variables of ONE linear predictor: its fixed-effect terms,
#' its smooth terms, its monotonic terms and its `mi()` terms. An `mo()`
#' or `mi()` column is a placeholder in the design matrix and its
#' variable never reaches `terms`, so what is stored with the term is
#' the only place the name survives.
#'
#' @noRd
ce_lp_vars <- function(lp) {
  v <- if (is.null(lp$terms)) {
    character(0)
  } else {
    all.vars(stats::delete.response(lp$terms))
  }
  v <- c(v, unlist(lapply(lp$smooths, function(si) si$sm$term)))
  for (m in lp$mo %||% list()) {
    v <- c(v, all.vars(m$expr), all.vars(m$mult_expr))
  }
  for (m in lp$mi %||% list()) {
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
#' @noRd
ce_plot_vars <- function(x, rspec, lp, resp) {
  v <- ce_lp_vars(lp)
  if (!is.null(lp$nl_body)) {
    v <- c(v, names(lp$data_list),
           setdiff(all.vars(lp$nl_body), rspec$nlpars))
    for (np in rspec$nlpars) {
      lpn <- x$frame$linpreds[[linpred_key(resp, np)]]
      if (!is.null(lpn)) v <- c(v, ce_lp_vars(lpn))
    }
  }
  unique(v)
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
    isTRUE(rspec$family$type %in% c("ordinal", "categorical"))
}

#' One conditional-effects grid: every predictor at its reference value
#' or at its `conditions` override, with the varied predictor(s) replaced
#' by their grid values. Split out because `band = "boot"` needs every
#' grid of the call BEFORE any refit happens, so that one bootstrap can
#' serve all of them.
#'
#' @noRd
ce_build_nd <- function(base, ev, v1, v2, cset, n, n2) {
  nd <- data.frame(.ce_row = seq_len(n))
  for (nm in names(base)) {
    val <- if (nm %in% names(cset)) {
      cnd <- cset[[nm]]
      if (is.factor(base[[nm]])) {
        factor(cnd, levels = levels(base[[nm]]))
      } else {
        cnd
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
  nd[[ev[1L]]] <- rep(v1, times = n2)
  if (length(ev) == 2L) nd[[ev[2L]]] <- rep(v2, each = length(v1))
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
ce_boot_one <- function(fit, nd, categorical, resp, dpar) {
  p <- if (categorical) {
    predict(fit, newdata = nd, type = "response", resp = resp,
            re.form = NA)
  } else {
    predict(fit, newdata = nd, type = "response", dpar = dpar,
            resp = resp, re.form = NA)
  }
  as.vector(p)
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
ce_boot_key <- function(grids, categorical, resp, dpar, lens) {
  serialize(list(
    resp = resp, dpar = dpar, categorical = categorical, lens = lens,
    grids = lapply(grids, function(g) {
      list(eff = g$eff, ci = g$ci, cset = g$cset, nd = as.list(g$nd))
    })
  ), NULL, xdr = FALSE)
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
ce_boot_draws <- function(x, grids, categorical, resp, dpar, boot, seed) {
  lens <- vapply(grids, function(g) {
    length(ce_boot_one(x, g$nd, categorical, resp, dpar))
  }, 1L)
  tot <- sum(lens)
  key <- ce_boot_key(grids, categorical, resp, dpar, lens)
  if (inherits(boot, "frmtmb_boot")) {
    if (!identical(boot$ce_key, key)) {
      stop("boot = was not produced by a conditional_effects(band = ",
           "\"boot\") call on this grid: its draws are ",
           if (is.null(boot$ce_key)) {
             "coefficients or another quantity"
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
    v <- tryCatch(
      unlist(lapply(grids, function(g) {
        ce_boot_one(f, g$nd, categorical, resp, dpar)
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
  list(bs = bs, lens = lens, offsets = cumsum(c(0L, lens)))
}

#' Component label per position of the outer (optimized) parameter
#' vector, in `obj$par` order. `outer_par_names()` gives the names; the
#' profile band needs to know which block each position belongs to so a
#' design row can be written as a `lincomb` over `obj$par`.
#'
#' @noRd
ce_outer_comp <- function(fit) {
  tpl <- fit$frame$par_template
  random <- c("b", "miss")
  if (fit$REML || isTRUE(fit$control$profile)) random <- c(random, "beta")
  comp <- character(0)
  for (cp in names(tpl)) {
    if (cp %in% random) next
    len <- length(tpl[[cp]])
    if (cp == "betad" && length(fit$frame$betad_fixed_idx)) {
      len <- len - length(fit$frame$betad_fixed_idx)
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
  if (!is.null(lp$nl_body)) {
    stop("band = \"profile\" is not available for a nonlinear ",
         "predictor: a grid value is not a linear combination of the ",
         "nonlinear parameters. Use band = \"boot\"", call. = FALSE)
  }
  if (!dpar_given && (!mean_is_mu(rspec$family) || has_trunc(rspec))) {
    stop("band = \"profile\" cannot cover the expected response of ",
         "family '", rspec$family$family, "': it runs through more ",
         "than one distributional parameter (zero inflation, a hurdle, ",
         "a dispersion) or through truncation bounds, so it is not a ",
         "single linear combination of the parameters. Use ",
         "band = \"boot\", or name one predictor with dpar =",
         call. = FALSE)
  }
  key <- linpred_key(lp$resp, lp$dpar)
  sm <- Filter(function(bk) {
    bk$covstruct %in% c("smooth", "gp", "hsgp") &&
      any(vapply(bk$components, function(cp) cp$lp_key == key, TRUE))
  }, x$frame$re_blocks)
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
  X <- as.matrix(ed$X)
  if (lp$par == "beta") {
    pos <- which(comp == "beta")[lp$idx]
    keep <- rep(TRUE, ncol(X))
  } else {
    tpl_len <- length(x$frame$par_template$betad)
    est_rank <- match(lp$idx, setdiff(seq_len(tpl_len),
                                      x$frame$betad_fixed_idx))
    keep <- !is.na(est_rank)
    pos <- which(comp == "betad")[est_rank[keep]]
  }
  X <- X[, keep, drop = FALSE]
  if (length(comp) != length(par) || anyNA(pos) || length(pos) != ncol(X)) {
    stop("band = \"profile\" cannot line this fit's coefficients up ",
         "with its outer parameter vector (a mapped, fixed or profiled ",
         "coefficient block). Use band = \"boot\"", call. = FALSE)
  }
  eta <- ed$eta
  n <- length(eta)
  target <- 0.5 * stats::qchisq(prob, df = 1)
  sel <- if (is.numeric(v1) && !is.matrix(v1) && n1 > profile_points) {
    unique(round(seq(1, n1, length.out = max(2L, profile_points))))
  } else {
    seq_len(n1)
  }
  rows <- as.vector(outer(sel, (seq_len(n2) - 1L) * n1, "+"))
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
      off <- (j - 1L) * n1
      fill <- function(y) {
        if (sum(is.finite(y)) < 2L) return(rep(NA_real_, n1))
        # na.rm = FALSE: an interval touching a failed point stays NA
        # rather than being bridged over silently
        stats::approx(xs[sel], y, xout = xs, na.rm = FALSE)$y
      }
      lo[off + seq_len(n1)] <- fill(lo[off + sel])
      up[off + seq_len(n1)] <- fill(up[off + sel])
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
#'   variable of the selected linear predictor.
#' @param resp,dpar Response and distributional parameter, as in
#'   [predict.frmtmb_fit()].
#' @param resolution Number of grid points for a varied numeric
#'   predictor.
#' @param prob Coverage of the confidence bands (brms spelling).
#' @param band How the confidence band is built: `"wald"` (default,
#'   the delta method on the link scale), `"profile"` (likelihood-root
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
#'   response. `"predict"`: prediction intervals - quantile bands from
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
#'   `method = "predict"`.
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
#' @param ... Passed to [predict.frmtmb_fit()].
#' @return A named list of data frames (one per effect) with the varied
#'   variable(s) plus `estimate__`, `se__` (link scale), `lower__`, and
#'   `upper__`; printing it draws the plots. An ordinal fit adds a
#'   `cats__` column and one block of rows per response category.
#'   `plot(ce, points = TRUE)` overlays the raw observations (the brms
#'   argument): all observations are shown regardless of `conditions`,
#'   and no points are drawn for a per-category ordinal display, a
#'   non-mean `dpar`, or a matrix response (a message says so).
#' @section Ordinal responses:
#' `cumulative()`, `sratio()`, `cratio()` and `acat()` have no mean, so
#' the display is per CATEGORY, as brms's `categorical = TRUE` is: each
#' effect data frame gains a `cats__` factor of the response's own
#' levels and carries the fitted category probability in `estimate__`,
#' with one curve per category in the plot (a second predictor gets a
#' panel of its own).
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
#' - `"wald"` (default, and free): the delta-method standard error on
#'   the link scale, back-transformed. Symmetric on the link scale by
#'   construction.
#' - `"profile"`: one likelihood-root search ([TMB::tmbroot()]) per grid
#'   point. A grid point's linear predictor is a linear combination of
#'   the coefficients, so the search inverts the likelihood ratio along
#'   that combination; the link is monotone, so the two endpoints map
#'   through it. The band is not symmetric and does not assume a
#'   quadratic log-likelihood, which is what makes it worth its cost
#'   near a boundary or at a small sample size.
#' - `"boot"`: percentiles of the grid predictions across the refits of
#'   ONE [frm_bootstrap()] (`re.form = NA`, the same population-level
#'   grid). One bootstrap serves every effect, condition set and ordinal
#'   category of the call; `attr(ce, "boot")` returns it, and passing it
#'   back as `boot =` costs no refits at all.
#'
#' What `se__` means follows the band: the Wald standard error (link
#' scale, or the probability scale on the ordinal display) for `"wald"`
#' and `"profile"` - the profile changes the endpoints, not the standard
#' error - and the standard deviation of the bootstrap draws, on the
#' displayed scale, for `"boot"`.
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
#' @section Which predictors are plotted by default:
#' Every variable of the selected linear predictor that the display can
#' vary: its fixed-effect terms, its smooth terms and its `mo()` terms
#' (whose design columns are placeholders, so the variable is read from
#' the term itself). On a nonlinear fit they live one level down - the
#' covariates the nonlinear formula reads, plus the terms of each
#' nonlinear parameter's own predictor - and all of those are collected
#' too. Matrix-valued columns are excluded, a grid over a matrix
#' covariate not being a curve. Naming `effects =` overrides the search.
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

#' @rdname conditional_effects
#' @export
conditional_effects.frmtmb_fit <- function(x, effects = NULL, resp = NULL,
                                           dpar = NULL, resolution = 100,
                                           prob = 0.95,
                                           method = c("epred", "predict"),
                                           band = c("wald", "profile",
                                                    "boot"),
                                           ndraws = 400, boot = NULL,
                                           profile_points = 25,
                                           seed = NULL,
                                           conditions = list(),
                                           surface = FALSE,
                                           data = NULL, ...) {
  method <- match.arg(method)
  band <- match.arg(band)
  resp <- resp %||% names(x$spec$responses)[1L]
  rspec <- x$spec$responses[[resp]]
  if (!is.null(rspec$family[["hmm"]])) {
    stop("conditional_effects() is not available for an hmm() fit: the ",
         "expected response weights the state means by posterior state ",
         "occupancies, which depend on the observed responses of a whole ",
         "sequence and are therefore undefined on the synthetic grid ",
         "this function builds. Plot one state's own predictor from ",
         "predict(dpar = \"mu2\"), or the occupancies from hmm_probs()",
         call. = FALSE)
  }
  if (isTRUE(surface)) {
    stop("conditional_effects(surface = TRUE) is not implemented: the ",
         "display draws curves with bands, not a fitted surface. Ask ",
         "for the two-variable effect instead, e.g. ",
         "effects = \"x1:x2\", which varies x1 over its range at three ",
         "values of x2 (or at its levels)", call. = FALSE)
  }
  # a per-category effect display is on the CATEGORIES, not the latent
  # scale; naming a dpar explicitly is the way back to the predictor
  categorical <- ce_cats_display(rspec, dpar)
  if (categorical && method == "predict") {
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
      !identical(rspec$family$type, "ordinal")) {
    stop("conditional_effects() has no analytic standard error for the ",
         "category probabilities of family '", rspec$family$family,
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
  dpar_given <- !is.null(dpar)
  dpar <- dpar %||% if ("mu" %in% names(rspec$dpars)) "mu" else
    rspec$primary_dpars[1]
  lp <- find_linpred(x, resp, dpar)
  if (band == "profile") {
    ce_profile_check(x, rspec, lp, dpar_given, categorical)
  }
  # a nonlinear predictor has no delta-method standard error, so its
  # only band is the one that refits
  if (!is.null(lp$nl_body) && band != "boot") {
    stop("conditional_effects() cannot put a ", band, " band on a ",
         "nonlinear predictor: predict() has no standard error for it. ",
         "Use band = \"boot\", which refits instead of differentiating, ",
         "or display one nonlinear parameter with dpar = \"",
         rspec$nlpars[1L], "\"", call. = FALSE)
  }
  base <- data %||% x$frame$data_frame

  vars <- ce_plot_vars(x, rspec, lp, resp)
  vars <- vars[vars %in% names(base)]
  vars <- vars[!vapply(vars, function(v) is.matrix(base[[v]]), TRUE)]
  if (is.null(effects)) {
    effects <- vars
    if (!length(effects)) {
      stop("No plottable predictors found for dpar '", dpar, "'",
           call. = FALSE)
    }
  }
  # one grid per effect: a repeated name would otherwise stack the same
  # grid twice inside its own data frame
  effects <- unique(effects)

  # a data-frame `conditions` defines one condition set per row (brms
  # style); a named list is a single condition set
  cond_sets <- if (is.data.frame(conditions)) {
    stats::setNames(lapply(seq_len(nrow(conditions)), function(r) {
      as.list(conditions[r, , drop = FALSE])
    }), rownames(conditions))
  } else {
    list(conditions)
  }

  z <- stats::qnorm(1 - (1 - prob) / 2)
  # every grid of the call is built before any band is: one bootstrap
  # covers all of them, which is the whole point of doing it here rather
  # than per effect
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
    v1 <- ce_grid_values(base[[ev[1L]]], resolution, ev[1L])
    v2 <- if (length(ev) == 2L) ce_second_values(base[[ev[2L]]])
    n1 <- length(v1)
    n2 <- max(1L, length(v2))
    for (ci in seq_along(cond_sets)) {
      grids[[length(grids) + 1L]] <- list(
        eff = eff, ev = ev, ci = ci, v1 = v1, n1 = n1, n2 = n2,
        n = n1 * n2, cset = cond_sets[[ci]],
        nd = ce_build_nd(base, ev, v1, v2, cond_sets[[ci]], n1 * n2, n2)
      )
    }
  }
  bd <- if (band == "boot") {
    ce_boot_draws(x, grids, categorical, resp, dpar, boot, seed)
  }
  pfail <- c(0L, 0L)

  out <- list()
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
    if (categorical) {
      if (identical(rspec$family$type, "ordinal")) {
        ed <- lp_eta_design(x, lp, nd, FALSE, FALSE)
        ps <- ord_prob_se(x, rspec, lp, ed, nd, FALSE)
      } else {
        # a nominal family has no thresholds, so the ordinal delta
        # method does not apply; under band = "boot" (the only band
        # allowed here) the draws supply the se and the bounds
        P <- predict(x, newdata = nd, type = "response", resp = resp,
                     re.form = NA)
        ps <- list(P = P, se = matrix(NA_real_, nrow(P), ncol(P)))
      }
      cats <- colnames(ps$P)
      df <- do.call(rbind, lapply(seq_along(cats), function(k) {
        d <- nd[ev]
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
        d$cats__ <- factor(cats[k], levels = cats)
        d
      }))
    } else {
      df <- nd[ev]
      if (band == "boot") {
        # the draws are the band AND the standard error here, so the
        # delta method is not asked for: that is what lets a nonlinear
        # predictor, which has no analytic se, reach a band at all
        df$estimate__ <- as.vector(predict(x, newdata = nd,
                                           type = "response", dpar = dpar,
                                           resp = resp, re.form = NA, ...))
        df$se__ <- NA_real_
        df$lower__ <- NA_real_
        df$upper__ <- NA_real_
      } else {
        p <- predict(x, newdata = nd, type = "link", dpar = dpar,
                     resp = resp, re.form = NA, se.fit = TRUE, ...)
        df$estimate__ <- lp$link$linkinv(p$fit)
        df$se__ <- p$se.fit
        df$lower__ <- lp$link$linkinv(p$fit - z * p$se.fit)
        df$upper__ <- lp$link$linkinv(p$fit + z * p$se.fit)
      }
      if (method == "predict") {
        fam <- rspec$family
        if (is.null(fam$sim)) {
          stop("method = 'predict' needs a family with a simulator",
               call. = FALSE)
        }
        dpv <- list()
        for (dnm in names(rspec$dpars)) {
          dpv[[dnm]] <- as.vector(predict(x, newdata = nd, dpar = dnm,
                                          resp = resp, re.form = NA,
                                          type = "response"))
        }
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
      lo <- lp$link$linkinv(pci$lower)
      up <- lp$link$linkinv(pci$upper)
      df$lower__ <- pmin(lo, up)
      df$upper__ <- pmax(lo, up)
      pfail <- pfail + c(pci$fails, pci$tried)
    }
    if (length(cond_sets) > 1L) {
      df$cond__ <- names(cond_sets)[ci] %||% as.character(ci)
    }
    dfs_by_eff[[g$eff]] <- c(dfs_by_eff[[g$eff]], list(df))
  }
  if (pfail[1L] > 0L) {
    warning("band = \"profile\": the likelihood-root search did not ",
            "converge at ", pfail[1L], " of ", pfail[2L],
            " profiled grid point(s); their bounds are NA",
            call. = FALSE)
  }

  for (eff in effects) {
    ev <- strsplit(eff, ":", fixed = TRUE)[[1L]]
    df <- do.call(rbind, dfs_by_eff[[eff]])
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
      attr(df, "points_df") <- pdf_
    }
    out[[eff]] <- df
  }
  out <- structure(out, class = "frmtmb_conditional_effects")
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
                                            ...) {
  ask <- ask %||% (length(x) > 1L && grDevices::dev.interactive())
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
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
      for (cv in unique(df$cond__)) {
        sub <- df[df$cond__ == cv, , drop = FALSE]
        for (a in c("effects", "response", "dpar", "points_df")) {
          attr(sub, a) <- attr(df, a)
        }
        ce_plot_one(sub, cond = cv, points = points)
      }
    } else {
      ce_plot_one(df, points = points)
    }
  }
  invisible(x)
}

#' Draw one conditional-effect panel: the estimate over the varied
#' predictor with its band, lines and a shaded band for a numeric
#' predictor, points and error bars for a discrete one, split by the
#' optional second predictor.
#'
#' @noRd
ce_plot_one <- function(df, cond = NULL, points = FALSE) {
  pts <- if (points) attr(df, "points_df")
  ev <- attr(df, "effects")
  if (!is.null(df[["cats__"]])) {
    # an ordinal display carries one curve per response category, so the
    # category takes the grouping slot; a second predictor then needs a
    # panel of its own rather than a second set of colors
    ylab <- paste0("P(", attr(df, "response"), ")")
    if (!is.null(cond)) ylab <- paste0(ylab, " | ", cond)
    ylim <- range(df$lower__, df$upper__, na.rm = TRUE)
    if (length(ev) == 2L) {
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
  ylim <- range(df$lower__, df$upper__, if (!is.null(pts)) pts$y,
                na.rm = TRUE)
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
    on.exit(grDevices::devAskNewPage(oask))
  }
  if (1L %in% which) {
    # fitted() is a K-column probability matrix on an ordinal fit, and
    # the panel needs one number per row: the expected category index,
    # which is what the residual on the vertical axis was taken against
    rspec <- x$spec$responses[[1L]]
    ordinal <- identical(rspec$family$type, "ordinal")
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
#' @param re.form Passed to [simulate()]; the default `NA` simulates new
#'   random effects.
#' @exportS3Method bayesplot::pp_check
#' @export
pp_check.frmtmb_fit <- function(object, type = "dens_overlay",
                                ndraws = 10, re.form = NA, ...) {
  rspec <- uni_resp(object, "pp_check()")
  y <- object$frame$y[[1L]]
  if (is.matrix(y)) {
    stop("pp_check() on a fit supports vector responses", call. = FALSE)
  }
  sims <- na_unpad(object, simulate(object, nsim = ndraws,
                                    re.form = re.form))
  # ordinal draws come back as ordered factors carrying the response's
  # levels; bayesplot compares them with y, which is the 1..K codes
  yrep <- if (identical(rspec$family$type, "ordinal")) {
    matrix(unlist(lapply(sims, as.integer), use.names = FALSE),
           nrow = nrow(sims))
  } else {
    as.matrix(sims)
  }
  yrep <- t(yrep)
  fun <- get(paste0("ppc_", type), envir = asNamespace("bayesplot"))
  fun(as.numeric(y), yrep, ...)
}
