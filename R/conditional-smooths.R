#' Display the smooth terms of a model
#'
#' Draws each smooth term (`s()`, `t2()`) of a model on its own: the
#' term's contribution to its linear predictor over a grid of its
#' covariates, without the intercept or any other term. This is brms's
#' `conditional_smooths()`, with brms's grid, keys and column layout.
#' [conditional_effects()] draws the expected response instead, where a
#' smooth is one term among the others.
#'
#' Every smooth term is drawn, a term indexed by a grouping factor (a
#' factor-smooth `s(x, g, bs = "fs")` or a `bs = "re"` smooth) included.
#' This is the rule `re_formula = NA` has everywhere in frmtmb: it drops
#' the `(... | g)` group-level terms and keeps every smooth, as brms
#' does.
#'
#' @section The grid:
#' The grid follows brms. A term's covariates vary and its `by`
#' variable is a second covariate or a condition:
#' * A numeric covariate takes `resolution` equally spaced values over
#'   its observed range. When it is the second covariate and `surface`
#'   is `FALSE`, it takes its mean and its mean plus or minus one
#'   standard deviation instead, as a moderator.
#' * A factor covariate takes each of its levels.
#' * A term of one covariate and a `by` variable takes the `by`
#'   variable as its second covariate. A term of more than two
#'   covariates varies the first two and holds the others at their mean
#'   and mean plus or minus one standard deviation (numeric) or at each
#'   level (factor), labeled in `cond__`.
#' * `int_conditions` gives the values of a covariate by name, as values
#'   or as a function of the observed column.
#' * With two numeric covariates and `surface = TRUE` (the default) the
#'   grid is `resolution` by `resolution` and the term is drawn as a
#'   surface. `too_far` then removes the grid points farther than that
#'   from every observation, on the unit square, with
#'   [mgcv::exclude.too.far()].
#'
#' @section The band:
#' On a maximum-likelihood fit `estimate__` is the term at the estimate
#' and `se__` is its delta-method standard error, over the covariance of
#' the term's coefficients: the null-space coefficients (fixed effects)
#' and the penalized ones (random effects, conditional on their modes
#' as a smooth's are). `lower__` and `upper__` are `estimate__` plus the
#' normal quantiles of `prob` times `se__`. The term is on the linear
#' predictor scale, so the band is symmetric. On draws from
#' `frmtmb.sample::frm_sample()` the columns are brms's: the median, the
#' median absolute deviation and the quantiles of the drawn curves.
#'
#' @param x A `frmtmb_fit`, or draws from `frmtmb.sample::frm_sample()`.
#' @param smooths Character vector of the terms to draw, as written in
#'   the formula (whitespace does not matter): `"s(x)"`,
#'   `"s(z, by = f)"`. `NULL` (the default) draws every smooth term.
#' @param int_conditions Named list of the values a covariate takes, as
#'   in [conditional_effects()].
#' @param prob Coverage of the band.
#' @param spaghetti `TRUE` adds one curve per posterior draw as the
#'   attribute `"spaghetti"`, for a term whose first covariate is
#'   numeric and that is not drawn as a surface. A maximum-likelihood
#'   fit has no draws and refuses it.
#' @param surface Whether a term of two numeric covariates is drawn as a
#'   surface over both. `FALSE` draws the second covariate as a
#'   moderator at three values.
#' @param resolution Number of grid points for a varied numeric
#'   covariate.
#' @param too_far For a surface: remove grid points farther than this
#'   from the data, on the unit square. `0` (the default) keeps all.
#' @param ndraws,draw_ids For draws: how many evenly spaced draws, or
#'   which draws, to use. A maximum-likelihood fit refuses both.
#' @param nsamples,subset brms's deprecated spellings of `ndraws` and
#'   `draw_ids`.
#' @param probs brms's deprecated spelling of the two band quantiles,
#'   which replaces `prob` with a warning.
#' @param ... Refused by name.
#' @return A named list of data frames, one per term and linear
#'   predictor, keyed as brms keys them: the distributional parameter,
#'   the response of a multivariate model and the nonlinear parameter,
#'   joined by `_`, then `": "` and the term, for example
#'   `"mu: s(x)"`, `"sigma: s(z)"` or `"mu_a: s(x)"`. A frame has the
#'   term's covariates and `by` variables, `effect1__` and (for two
#'   covariates) `effect2__`, `cond__`, then `estimate__`, `se__`,
#'   `lower__` and `upper__`. Its attributes are `response` (the key),
#'   `effects`, `surface`, `spaghetti` and `points` (the observed
#'   covariate values). The list has class
#'   `"frmtmb_conditional_effects"` and the attribute
#'   `smooths_only = TRUE`, so `plot()` draws it.
#' @seealso [conditional_effects()] for the expected response,
#'   [make_conditions()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(150), z = runif(150))
#' dd$y <- sin(2 * pi * dd$x) + dd$z + rnorm(150, 0, 0.3)
#' fit <- frm(bf(y ~ z + s(x)), family = gaussian(), data = dd)
#' cs <- conditional_smooths(fit)
#' head(cs[["mu: s(x)"]])
#' plot(cs, ask = FALSE)
#' @export
conditional_smooths <- function(x, ...) UseMethod("conditional_smooths")

#' @rdname conditional_smooths
#' @exportS3Method brms::conditional_smooths
#' @export
conditional_smooths.frmtmb_fit <- function(x, smooths = NULL,
                                           int_conditions = NULL,
                                           prob = 0.95, spaghetti = FALSE,
                                           surface = TRUE,
                                           resolution = 100, too_far = 0,
                                           ndraws = NULL, draw_ids = NULL,
                                           nsamples = NULL, subset = NULL,
                                           probs = NULL, ...) {
  frm_check_dots(...)
  require_fitted(x, "conditional_smooths()")
  check_flag(spaghetti, "spaghetti")
  band_p <- cs_probs(prob, probs)
  # a maximum-likelihood fit has one parameter vector, so an argument
  # that picks or thins draws has nothing to act on; ignoring it would
  # return a band the caller did not ask for
  given <- c(ndraws = !is.null(ndraws), draw_ids = !is.null(draw_ids),
             nsamples = !is.null(nsamples), subset = !is.null(subset),
             spaghetti = isTRUE(spaghetti))
  if (any(given)) {
    frm_stop("conditional_smooths() on a maximum-likelihood fit has no ",
             "draws, so `", names(given)[given][1L], "` has nothing to act ",
             "on. Draw from the posterior with frmtmb.sample::frm_sample() ",
             "and call conditional_smooths() on the draws", call. = FALSE)
  }
  terms <- cs_build(x, smooths, int_conditions, surface, resolution,
                    too_far)
  if (!length(terms)) {
    frm_stop("No valid smooth terms found in the model.", call. = FALSE)
  }
  jc <- get_joint_cov(x)
  out <- lapply(terms, function(tm) {
    pos <- cs_joint_pos(x, tm, jc)
    est <- as.vector(tm$A %*% cs_coef(x, tm))
    keep <- !is.na(pos)
    Ak <- tm$A[, keep, drop = FALSE]
    V <- jc$V[pos[keep], pos[keep], drop = FALSE]
    se <- sqrt(pmax(rowSums((Ak %*% V) * Ak), 0))
    cs_frame(tm, est, se, est + stats::qnorm(band_p[1L]) * se,
             est + stats::qnorm(band_p[2L]) * se)
  })
  cs_finalize(out)
}

#' The two band quantiles from brms's `prob` and its deprecated
#' `probs`.
#'
#' @noRd
cs_probs <- function(prob, probs = NULL) {
  check_probability(prob, "prob")
  if (is.null(probs)) return(c((1 - prob) / 2, 1 - (1 - prob) / 2))
  frm_warning("Argument 'probs' is deprecated. Please use 'prob' instead.",
              call. = FALSE)
  if (!is.numeric(probs) || length(probs) != 2L || anyNA(probs) ||
      any(probs < 0 | probs > 1)) {
    frm_stop("`probs` must be two probabilities, not ", arg_desc(probs),
             call. = FALSE)
  }
  sort(as.numeric(probs))
}

#' Every smooth term of a fit as one display grid each, in brms's order:
#' response by response, the distributional parameters, then the
#' nonlinear ones, and the terms of each in formula order.
#'
#' A term is a written `s()` or `t2()` call. A factor `by` gives it one
#' basis per level, all drawn together, which is why the grouping is by
#' the written term and not by mgcv's label.
#'
#' Each grid carries its design `A` over the term's coefficients, which
#' does not depend on their values: the fit method takes one product
#' with the estimate, a sampler one per draw ([cs_coef()]).
#'
#' @noRd
cs_build <- function(fit, smooths = NULL, int_conditions = NULL,
                     surface = TRUE, resolution = 100, too_far = 0) {
  check_flag(surface, "surface")
  check_count(resolution, "resolution", min = 1L)
  if (!is.numeric(too_far) || length(too_far) != 1L || is.na(too_far) ||
      too_far < 0) {
    frm_stop("`too_far` must be one number of at least 0, not ",
             arg_desc(too_far), call. = FALSE)
  }
  if (!is.null(int_conditions)) {
    check_named_list(int_conditions, "int_conditions",
                     "int_conditions = list(z = c(-1, 0, 1))")
  }
  if (!is.null(smooths) && !is.character(smooths)) {
    frm_stop("`smooths` must be a character vector of smooth terms, ",
             "such as \"s(x)\", not ", arg_desc(smooths), call. = FALSE)
  }
  smooths <- gsub("[ \t\r\n]+", "", as.character(smooths))
  mf <- fit$frame[["data_frame"]]
  multi <- length(fit$spec$responses) > 1L
  out <- list()
  for (resp in names(fit$spec$responses)) {
    rspec <- fit$spec$responses[[resp]]
    nlp <- rspec$nlpars %||% character(0)
    for (par in c(setdiff(names(rspec$dpars), nlp), nlp)) {
      lp <- fit$frame[["linpreds"]][[linpred_key(rspec$resp_name, par)]]
      sis <- lp[["smooths"]] %||% list()
      if (!length(sis)) next
      prefix <- if (par %in% nlp) c("mu", if (multi) resp, par) else
        c(par, if (multi) resp)
      prefix <- paste(prefix, collapse = "_")
      tlab <- vapply(sis, function(si) si[["term"]] %||% si$label, "")
      for (term in unique(tlab)) {
        if (length(smooths) && !term %in% smooths) next
        tm <- cs_term_grid(fit, lp, sis[tlab == term], term, mf,
                           int_conditions, surface, resolution, too_far)
        tm$key <- paste0(prefix, ": ", term)
        out[[tm$key]] <- tm
      }
    }
  }
  out
}

#' One term's grid and design. The grid rules are brms's
#' `conditional_smooths.btl()`, line for line, because a ported script
#' indexes the frame by row.
#'
#' @noRd
cs_term_grid <- function(fit, lp, sis, term, mf, int_conditions, surface,
                         resolution, too_far) {
  sm <- sis[[1L]]$sm
  covars <- unique(unlist(lapply(sm$term, function(v) {
    all.vars(str2lang(v))
  })))
  byvars <- if (!identical(sm$by, "NA")) all.vars(str2lang(sm$by))
  if (length(covars) > 2L) {
    byvars <- c(covars[-(1:2)], byvars)
    covars <- covars[1:2]
  } else if (length(covars) == 1L && length(byvars)) {
    covars <- c(covars, byvars[1L])
    byvars <- byvars[-1L]
  }
  vars <- c(covars, byvars)
  miss <- setdiff(vars, names(mf))
  if (length(miss)) {
    frm_stop("conditional_smooths(): the term ", term, " reads ",
             paste0("`", miss, "`", collapse = ", "), ", which the model ",
             "data does not store", call. = FALSE)
  }
  int_value <- function(v) {
    ic <- int_conditions[[v]]
    if (is.function(ic)) ic(mf[[v]]) else ic
  }
  mod_values <- function(v) {
    if (is.numeric(mf[[v]])) {
      (-1:1) * stats::sd(mf[[v]], na.rm = TRUE) + mean(mf[[v]], na.rm = TRUE)
    } else {
      levels(factor(mf[[v]]))
    }
  }
  values <- stats::setNames(vector("list", length(vars)), vars)
  is_numeric <- stats::setNames(rep(FALSE, length(covars)), covars)
  for (cv in covars) {
    is_numeric[cv] <- is.numeric(mf[[cv]])
    values[[cv]] <- if (cv %in% names(int_conditions)) {
      int_value(cv)
    } else if (is_numeric[cv]) {
      if (!surface && identical(cv, covars[2L])) {
        mod_values(cv)
      } else {
        seq(min(mf[[cv]]), max(mf[[cv]]), length.out = resolution)
      }
    } else {
      levels(factor(mf[[cv]]))
    }
  }
  for (cv in byvars) {
    values[[cv]] <- if (cv %in% names(int_conditions)) {
      int_value(cv)
    } else {
      mod_values(cv)
    }
  }
  newdata <- expand.grid(values)
  show_surface <- surface && length(covars) == 2L && all(is_numeric)
  if (show_surface && too_far > 0) {
    far <- mgcv::exclude.too.far(g1 = newdata[[covars[1L]]],
                                 g2 = newdata[[covars[2L]]],
                                 d1 = mf[[covars[1L]]],
                                 d2 = mf[[covars[2L]]], dist = too_far)
    newdata <- newdata[!far, , drop = FALSE]
  }
  # a factor column of the grid holds the level names expand.grid()
  # made; the basis needs the fitted levels, in the fitted order
  pdata <- newdata
  for (v in names(pdata)) {
    if (is.factor(mf[[v]]) || is.character(mf[[v]])) {
      pdata[[v]] <- factor(as.character(pdata[[v]]),
                           levels = levels(factor(mf[[v]])))
    }
  }
  A <- NULL
  parts <- list()
  for (si in sis) {
    smooth_newdata_check(si, pdata, FALSE)
    M <- smooth_basis_at(si, pdata)
    A <- cbind(A, M)
    parts[[length(parts) + 1L]] <- si
  }
  effects <- stats::na.omit(covars[1:2])
  effects <- as.character(effects)
  cond_data <- newdata[, vars, drop = FALSE]
  for (i in seq_along(effects)) {
    cond_data[[paste0("effect", i, "__")]] <- cond_data[[effects[i]]]
  }
  if (isTRUE(is_numeric[2L]) && !surface) {
    mde2 <- round(cond_data[[effects[2L]]], 2)
    levels2 <- sort(unique(mde2), TRUE)
    cond_data$effect2__ <- factor(mde2, levels = levels2)
    labels2 <- names(int_conditions[[effects[2L]]])
    if (length(labels2) == length(levels2)) {
      levels(cond_data$effect2__) <- labels2
    }
  }
  cond_data$cond__ <- if (length(byvars)) {
    brms_rows2labels(cond_data[, byvars, drop = FALSE])
  } else {
    factor(1)
  }
  points <- mf[, vars, drop = FALSE]
  for (i in seq_along(covars)) {
    points[[paste0("effect", i, "__")]] <- points[[covars[i]]]
  }
  list(term = term, lp = lp, sis = parts, A = A, cond_data = cond_data,
       effects = effects, covars = covars,
       spaghetti_ok = !show_surface && isTRUE(is_numeric[1L]),
       surface = show_surface, points = points)
}

#' The term's coefficients at the parameter values `fit` carries, in the
#' column order of the grid's `A`: per basis, the penalized blocks, then
#' the null-space columns. A sampler passes one fit per draw.
#'
#' @noRd
cs_coef <- function(fit, tm) {
  est <- fit$estimates
  cvec <- coef_b(fit)
  lp <- tm$lp
  unlist(lapply(tm$sis, function(si) {
    c(unlist(lapply(si$block_ids, function(k) {
      cvec[fit$frame[["re_blocks"]][[k]][["c_idx"]]]
    })), est[[lp[["par"]]]][lp[["idx"]][si$xf_idx]])
  }), use.names = FALSE)
}

#' Rows of the joint covariance for the columns of the grid's `A`, NA
#' for a coefficient the fit holds fixed. The fixed-effect rows follow
#' `lp_delta_A()`'s reading of `betad`, whose held entries have no row.
#'
#' @noRd
cs_joint_pos <- function(fit, tm, jc) {
  rn <- jc$names
  lp <- tm$lp
  b_pos <- which(rn == "b")
  if (!length(b_pos) && length(fit$frame[["re_blocks"]])) {
    warn_modes_conditional_se()
  }
  xpos <- if (lp[["par"]] == "beta") {
    which(rn == "beta")[lp[["idx"]]]
  } else {
    tpl_len <- length(fit$frame[["par_template"]][["betad"]])
    rank <- match(lp[["idx"]], setdiff(seq_len(tpl_len),
                                       fit$frame[["betad_fixed_idx"]]))
    which(rn == "betad")[rank]
  }
  unlist(lapply(tm$sis, function(si) {
    c(unlist(lapply(si$block_ids, function(k) {
      bi <- fit$frame[["re_blocks"]][[k]][["b_idx"]]
      if (length(b_pos)) b_pos[bi] else rep(NA_integer_, length(bi))
    })), xpos[si$xf_idx])
  }), use.names = FALSE)
}

#' One term's frame in brms's layout, from the summary columns a method
#' computed.
#'
#' @noRd
cs_frame <- function(tm, estimate, se, lower, upper, spaghetti = NULL) {
  out <- tm$cond_data
  out$estimate__ <- estimate
  out$se__ <- se
  out$lower__ <- lower
  out$upper__ <- upper
  attr(out, "response") <- tm$key
  attr(out, "effects") <- tm$effects
  attr(out, "surface") <- tm$surface
  attr(out, "spaghetti") <- spaghetti
  attr(out, "points") <- tm$points
  out
}

#' brms's spaghetti frame: the grid once per draw, stacked, with each
#' draw's curve in `estimate__` and its id in `sample__`. `m` has one
#' row per draw and one column per row of `grid`.
#'
#' @noRd
ce_spaghetti <- function(grid, m, effects) {
  sample <- rep(seq_len(nrow(m)), each = ncol(m))
  if (length(effects) == 2L) {
    sample <- paste0(sample, "_", grid[[effects[2L]]])
  }
  spag <- data.frame(estimate__ = as.numeric(t(m)),
                     sample__ = factor(sample))
  rep_grid <- grid[rep(seq_len(nrow(grid)), nrow(m)), , drop = FALSE]
  rownames(rep_grid) <- NULL
  cbind(rep_grid, spag)
}

#' The list the methods return.
#'
#' @noRd
cs_finalize <- function(frames) {
  structure(frames, class = "frmtmb_conditional_effects",
            smooths_only = TRUE)
}
