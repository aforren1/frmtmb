# Two brms post-fit functions on draws: conditional_smooths(), whose
# grid and design are core's (so a draw's curve is the fit's curve at
# that draw), and posterior_average(), brms's model-averaged draws.

#' @rdname sample-conditional-smooths
#' @title Smooth terms of draws
#' @description
#' [frmtmb::conditional_smooths()] on draws from [frm_sample()]: each
#' smooth term over its grid once per draw, summarized the way brms
#' summarizes it. `estimate__` is the median of the drawn curves,
#' `se__` their median absolute deviation, and `lower__` and `upper__`
#' their quantiles at `prob`. The grid, the keys and the columns are
#' core's, so the frame has the layout of the fit method's.
#'
#' `spaghetti = TRUE` adds brms's `"spaghetti"` attribute: the grid once
#' per draw, stacked, with that draw's curve in `estimate__` and the draw
#' in `sample__`. brms adds it only for a term whose first covariate is
#' numeric and that is not drawn as a surface, and so does this method.
#'
#' Draws from `frm_sample(laplace = TRUE)` are refused: they hold no
#' penalized smooth coefficients.
#' @param x A `frmtmb_draws` object.
#' @inheritParams frmtmb::conditional_smooths
#' @return As [frmtmb::conditional_smooths()].
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE) &&
#'     !frmtmb.sample:::tmbstan_build_broken()) {
#'   set.seed(1)
#'   dd <- data.frame(x = runif(80))
#'   dd$y <- sin(2 * pi * dd$x) + rnorm(80, 0, 0.3)
#'   ds <- frm_sample(bf(y ~ s(x)), family = gaussian(), data = dd,
#'                    chains = 1, iter = 400, seed = 1, refresh = 0)
#'   cs <- conditional_smooths(ds, spaghetti = TRUE, ndraws = 20)
#'   nrow(attr(cs[["mu: s(x)"]], "spaghetti"))
#' }
#' }
#' @exportS3Method brms::conditional_smooths
#' @export
conditional_smooths.frmtmb_draws <- function(x, smooths = NULL,
                                             int_conditions = NULL,
                                             prob = 0.95, spaghetti = FALSE,
                                             surface = TRUE,
                                             resolution = 100, too_far = 0,
                                             ndraws = NULL, draw_ids = NULL,
                                             nsamples = NULL, subset = NULL,
                                             probs = NULL, ...) {
  frm_check_dots(...)
  check_flag(spaghetti, "spaghetti")
  band_p <- cs_probs(prob, probs)
  # brms's deprecated spellings, with brms's rule that the new name wins
  ndraws <- ndraws %||% nsamples
  draw_ids <- draw_ids %||% subset
  fit <- draws_base_fit(x)
  if (length(fit$frame[["re_blocks"]]) && draws_is_laplace(x)) {
    frm_stop("conditional_smooths() on draws from frm_sample(laplace = ",
             "TRUE) has no smooth coefficients to draw: the penalized ",
             "coefficients were integrated out. Resample without ",
             "laplace = TRUE, or call conditional_smooths() on the fit",
             call. = FALSE)
  }
  terms <- cs_build(fit, smooths, int_conditions, surface, resolution,
                    too_far)
  if (!length(terms)) {
    frm_stop("No valid smooth terms found in the model.", call. = FALSE)
  }
  rows <- draws_subsample(x, ndraws, draw_ids)
  idx <- draws_par_index(x$fit)
  # one parameter vector per draw, shared by every term
  coefs <- lapply(rows, function(i) {
    fi <- draws_fit_at(x, i, idx)
    lapply(terms, function(tm) cs_coef(fi, tm))
  })
  out <- lapply(names(terms), function(k) {
    tm <- terms[[k]]
    eta <- do.call(rbind, lapply(coefs, function(cf) {
      as.vector(tm$A %*% cf[[k]])
    }))
    spag <- if (spaghetti && tm$spaghetti_ok) {
      ce_spaghetti(tm$cond_data, eta, tm$effects)
    }
    cs_frame(tm, apply(eta, 2, stats::median), apply(eta, 2, stats::mad),
             ce_pctl(eta, band_p[1L]), ce_pctl(eta, band_p[2L]), spag)
  })
  names(out) <- names(terms)
  cs_finalize(out)
}

#' Posterior draws averaged over models
#'
#' brms's `posterior_average()` on draws from [frm_sample()]: a sample
#' of posterior draws taken from several models in proportion to their
#' weights. The number of draws taken from each model is the weight
#' times `ndraws`, rounded by the largest remainder, and the draws of a
#' model are taken without replacement. The result equals brms's on the
#' same draws, the same weights and the same seed.
#'
#' @section Weights:
#' `weights` is a numeric vector with one weight per model, or the name
#' of a method that computes them:
#' * `"stacking"` (the default) and `"pseudobma"`:
#'   [loo::loo_model_weights()] over the [loo()] of each model.
#' * `"loo"` and `"waic"`: Akaike-type weights,
#'   `exp(-delta / 2)` normalized, over each model's LOOIC or WAIC.
#' * `"kfold"` and `"bma"` are refused: frmtmb.sample has no `kfold()`
#'   and no marginal likelihood (see [kfold()] and [post_prob()]).
#' @param x,... Draws objects from [frm_sample()], one per model. Every
#'   argument in `...` must be a model; name one to label it.
#' @param variable Names of the variables to average. `NULL` takes the
#'   variables that every model has.
#' @param pars brms's deprecated spelling of `variable`.
#' @param weights Numeric weights, one per model, or a method name. See
#'   the weights section.
#' @param ndraws Total number of draws. `NULL` takes the number of draws
#'   of the first model.
#' @param nsamples brms's deprecated spelling of `ndraws`.
#' @param missing What a model that lacks a variable contributes for it:
#'   one number for every variable, or a named list with one number per
#'   missing variable. `NULL` (the default) refuses a variable that some
#'   model lacks.
#' @param model_names Labels for the models. The default deparses the
#'   arguments, as brms does.
#' @param control A list of arguments passed on to the weight method.
#' @param seed Seed for the draw selection.
#' @return A data frame with one row per draw and one column per
#'   variable, the draws of the first model first. The attribute
#'   `weights` holds the normalized weights and the attribute `ndraws`
#'   the number of draws taken from each model, both named by model.
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE) &&
#'     !frmtmb.sample:::tmbstan_build_broken()) {
#'   set.seed(1)
#'   dd <- data.frame(x = rnorm(60))
#'   dd$y <- rnorm(60, 1 + 0.5 * dd$x)
#'   d1 <- frm_sample(bf(y ~ x), family = gaussian(), data = dd,
#'                    chains = 1, iter = 400, seed = 1, refresh = 0)
#'   d2 <- frm_sample(bf(y ~ 1), family = gaussian(), data = dd,
#'                    chains = 1, iter = 400, seed = 1, refresh = 0)
#'   pa <- posterior_average(d1, d2, variable = "sigma",
#'                           weights = c(0.7, 0.3), seed = 1)
#'   attr(pa, "ndraws")
#' }
#' }
#' @export
posterior_average <- function(x, ...) UseMethod("posterior_average")

#' @rdname posterior_average
#' @exportS3Method brms::posterior_average
#' @export
posterior_average.frmtmb_draws <- function(x, ..., variable = NULL,
                                           pars = NULL,
                                           weights = "stacking",
                                           ndraws = NULL, nsamples = NULL,
                                           missing = NULL,
                                           model_names = NULL,
                                           control = list(), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  variable <- variable %||% pars
  ndraws <- ndraws %||% nsamples
  models <- c(list(x), list(...))
  labels <- vapply(as.list(substitute(list(x, ...)))[-1L],
                   function(e) paste(deparse(e), collapse = ""), "")
  given <- names(models) %||% rep("", length(models))
  labels[nzchar(given)] <- given[nzchar(given)]
  bad <- which(!vapply(models, inherits, NA, "frmtmb_draws"))
  if (length(bad)) {
    frm_stop("posterior_average() averages draws objects, and argument ",
             bad[1L], " is a ", paste(class(models[[bad[1L]]]),
                                      collapse = "/"),
             ". Draw from each model with frm_sample() first",
             call. = FALSE)
  }
  if (!is.null(model_names)) {
    if (length(model_names) != length(models) ||
        anyDuplicated(model_names)) {
      frm_stop("`model_names` must give one distinct name per model",
               call. = FALSE)
    }
    labels <- as.character(model_names)
  }
  # a repeated label stays repeated, as brms leaves it: the attributes
  # are read by position, and brms's own output names them that way
  names(models) <- labels
  vars_list <- lapply(models, variables)
  all_vars <- unique(unlist(vars_list))
  if (is.null(missing)) {
    common <- Reduce(intersect, vars_list)
    variable <- as.character(variable %||% setdiff(common, "lp__"))
    inv <- setdiff(variable, common)
    if (length(inv)) {
      frm_stop("Parameters ", paste(inv, collapse = ", "), " cannot be ",
               "found in all of the models. Consider using argument ",
               "'missing'.", call. = FALSE)
    }
  } else {
    variable <- as.character(variable %||% setdiff(all_vars, "lp__"))
    inv <- setdiff(variable, all_vars)
    if (length(inv)) {
      frm_stop("Parameters ", paste(inv, collapse = ", "), " cannot be ",
               "found in any of the models.", call. = FALSE)
    }
    # brms's as_one_numeric(allow_na = TRUE): one number, or one NA
    one_num <- function(v) {
      if (length(v) != 1L || !(is.numeric(v) || is.na(v))) {
        frm_stop("`missing` must hold single numbers", call. = FALSE)
      }
      as.numeric(v)
    }
    if (is.list(missing)) {
      need <- unique(unlist(lapply(vars_list, function(v) {
        setdiff(variable, v)
      })))
      inv <- setdiff(need, names(missing))
      if (length(inv)) {
        frm_stop("Argument 'missing' has no value for parameters ",
                 paste(inv, collapse = ", "), ".", call. = FALSE)
      }
      missing <- lapply(missing, one_num)
    } else {
      missing <- stats::setNames(rep(list(one_num(missing)),
                                     length(variable)), variable)
    }
  }
  ndraws <- ndraws %||% ndraws(models[[1L]])
  if (!is.numeric(ndraws) || length(ndraws) != 1L || is.na(ndraws) ||
      ndraws < 0 || ndraws != round(ndraws)) {
    frm_stop("`ndraws` must be one whole number of at least 0",
             call. = FALSE)
  }
  weights <- pa_weights(weights, models, control)
  nd <- pa_round_largest_remainder(weights * ndraws)
  names(weights) <- names(nd) <- labels
  out <- stats::setNames(vector("list", length(models)), labels)
  for (i in seq_along(models)) {
    if (nd[i] <= 0) next
    avail <- ndraws(models[[i]])
    if (nd[i] > avail) {
      frm_stop("posterior_average() needs ", nd[i], " draws from model '",
               labels[i], "', which has ", avail, ". Lower `ndraws`",
               call. = FALSE)
    }
    draw <- sort(sample(seq_len(avail), nd[i]))
    found <- intersect(variable, vars_list[[i]])
    out[[i]] <- if (length(found)) {
      as.data.frame(models[[i]], variable = found, draw = draw)
    } else {
      as.data.frame(matrix(numeric(0), nrow = nd[i], ncol = 0))
    }
    if (!is.null(missing)) {
      miss <- setdiff(variable, names(out[[i]]))
      if (length(miss)) out[[i]][miss] <- missing[miss]
    }
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  attr(out, "weights") <- weights
  attr(out, "ndraws") <- nd
  out
}

#' Model weights for posterior_average(), normalized: numeric weights as
#' given, or brms's model_weights() methods over draws.
#'
#' @noRd
pa_weights <- function(weights, models, control = list()) {
  if (is.numeric(weights)) {
    if (length(weights) != length(models)) {
      frm_stop("If numeric, 'weights' must have the same length as the ",
               "number of models.", call. = FALSE)
    }
    if (anyNA(weights) || any(weights < 0)) {
      frm_stop("If numeric, 'weights' must be positive.", call. = FALSE)
    }
    if (sum(weights) <= 0) {
      frm_stop("If numeric, 'weights' must not all be zero",
               call. = FALSE)
    }
    return(weights / sum(weights))
  }
  if (!is.character(weights) || length(weights) != 1L || is.na(weights)) {
    frm_stop("`weights` must be numeric or one method name, not a ",
             paste(class(weights), collapse = "/"), " of length ",
             length(weights), call. = FALSE)
  }
  method <- tolower(weights)
  if (method == "loo2") method <- "stacking"
  if (method == "marglik") method <- "bma"
  method <- frm_match_arg(method, c("loo", "waic", "kfold", "stacking",
                                    "pseudobma", "bma"))
  if (method %in% c("kfold", "bma")) {
    frm_stop("weights = \"", method, "\" needs ",
             if (method == "kfold") "kfold()" else "a marginal likelihood",
             ", which draws from frm_sample() do not have. Use ",
             "\"stacking\", \"pseudobma\", \"loo\" or \"waic\", or give ",
             "numeric weights", call. = FALSE)
  }
  if (!requireNamespace("loo", quietly = TRUE)) {
    frm_stop("weights = \"", method, "\" needs the 'loo' package; install ",
             "it, or give numeric weights", call. = FALSE)
  }
  w <- if (method %in% c("loo", "waic")) {
    ics <- vapply(models, function(m) {
      ic <- if (method == "loo") loo(m) else waic(m)
      ic$estimates[3L, 1L]
    }, 1)
    exp(-(ics - min(ics)) / 2)
  } else {
    crit <- lapply(models, loo)
    as.numeric(do.call(loo::loo_model_weights,
                       c(list(crit), list(method = method), control)))
  }
  w / sum(w)
}

#' brms's round_largest_remainder(): whole numbers that sum to the
#' rounded total and differ from `x` by less than one each.
#'
#' @noRd
pa_round_largest_remainder <- function(x) {
  x <- as.numeric(x)
  total <- round(sum(x))
  out <- floor(x)
  j <- order(x - out, decreasing = TRUE)
  i <- seq_len(total - floor(sum(out)))
  out[j[i]] <- out[j[i]] + 1
  out
}
