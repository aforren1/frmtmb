# Conditional-effects displays, diagnostic plot method, pp_check.

# Addition-term values for the conditional-effects grid. A grid row is
# an artificial observation, so an aterm that changes the predictive
# distribution (trials, se, truncation bounds) must not be taken at a
# reference value: the mean number of trials is rarely a whole number,
# and a mean truncation bound is nobody's bound. Those terms are read
# only from variables the user pinned in `conditions`; literal bounds
# apply as written. Everything else (vint/vreal payloads a custom
# family needs) is evaluated against the grid when it can be.
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

# Reference value a predictor is held at when it is not varied.
ce_ref_value <- function(col) {
  if (is.matrix(col)) {
    matrix(colMeans(col), 1, ncol(col))
  } else if (is.factor(col)) {
    factor(levels(col)[1L], levels = levels(col))
  } else if (is.numeric(col)) {
    mean(col)
  } else if (is.logical(col)) {
    FALSE
  } else {
    sort(unique(col))[1L]
  }
}

# Grid of values for the varied (first) predictor.
ce_grid_values <- function(col, resolution) {
  if (is.factor(col)) {
    factor(levels(col), levels = levels(col))
  } else if (is.numeric(col)) {
    seq(min(col), max(col), length.out = resolution)
  } else if (is.logical(col)) {
    c(FALSE, TRUE)
  } else {
    sort(unique(col))
  }
}

# Values for the second predictor of an "x:z" effect.
ce_second_values <- function(col) {
  if (is.numeric(col) && !is.matrix(col)) {
    signif(mean(col) + c(-1, 0, 1) * stats::sd(col), 3)
  } else {
    ce_grid_values(col, resolution = 0)
  }
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
#' @param data The original model data. Only needed when the model frame
#'   does not store a raw variable (e.g. a variable used only inside
#'   `poly()`).
#' @param ... Passed to [predict.frmtmb_fit()].
#' @return A named list of data frames (one per effect) with the varied
#'   variable(s) plus `estimate__`, `se__` (link scale), `lower__`, and
#'   `upper__`; printing it draws the plots.
#' @examples
#' set.seed(5)
#' dd <- data.frame(x = rnorm(120), f = factor(rep(c("a", "b"), 60)))
#' dd$y <- rnorm(120, 1 + 0.5 * dd$x + (dd$f == "b"), 1)
#' fit <- frm(bf(y ~ x * f) + gaussian(), data = dd)
#' ce <- conditional_effects(fit, effects = c("x", "x:f"))
#' plot(ce, ask = FALSE)
#' # prediction intervals instead of epred bands
#' ce_p <- conditional_effects(fit, effects = "x", method = "predict")
#' @export
conditional_effects <- function(x, ...) UseMethod("conditional_effects")

#' @rdname conditional_effects
#' @export
conditional_effects.frmtmb_fit <- function(x, effects = NULL, resp = NULL,
                                           dpar = NULL, resolution = 100,
                                           prob = 0.95,
                                           method = c("epred", "predict"),
                                           ndraws = 400,
                                           conditions = list(),
                                           data = NULL, ...) {
  method <- match.arg(method)
  resp <- resp %||% names(x$spec$responses)[1L]
  rspec <- x$spec$responses[[resp]]
  dpar <- dpar %||% if ("mu" %in% names(rspec$dpars)) "mu" else
    rspec$primary_dpars[1]
  lp <- find_linpred(x, resp, dpar)
  base <- data %||% x$frame$data_frame

  fixed_vars <- all.vars(stats::delete.response(lp$terms))
  sm_vars <- unlist(lapply(lp$smooths, function(si) si$sm$term))
  vars <- unique(c(fixed_vars, sm_vars))
  vars <- vars[vars %in% names(base)]
  vars <- vars[!vapply(vars, function(v) is.matrix(base[[v]]), TRUE)]
  if (is.null(effects)) {
    effects <- vars
    if (!length(effects)) {
      stop("No plottable predictors found for dpar '", dpar, "'",
           call. = FALSE)
    }
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

  z <- stats::qnorm(1 - (1 - prob) / 2)
  out <- list()
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
    v1 <- ce_grid_values(base[[ev[1L]]], resolution)
    v2 <- if (length(ev) == 2L) ce_second_values(base[[ev[2L]]])
    n1 <- length(v1)
    n2 <- max(1L, length(v2))
    n <- n1 * n2

    dfs <- list()
    for (ci in seq_along(cond_sets)) {
      cset <- cond_sets[[ci]]
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
      if (length(ev) == 2L) nd[[ev[2L]]] <- rep(v2, each = n1)
      nd$.ce_row <- NULL

      p <- predict(x, newdata = nd, type = "link", dpar = dpar,
                   resp = resp, re.form = NA, se.fit = TRUE, ...)
      df <- nd[ev]
      df$estimate__ <- lp$link$linkinv(p$fit)
      df$se__ <- p$se.fit
      df$lower__ <- lp$link$linkinv(p$fit - z * p$se.fit)
      df$upper__ <- lp$link$linkinv(p$fit + z * p$se.fit)
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
      if (length(cond_sets) > 1L) {
        df$cond__ <- names(cond_sets)[ci] %||% as.character(ci)
      }
      dfs[[ci]] <- df
    }
    df <- do.call(rbind, dfs)
    attr(df, "effects") <- ev
    attr(df, "response") <- resp
    attr(df, "dpar") <- dpar
    out[[eff]] <- df
  }
  structure(out, class = "frmtmb_conditional_effects")
}

#' @export
print.frmtmb_conditional_effects <- function(x, ...) {
  plot(x, ...)
  invisible(x)
}

#' @export
plot.frmtmb_conditional_effects <- function(x, ask = NULL, ...) {
  ask <- ask %||% (length(x) > 1L && grDevices::dev.interactive())
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }
  for (nm in names(x)) {
    df <- x[[nm]]
    if (!is.null(df$cond__) && length(unique(df$cond__)) > 1L) {
      for (cv in unique(df$cond__)) {
        sub <- df[df$cond__ == cv, , drop = FALSE]
        for (a in c("effects", "response", "dpar")) {
          attr(sub, a) <- attr(df, a)
        }
        ce_plot_one(sub, cond = cv)
      }
    } else {
      ce_plot_one(df)
    }
  }
  invisible(x)
}

ce_plot_one <- function(df, cond = NULL) {
  ev <- attr(df, "effects")
  ylab <- paste0(attr(df, "response"), " (", attr(df, "dpar"), ")")
  if (!is.null(cond)) ylab <- paste0(ylab, " | ", cond)
  v1 <- df[[ev[1L]]]
  grp <- if (length(ev) == 2L) factor(df[[ev[2L]]])
  ylim <- range(df$lower__, df$upper__)

  if (is.numeric(v1)) {
    graphics::plot(range(v1), ylim, type = "n", xlab = ev[1L],
                   ylab = ylab)
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
    ft <- fitted(x)
    graphics::plot(ft, r, xlab = "Fitted values",
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
  yrep <- t(as.matrix(na_unpad(
    object, simulate(object, nsim = ndraws, re.form = re.form))))
  fun <- get(paste0("ppc_", type), envir = asNamespace("bayesplot"))
  fun(as.numeric(y), yrep, ...)
}
