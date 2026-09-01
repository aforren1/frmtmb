# Deletion diagnostics: leave-one-group (or one-observation) refits.

#' Influence measures by case deletion
#'
#' Refits the model with one group (or observation) left out at a time,
#' warm-started at the full-data estimates, and collects the
#' fixed-effect and covariance-parameter changes. [cooks.distance()] on
#' the result gives the scaled fixed-effect displacement (calling it on
#' the fit itself runs `influence()` first); `dfbeta()` and `dfbetas()`
#' give the per-unit coefficient changes, raw and scaled by the
#' coefficient standard errors (the lme4 influence surface).
#'
#' @param model A `frmtmb_fit`.
#' @param groups Name of a random-effect grouping factor (see
#'   [ngrps()]) to delete level-wise; `NULL` deletes single
#'   observations (refuses for large n unless `force = TRUE`).
#' @param data The original model data; defaults to the stored model
#'   frame, which works unless the formula uses variables that are not
#'   stored raw (e.g. inside `poly()`).
#' @param force Allow observation-wise deletion for n > 500.
#' @param ... Unused.
#' @return A `frmtmb_influence` object: `fixed` and `theta` matrices
#'   (one row per deleted unit) plus the full-data reference.
#' @examples
#' \donttest{
#' set.seed(7)
#' dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
#' dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#' infl <- influence(fit, groups = "g")
#' cooks.distance(infl)
#' }
#' @export
influence.frmtmb_fit <- function(model, groups = NULL, data = NULL,
                                 force = FALSE, ...) {
  data <- data %||% model$frame$data_frame
  if (!is.null(groups)) {
    gv <- data[[groups]]
    if (is.null(gv)) {
      stop("Grouping variable '", groups, "' not found in the data",
           call. = FALSE)
    }
    units <- as.character(unique(gv))
    unit_rows <- lapply(units, function(u) which(as.character(gv) != u))
  } else {
    if (nrow(data) > 500 && !force) {
      stop("Observation-wise influence for n > 500 is expensive; pass ",
           "groups = or force = TRUE", call. = FALSE)
    }
    units <- as.character(seq_len(nrow(data)))
    unit_rows <- lapply(seq_len(nrow(data)), function(i) {
      setdiff(seq_len(nrow(data)), i)
    })
  }

  # get_coef (not unlist(fixef)): aligned with vcov() rows, so mapped
  # constant-dpar coefficients stay out and cooks.distance/dfbetas
  # conform with the covariance matrix
  full_fe <- get_coef.frmtmb_fit(model)
  full_th <- model$estimates$theta
  ctl <- model$control %||% frmtmb_control()
  # one refit per unit: inheriting a verbose fit's control here would
  # print hundreds of stage lines for a result that is already a table
  ctl$verbose <- FALSE
  fe <- matrix(NA_real_, length(units), length(full_fe),
               dimnames = list(units, names(full_fe)))
  th <- matrix(NA_real_, length(units), length(full_th),
               dimnames = list(units, names(outer_theta_names(model))))

  for (i in seq_along(units)) {
    fit_i <- tryCatch(suppressWarnings({
      frame_i <- assemble_frame(model$spec,
                                data[unit_rows[[i]], , drop = FALSE],
                                sparse_x = isTRUE(ctl$sparse_x))
      tpl <- frame_i$par_template
      for (cp in setdiff(names(tpl), "b")) {
        if (length(model$estimates[[cp]]) == length(tpl[[cp]])) {
          tpl[[cp]][] <- model$estimates[[cp]]
        }
      }
      fit_assembled(model$spec, frame_i, model$bform, model$call,
                    REML = model$REML, start = NULL, control = ctl,
                    se = FALSE, lower = model$lower, upper = model$upper,
                    priors = model$priors,
                    quadrature = isTRUE(model$quadrature),
                    template = tpl)
    }), error = function(e) NULL)
    if (is.null(fit_i)) next
    fe_i <- get_coef.frmtmb_fit(fit_i)
    fe[i, names(fe_i)] <- fe_i
    if (length(full_th) &&
        length(fit_i$estimates$theta) == length(full_th)) {
      th[i, ] <- fit_i$estimates$theta
    }
  }
  structure(list(fixed = fe, theta = th, fixed_full = full_fe,
                 theta_full = full_th, groups = groups, fit = model),
            class = "frmtmb_influence")
}

# theta labels for the influence table
outer_theta_names <- function(fit) {
  th <- fit$estimates$theta
  if (!length(th)) return(character(0))
  stats::setNames(th, paste0("theta_", seq_along(th)))
}

#' @export
print.frmtmb_influence <- function(x, n = 6, ...) {
  cd <- cooks.distance(x)
  cat("Case-deletion influence over ",
      if (is.null(x$groups)) "observations" else
        paste0("levels of '", x$groups, "'"),
      " (", nrow(x$fixed), " refits)\n\n", sep = "")
  ord <- order(cd, decreasing = TRUE)
  top <- utils::head(ord, n)
  print(data.frame(unit = rownames(x$fixed)[top],
                   cooks_d = signif(cd[top], 4)), row.names = FALSE)
  invisible(x)
}

#' @rdname influence.frmtmb_fit
#' @export
cooks.distance.frmtmb_fit <- function(model, ...) {
  cooks.distance(influence(model, ...))
}

#' @rdname influence.frmtmb_fit
#' @export
dfbeta.frmtmb_influence <- function(model, ...) {
  # stats convention: the change when the unit is deleted,
  # full-data estimate minus leave-one-out estimate
  -sweep(model$fixed, 2, model$fixed_full)
}

#' @rdname influence.frmtmb_fit
#' @export
dfbetas.frmtmb_influence <- function(model, ...) {
  sweep(dfbeta(model), 2, sqrt(diag(vcov(model$fit))), `/`)
}

#' @export
cooks.distance.frmtmb_influence <- function(model, ...) {
  V <- vcov(model$fit)
  p <- ncol(V)
  Vi <- solve(V)
  d <- sweep(model$fixed, 2, model$fixed_full)
  out <- vapply(seq_len(nrow(d)), function(i) {
    di <- d[i, ]
    if (anyNA(di)) return(NA_real_)
    drop(t(di) %*% Vi %*% di) / p
  }, numeric(1))
  stats::setNames(out, rownames(model$fixed))
}
