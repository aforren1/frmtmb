# The draws method's DISPLAY is core's display. Every column, every
# grid decision and every default here comes from the shared engine
# (`?frmtmb::"frmtmb-sampling-api"`), because what differs between a
# fit and a posterior is where the curves come from, not what is being
# drawn. Assembling a second frame here is what left the two surfaces
# disagreeing in both directions: this method silently inherited every
# grid change core made through ce_grids_build() while keeping its own
# older column set, cond__ rule and ordinal key.

#' The expected category NUMBER of each grid row, per draw.
#'
#' brms's ordinal default display, `sum_k k p_k`. Taken from the drawn
#' probability matrices rather than from a delta method, so its band is
#' a posterior quantile like every other band on this surface.
#'
#' @noRd
ce_cat_mean <- function(seg, n, K) {
  # ce_boot_one() flattens the n x K matrix column-major, so category k
  # occupies columns (k - 1) n + 1:n of every draw's row
  out <- matrix(0, nrow(seg), n)
  for (k in seq_len(K)) {
    out <- out + k * seg[, (k - 1L) * n + seq_len(n), drop = FALSE]
  }
  out
}

#' @exportS3Method brms::conditional_effects
#' @export
conditional_effects.frmtmb_draws <- function(x, effects = NULL,
                                             resp = NULL, dpar = NULL,
                                             resolution = 100,
                                             prob = 0.95, ndraws = NULL,
                                             re_formula = NA,
                                             conditions = list(),
                                             data = NULL,
                                             int_conditions = list(),
                                             categorical = NULL,
                                             seed = NULL, ...) {
  # same recycling trap as the fit method: prob and resolution index a
  # grid, and a length-2 value silently alternated or truncated it
  check_probability(prob, "prob")
  check_count(resolution, "resolution", min = 1L)
  if (!is.null(ndraws)) check_count(ndraws, "ndraws", min = 1L)
  check_named_list(conditions, "conditions", "conditions = list(z = 0)")
  check_named_list(int_conditions, "int_conditions",
                   "int_conditions = list(z = c(-1, 0, 1))")
  dots <- list(...)
  re_formula <- ce_re_formula(re_formula, dots)
  if (!is.null(dots$method)) {
    stop("conditional_effects() on draws has no method =: the curves ",
         "ARE posterior expected-response draws. For predictive bands, ",
         "quantile posterior_predict() over your own grid", call. = FALSE)
  }
  if (!is.null(dots$band)) {
    stop("conditional_effects() on draws has no band =: the band IS ",
         "the posterior quantile band of the drawn curves, so there is ",
         "no wald/profile/boot choice to make", call. = FALSE)
  }
  # core's own dots reader, not a second one: it accepts the two
  # allow_new_levels spellings the fit method accepts and reports the
  # rest. Hand-rolling the check here made the SAME argument work on a
  # fit and warn on draws, which is the argument-surface divergence this
  # method exists to remove.
  allow_new_levels <- ce_dots(dots)
  fit <- draws_base_fit(x)
  resp <- resp %||% names(fit$spec$responses)[1L]
  rspec <- fit$spec$responses[[resp]]
  ce_structure_check(rspec)
  if (length(fit$frame[["re_blocks"]]) &&
      !any(startsWith(colnames(x$draws), "b["))) {
    stop("conditional_effects() on draws from frm_sample(laplace = ",
         "TRUE) cannot rebuild the per-draw parameter vectors: the ",
         "inner parameters were integrated out, so the draws columns ",
         "do not align with the model's parameter template. Resample ",
         "without laplace = TRUE, or call conditional_effects() on the ",
         "fit itself", call. = FALSE)
  }
  # the three displays, resolved core's way: `categorical =` used to be
  # neither honored nor a formal here, so the expected category number
  # could not be asked for at all
  kind <- ce_display_kind(rspec, dpar, categorical)
  categorical <- identical(kind, "cats")
  cats_mean <- identical(kind, "cats_mean")
  poly <- categorical || cats_mean
  dpar_given <- !is.null(dpar)
  dpar <- dpar %||% if ("mu" %in% names(rspec$dpars)) "mu" else
    rspec$primary_dpars[1]
  lp <- find_linpred(fit, resp, dpar)
  # WHAT IS PREDICTED, which is not always the dpar the display is
  # labeled with. Passing the resolved `dpar` here instead drew the mu
  # predictor on every family whose mean is not the inverse link of mu:
  # (1 - zi) too little on a zero-inflated fit, one component's mean on
  # a mixture, and a negative mean on a trunc(lb = 0) response. Core's
  # fit method has made this call since 0.52.0; the decision is shared
  # rather than re-derived so the two surfaces cannot drift again.
  pred_dpar <- ce_pred_dpar(rspec, dpar, dpar_given, categorical,
                            cats_mean)
  # re_formula = NULL means a NEW group here exactly as it does on the
  # fit side: the grid's grouping column says so with NA. Where the fit
  # method draws that group's effects per BOOTSTRAP replicate, this one
  # draws them per POSTERIOR DRAW, from that draw's own theta - which
  # is brms's construction, and the reason its curve for an unobserved
  # group is not the population curve under another name.
  pop_level <- !inherits(re_formula, "formula") &&
    length(re_formula) == 1L && is.na(re_formula)
  na_vars <- if (!pop_level) ce_group_vars(fit) else character(0)
  gb <- ce_grids_build(fit, rspec, lp, effects, resp, dpar, resolution,
                       conditions, data, int_conditions, na_vars)
  anl <- allow_new_levels ||
    length(setdiff(na_vars, names(conditions))) > 0L
  nspec <- if (length(na_vars)) {
    ce_new_level_spec(fit, na_vars, gb$base)
  } else {
    list()
  }
  # the returned frame reports the NA level; the PREDICTION runs on the
  # placeholder level the design can map, whose coefficients each draw
  # overwrites
  egrids <- if (length(nspec)) {
    ce_boot_grids(gb$grids, nspec, gb$base)
  } else {
    gb$grids
  }
  if (!is.null(seed) && length(nspec)) {
    if (!exists(".Random.seed", envir = globalenv())) stats::runif(1)
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
            add = TRUE)
    set.seed(seed)
  }

  # the same per-parameter-vector grid evaluation the bootstrap band
  # runs per refit, here run per posterior draw. The estimate and the
  # band both come from the draws, so a nonlinear predictor or a
  # nominal category display needs no delta method here
  idx <- draws_par_index(x$fit)
  rows <- draws_subsample(x, ndraws)
  f1 <- draws_fit_at(x, rows[1L], idx)
  cats <- if (poly) {
    colnames(predict(f1, newdata = egrids[[1L]]$nd, type = "response",
                     resp = resp, re.form = re_formula,
                     allow_new_levels = anl))
  }
  lens <- vapply(egrids, function(g) {
    length(ce_boot_one(f1, g$nd, poly, resp, pred_dpar, re_formula, anl))
  }, 1L)
  offsets <- cumsum(c(0L, lens))
  M <- matrix(NA_real_, length(rows), sum(lens))
  for (i in seq_along(rows)) {
    fi <- if (i == 1L) f1 else draws_fit_at(x, rows[i], idx)
    if (length(nspec)) fi <- ce_draw_new_levels(fi, nspec)
    M[i, ] <- unlist(lapply(egrids, function(g) {
      ce_boot_one(fi, g$nd, poly, resp, pred_dpar, re_formula, anl)
    }), use.names = FALSE)
  }

  # brms's frame always carries cond__, with one level when there is
  # one condition set: its own plot() facets on the column, so a ported
  # faceting call has nothing to facet on without it
  clev <- names(gb$cond_sets) %||% as.character(seq_along(gb$cond_sets))
  dfs_by_eff <- list()
  for (gi in seq_along(gb$grids)) {
    g <- gb$grids[[gi]]
    seg <- M[, offsets[gi] + seq_len(lens[gi]), drop = FALSE]
    if (cats_mean) seg <- ce_cat_mean(seg, g$n, length(cats))
    est <- colMeans(seg)
    lo <- ce_pctl(seg, (1 - prob) / 2)
    up <- ce_pctl(seg, 1 - (1 - prob) / 2)
    se <- apply(seg, 2, stats::sd)
    cond <- factor(clev[g$ci], levels = clev)
    if (categorical) {
      df <- do.call(rbind, lapply(seq_along(cats), function(k) {
        d <- ce_frame(g$nd, g$ev, g$v2, cond,
                      cats = factor(cats[k], levels = cats))
        kk <- (k - 1L) * g$n + seq_len(g$n)
        d$estimate__ <- est[kk]
        d$se__ <- se[kk]
        d$lower__ <- lo[kk]
        d$upper__ <- up[kk]
        d
      }))
    } else {
      df <- ce_frame(g$nd, g$ev, g$v2, cond)
      df$estimate__ <- est
      df$se__ <- se
      df$lower__ <- lo
      df$upper__ <- up
    }
    dfs_by_eff[[g$eff]] <- c(dfs_by_eff[[g$eff]], list(df))
  }
  ce_finalize(dfs_by_eff, gb$effects, rspec, resp, dpar, "posterior",
              gb$base, categorical, gb$cond_sets, gb$groups, categorical)
}
