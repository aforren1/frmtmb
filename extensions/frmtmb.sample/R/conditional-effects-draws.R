#' @exportS3Method brms::conditional_effects
#' @export
conditional_effects.frmtmb_draws <- function(x, effects = NULL,
                                             resp = NULL, dpar = NULL,
                                             resolution = 100,
                                             prob = 0.95, ndraws = NULL,
                                             re_formula = NA,
                                             conditions = list(),
                                             data = NULL, ...) {
  # same recycling trap as the fit method: prob and resolution index a
  # grid, and a length-2 value silently alternated or truncated it
  check_probability(prob, "prob")
  check_count(resolution, "resolution", min = 1L)
  if (!is.null(ndraws)) check_count(ndraws, "ndraws", min = 1L)
  check_named_list(conditions, "conditions", "conditions = list(z = 0)")
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
  fit <- draws_base_fit(x)
  resp <- resp %||% names(fit$spec$responses)[1L]
  rspec <- fit$spec$responses[[resp]]
  ce_structure_check(rspec)
  if (length(fit$frame$re_blocks) &&
      !any(startsWith(colnames(x$draws), "b["))) {
    stop("conditional_effects() on draws from frm_sample(laplace = ",
         "TRUE) cannot rebuild the per-draw parameter vectors: the ",
         "inner parameters were integrated out, so the draws columns ",
         "do not align with the model's parameter template. Resample ",
         "without laplace = TRUE, or call conditional_effects() on the ",
         "fit itself", call. = FALSE)
  }
  categorical <- ce_cats_display(rspec, dpar)
  dpar <- dpar %||% if ("mu" %in% names(rspec$dpars)) "mu" else
    rspec$primary_dpars[1]
  lp <- find_linpred(fit, resp, dpar)
  gb <- ce_grids_build(fit, rspec, lp, effects, resp, dpar, resolution,
                       conditions, data)

  # the same per-parameter-vector grid evaluation the bootstrap band
  # runs per refit, here run per posterior draw. The estimate and the
  # band both come from the draws, so a nonlinear predictor or a
  # nominal category display needs no delta method here
  idx <- draws_par_index(x$fit)
  rows <- draws_subsample(x, ndraws)
  f1 <- draws_fit_at(x, rows[1L], idx)
  lens <- vapply(gb$grids, function(g) {
    length(ce_boot_one(f1, g$nd, categorical, resp, dpar, re_formula))
  }, 1L)
  offsets <- cumsum(c(0L, lens))
  M <- matrix(NA_real_, length(rows), sum(lens))
  for (i in seq_along(rows)) {
    fi <- if (i == 1L) f1 else draws_fit_at(x, rows[i], idx)
    M[i, ] <- unlist(lapply(gb$grids, function(g) {
      ce_boot_one(fi, g$nd, categorical, resp, dpar, re_formula)
    }), use.names = FALSE)
  }

  dfs_by_eff <- list()
  for (gi in seq_along(gb$grids)) {
    g <- gb$grids[[gi]]
    seg <- M[, offsets[gi] + seq_len(lens[gi]), drop = FALSE]
    est <- colMeans(seg)
    lo <- ce_pctl(seg, (1 - prob) / 2)
    up <- ce_pctl(seg, 1 - (1 - prob) / 2)
    se <- apply(seg, 2, stats::sd)
    if (categorical) {
      # ce_boot_one() flattened the n x categories probability matrix
      # column-major, so category k occupies positions (k-1)*n + 1:n
      cats <- colnames(predict(f1, newdata = g$nd, type = "response",
                               resp = resp, re.form = re_formula))
      df <- do.call(rbind, lapply(seq_along(cats), function(k) {
        d <- g$nd[g$ev]
        kk <- (k - 1L) * g$n + seq_len(g$n)
        d$estimate__ <- est[kk]
        d$se__ <- se[kk]
        d$lower__ <- lo[kk]
        d$upper__ <- up[kk]
        d$cats__ <- factor(cats[k], levels = cats)
        d
      }))
    } else {
      df <- g$nd[g$ev]
      df$estimate__ <- est
      df$se__ <- se
      df$lower__ <- lo
      df$upper__ <- up
    }
    if (length(gb$cond_sets) > 1L) {
      df$cond__ <- names(gb$cond_sets)[g$ci] %||% as.character(g$ci)
    }
    dfs_by_eff[[g$eff]] <- c(dfs_by_eff[[g$eff]], list(df))
  }
  ce_finalize(dfs_by_eff, gb$effects, rspec, resp, dpar, "posterior",
              gb$base, categorical, gb$cond_sets, gb$groups)
}
