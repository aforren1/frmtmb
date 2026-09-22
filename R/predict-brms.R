# brms's predict(): a summary of the PREDICTIVE distribution (item
# 2.6d).
#
# brms summarizes the posterior predictive distribution, so observation
# noise is in it and the interval is much wider than the one around a
# fitted value. frmtmb's predict() used to be glmmTMB's: the LINEAR
# PREDICTOR by default, with a `type` vocabulary. On a lognormal fit
# that returned a different scale from brms with nothing said, which is
# why the two could not both keep the name. brms is the tiebreaker, so
# `predict()` is brms's and the linear predictor moved to
# `frm_linpred()`, which is the old function under a name that says what
# it returns.
#
# The frequentist construction. A maximum-likelihood fit has no
# posterior, so the draws are SIMULATED: each replicate draws the outer
# parameter vector from its asymptotic normal law,
# `N(theta_hat, vcov(fit, full = TRUE))`, evaluates every distributional
# parameter there, and draws a response from the family's own
# simulator, the one `simulate()` and `posterior_predict()` use. That
# is a parametric bootstrap of the predictive distribution, and it
# carries both sources of spread brms's carries: the observation noise
# and the uncertainty in the estimates. `param_uncertainty = FALSE`
# drops the second and simulates at the estimates alone, which is the
# plug-in predictive distribution; `dev/shapes-findings.md` reports what
# each one covers.
#
# The random effects stay at their conditional modes, as they do
# everywhere else in this package (`fitted()`, `frm_linpred()` and
# `residuals()` are all conditional on the modes). Their own conditional
# variance is therefore NOT in the interval; on a design with few
# observations per group that is the term the coverage measurement is
# sensitive to, and the section "Which uncertainty is in the interval"
# of `?predict.frmtmb_fit` says so.

#' Arguments `predict()` carried before it became brms's, and the brms
#' arguments it cannot answer.
#'
#' `type` is named rather than left to the unknown-argument message
#' because it was the CENTRE of the old vocabulary: a script written
#' against 0.60.0 passes it, and silently getting a predictive summary
#' back where the linear predictor was asked for is the failure this
#' whole item exists to stop.
#'
#' @noRd
predict_type_retired <- c(
  type = paste("predict() is brms's predictive summary now and has no",
               "`type`. The linear predictor and the glmmTMB scale",
               "vocabulary moved to frm_linpred(), unchanged:",
               "frm_linpred(object, type = \"link\") is what",
               "predict(object) used to return, and fitted(object) is",
               "the expected response"),
  se.fit = paste("predict() summarizes a predictive distribution and",
                 "reports its spread in the Est.Error column. For the",
                 "delta-method standard error of a linear predictor use",
                 "frm_linpred(object, se.fit = TRUE)"),
  dpar = paste("brms's predict() draws the RESPONSE, so it has no dpar.",
               "fitted(object, dpar = ) reports a distributional",
               "parameter and frm_linpred(object, dpar = ) its linear",
               "predictor"),
  scale = paste("brms spells this `type` on fitted() and has neither on",
                "predict(). frm_linpred() takes the scale"),
  re.form = paste("lme4's spelling, no longer accepted. brms is the",
                  "tiebreaker on a name, so this setting is",
                  "`re_formula` here as it is in brms. Pass",
                  "re_formula ="),
  allow.new.levels = paste("lme4's spelling, no longer accepted. brms",
                           "spells it `allow_new_levels`, and so does",
                           "this. Pass allow_new_levels ="),
  nlpar = paste("a non-linear parameter is reached through `dpar` in",
                "this package, on fitted() and frm_linpred()")
)

#' Predictions from a frmtmb fit
#'
#' @description
#' brms's `predict()`: a summary of the PREDICTIVE distribution of the
#' response, with observation noise included, in the columns `Estimate`,
#' `Est.Error`, `Q2.5` and `Q97.5`. It is a much wider interval than
#' [fitted()]'s, which summarizes the expected response.
#'
#' A maximum-likelihood fit has no posterior, so the draws are
#' simulated. Each replicate draws the outer parameter vector from its
#' asymptotic normal law, `N(theta_hat, vcov(object, full = TRUE))`,
#' evaluates every distributional parameter there, and draws a response
#' from the family's own simulator, the one [simulate()] uses. The
#' interval is the empirical quantile of those draws, so it carries the
#' family's skew and its discreteness: a poisson predictive interval is
#' on the counts.
#'
#' @section Which uncertainty is in the interval:
#' Two of the three sources.
#'
#' * Observation noise, from the family's simulator. This is almost
#'   always the largest term.
#' * Uncertainty in the estimates, from the asymptotic covariance of
#'   the outer parameter vector. `param_uncertainty = FALSE` drops it
#'   and simulates at the estimates alone.
#' * Uncertainty in the random-effect modes of a level the fit SAW is
#'   NOT included. The draw is conditional on them, the convention
#'   every other method here follows. On a design with few
#'   observations per grouping level that term is not small, and the
#'   interval is then narrower than a fully marginal one.
#'   `dev/shapes-findings.md` reports the measured coverage for a
#'   design with and without a grouping factor.
#'
#' A level the fit did NOT see is different: there is no mode to
#' condition on, so `allow_new_levels = TRUE` draws its effect from the
#' block's own estimated covariance, once per replicate, at that
#' replicate's parameters. This is brms's `sample_new_levels =
#' "gaussian"` and it is the whole between-group variance, so the
#' interval at an unseen level is wider than at a known one.
#'
#' A fit whose covariance could not be recovered from the Hessian has no
#' second term to draw: the simulation falls back to the estimates and
#' warns, rather than drawing from a matrix of `NaN`.
#'
#' Under `REML = TRUE` (or `frmtmb_control(profile = TRUE)`) the fixed
#' effects are integrated out of the outer problem, so they are not in
#' `vcov(object, full = TRUE)` and are not drawn either. The interval
#' there carries the covariance parameters' uncertainty and the
#' observation noise, and holds the coefficients fixed.
#'
#' @section Ordinal and categorical responses:
#' There is no mean to summarize, so `predict()` returns the simulated
#' proportion of each category, one column per category named
#' `P(Y = k)`, which is brms's shape. [fitted()] returns the MODELLED
#' category probabilities with their standard errors instead, which is
#' the smoother quantity.
#'
#' @param object A `frmtmb_fit`.
#' @param newdata Optional data frame to predict on. Defaults to the
#'   training data.
#' @param re_formula `NULL` (default) keeps the random effects, so the
#'   draw is conditional on the modes; `NA` or `~0` draws at the
#'   population level.
#' @param transform A function applied to the draws before they are
#'   summarized, as in brms.
#' @param resp For multivariate fits: which response, or `NULL`
#'   (default) for all of them, which is brms's `nrow x 4 x nresp`
#'   array with the responses named on the third dimension.
#' @param negative_rt brms's sign convention for its `wiener` family.
#'   Refused, with the reason: this package's evidence-accumulation
#'   families return the response and the time as they declare them.
#' @param ndraws Number of simulated replicates. Defaults to 1000.
#'   brms thins its posterior with this; here it sets how many draws
#'   are taken, so a larger number narrows the Monte Carlo error of the
#'   quantile columns and costs proportionally more.
#' @param draw_ids,sort,cores Refused by name. There are no stored
#'   draws to index, rows always come back in the order of the data,
#'   and the simulation is not parallelized.
#' @param ntrys Rejection-sampling attempts per row for a truncated
#'   response, brms's spelling of the family simulator's `max_iter`.
#' @param summary If `FALSE`, the `ndraws x nrow` matrix of simulated
#'   draws instead of the summary.
#' @param robust If `TRUE`, the median and the median absolute
#'   deviation of the draws instead of the mean and the standard
#'   deviation. brms's argument, and it is answerable here because these
#'   draws exist.
#' @param probs Probabilities of the quantile columns.
#' @param allow_new_levels Predict rows whose grouping-factor level the
#'   fit never saw, instead of erroring. Each replicate draws that
#'   level's effect from the block's estimated covariance.
#' @param sample_new_levels brms's argument. `"gaussian"`, which is
#'   what happens, or `NULL`. `"uncertainty"` and `"old_levels"`
#'   resample the posterior draws of the levels that were seen, and a
#'   maximum-likelihood fit has none, so both are refused by name.
#' @param param_uncertainty If `FALSE`, simulate at the estimates alone
#'   (the plug-in predictive distribution) instead of drawing the
#'   parameters first.
#' @param ... Refused. An argument this method does not have is an
#'   error naming it, and `type`, `se.fit`, `dpar` and `scale` are
#'   refused with the name of the function that took over each one.
#' @return With `summary = TRUE` (the default) an `nrow x 4` matrix with
#'   the columns `Estimate`, `Est.Error` and one per entry of `probs`;
#'   for an ordinal or categorical response an `nrow x K` matrix of
#'   simulated category proportions; for a multivariate fit an
#'   `nrow x 4 x nresp` array with the responses named. With
#'   `summary = FALSE` the `ndraws x nrow` matrix of draws, or
#'   `ndraws x nrow x nresp` for a multivariate fit.
#' @seealso [fitted.frmtmb_fit()] for the expected response,
#'   [frm_linpred()] for the linear predictor and the glmmTMB scale
#'   vocabulary, [simulate.frmtmb_fit()] for the draws themselves, and
#'   [frmtmb-scales].
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100))
#' dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x))
#' fit <- frm(bf(y ~ x) + poisson(), data = dd)
#'
#' # the predictive summary: observation noise is in it
#' head(predict(fit, ndraws = 200))
#' # against the expected response, which is much tighter
#' head(fitted(fit))
#'
#' # the draws themselves
#' dim(predict(fit, summary = FALSE, ndraws = 50))
#'
#' # the linear predictor is frm_linpred() now
#' head(frm_linpred(fit, type = "link"))
#' @export
predict.frmtmb_fit <- function(object, newdata = NULL, re_formula = NULL,
                               transform = NULL, resp = NULL,
                               negative_rt = FALSE, ndraws = NULL,
                               draw_ids = NULL, sort = FALSE, ntrys = 100L,
                               cores = NULL, summary = TRUE,
                               robust = FALSE, probs = c(0.025, 0.975),
                               ..., allow_new_levels = FALSE,
                               sample_new_levels = NULL,
                               param_uncertainty = TRUE) {
  frm_check_dots(..., .unsupported = predict_type_retired)
  require_fitted(object, "predict()")
  check_flag(summary, "summary")
  check_flag(robust, "robust")
  check_flag(negative_rt, "negative_rt")
  check_flag(allow_new_levels, "allow_new_levels")
  check_flag(param_uncertainty, "param_uncertainty")
  check_re_form(re_formula)
  if (negative_rt) {
    frm_stop("predict(negative_rt = TRUE) is brms's sign convention for ",
             "its wiener family, which codes the lower boundary as a ",
             "negative reaction time. frmtmb's evidence-accumulation ",
             "families return the response and the time as they declare ",
             "them; see the family's own documentation in frmtmb.eam",
             call. = FALSE)
  }
  if (!is.null(draw_ids)) {
    frm_stop("predict() cannot honor `draw_ids`: there are no stored draws ",
             "to index. The draws here are simulated by the call itself, so ",
             "`ndraws` sets how many are taken", call. = FALSE)
  }
  if (!isFALSE(sort)) {
    frm_stop("predict() cannot honor `sort`: rows come back in the order ",
             "of the data, always", call. = FALSE)
  }
  if (!is.null(cores)) {
    frm_stop("predict() cannot honor `cores`: the simulation is not ",
             "parallelized. Lower `ndraws` if a call is too slow",
             call. = FALSE)
  }
  if (!is.null(sample_new_levels) &&
      !identical(sample_new_levels, "gaussian")) {
    frm_stop("predict() honors sample_new_levels = \"gaussian\" only. ",
             "An unseen level's effect is drawn from its block's own ",
             "estimated covariance, which is what brms's \"gaussian\" ",
             "does. \"uncertainty\" and \"old_levels\" resample the ",
             "POSTERIOR draws of the levels that were seen, and a ",
             "maximum-likelihood fit has no such draws", call. = FALSE)
  }
  if (!is.null(newdata) && !is.data.frame(newdata)) {
    frm_stop("predict(): `newdata` must be a data frame, or NULL to draw ",
             "for the training rows, not ", arg_desc(newdata),
             call. = FALSE)
  }
  # check BEFORE the coercion: as.integer(2.5) is 2, so checking the
  # coerced value accepted a fraction and silently halved the request
  ndraws <- ndraws %||% 1000L
  check_count(ndraws, "ndraws", min = 1L)
  ndraws <- as.integer(ndraws)
  check_count(ntrys, "ntrys", min = 1L)
  # brms answers a multivariate fit for EVERY response, in an
  # n x 4 x nresp array with the responses named. Taking the first one
  # and saying nothing is the failure mode this whole item exists to
  # stop, so `resp` is only a FILTER now, not a silent default.
  rnames <- names(object$spec$responses)
  if (is.null(resp)) {
    resp <- rnames
  } else {
    if (!is.character(resp) || anyNA(resp)) {
      frm_stop("predict(): `resp` must name one or more responses, or NULL ",
               "for all of them, not ", arg_desc(resp), call. = FALSE)
    }
    for (r in resp) {
      if (!r %in% rnames) stop_unknown_response(object, r)
    }
  }
  rspecs <- object$spec$responses[resp]
  if (length(rspecs) > 1L) predict_mv_refuse(object, rspecs)
  ds <- predict_simulate(object, rspecs, newdata, re_formula,
                         allow_new_levels, ndraws, ntrys, param_uncertainty)
  if (!is.null(transform)) ds <- lapply(ds, match.fun(transform))
  if (!summary) {
    return(if (length(ds) == 1L) ds[[1L]] else predict_stack(ds))
  }
  per <- lapply(seq_along(ds), function(i) {
    predict_summarize_one(object, rspecs[[i]], ds[[i]], probs, robust)
  })
  if (length(per) == 1L) return(per[[1L]])
  out <- array(unlist(per), c(nrow(per[[1L]]), ncol(per[[1L]]), length(per)),
               dimnames = list(rownames(per[[1L]]), colnames(per[[1L]]),
                               names(ds)))
  out
}

#' brms's summary of one response's simulated draws.
#'
#' @noRd
predict_summarize_one <- function(object, rspec, d, probs, robust) {
  if (fam_is_category_valued(rspec$family)) {
    return(predict_category_props(object, rspec, d))
  }
  if (length(dim(d)) == 3L) {
    # one summary layer per column of a matrix-valued response
    per <- lapply(seq_len(dim(d)[3L]), function(k) {
      brms_summarize_draws(d[, , k, drop = TRUE], probs, robust)
    })
    return(array(unlist(per), c(nrow(per[[1L]]), ncol(per[[1L]]),
                                length(per)),
                 dimnames = list(dimnames(d)[[2L]], colnames(per[[1L]]),
                                 dimnames(d)[[3L]])))
  }
  brms_summarize_draws(d, probs, robust)
}

#' The responses a multivariate `predict()` cannot stack into one
#' array, and the spelling that answers for them one at a time.
#'
#' A category-valued response summarizes to `n x K` proportions and a
#' matrix-valued one to `n x 4 x ncol`, so neither has the `n x 4` face
#' brms's third dimension is made of. Refusing names the way through
#' rather than returning an array whose layers mean different things.
#'
#' @noRd
predict_mv_refuse <- function(object, rspecs) {
  bad <- names(rspecs)[vapply(rspecs, function(rs) {
    fam_is_category_valued(rs$family)
  }, TRUE)]
  if (length(bad)) {
    frm_stop("predict() summarizes a multivariate fit into one ",
             "nrow x 4 x nresp array, and the response(s) ",
             paste0("'", bad, "'", collapse = ", "), " summarize to ",
             "category proportions instead, which do not stack with it. ",
             "Ask for one response at a time: predict(object, resp = \"",
             bad[1L], "\")", call. = FALSE)
  }
  invisible(NULL)
}

#' The `ndraws x nrow x nresp` array brms returns for a multivariate
#' `predict(summary = FALSE)`.
#'
#' @noRd
predict_stack <- function(ds) {
  if (any(vapply(ds, function(d) length(dim(d)) == 3L, TRUE))) {
    frm_stop("predict(summary = FALSE) on a multivariate fit with a ",
             "matrix-valued response would need four dimensions. Ask for ",
             "one response at a time with resp =", call. = FALSE)
  }
  array(unlist(ds), c(nrow(ds[[1L]]), ncol(ds[[1L]]), length(ds)),
        dimnames = list(NULL, colnames(ds[[1L]]), names(ds)))
}

#' Whether the response is a set of unordered or ordered CATEGORIES, so
#' that a predictive summary is a set of proportions rather than a mean.
#'
#' @noRd
fam_is_category_valued <- function(fam) {
  fam[["type"]] %in% c("ordinal", "categorical")
}

#' brms's predictive summary of a category-valued response: the
#' simulated proportion of each category, one column per category.
#'
#' @noRd
predict_category_props <- function(object, rspec, d) {
  lv <- object$frame[["y_levels"]][[rspec$resp_name]]
  K <- if (!is.null(lv)) {
    length(lv)
  } else if (identical(rspec$family[["type"]], "ordinal")) {
    ordinal_ncat(object)
  } else {
    max(c(as.integer(d), object$frame[["y"]][[rspec$resp_name]]))
  }
  out <- t(apply(d, 2L, function(col) {
    # over the replicates that are THERE: a replicate that drew a
    # non-finite dpar is NA, and dividing by the full count would
    # report proportions that do not sum to one
    ok <- col[!is.na(col)]
    if (!length(ok)) return(rep(NA_real_, K))
    tabulate(as.integer(ok), K) / length(ok)
  }))
  colnames(out) <- brms_category_labels(lv, K)
  rownames(out) <- NULL
  out
}

#' The `ndraws x nrow` matrix of simulated responses behind
#' `predict()`, one per response.
#'
#' One code path with `posterior_predict()` on draws: both evaluate the
#' distributional parameters at one parameter vector and hand them to
#' the family's own simulator. What differs is where the parameter
#' vector comes from, a posterior draw there and the asymptotic normal
#' law here.
#'
#' Every response of a multivariate fit is drawn at the SAME parameter
#' draw, so the parameter uncertainty they share is shared in the
#' output too. Each response's own observation noise comes from its own
#' family simulator, so a residual correlation BETWEEN responses is not
#' in the joint draws; each response's marginal summary, which is what
#' `predict()` returns, does not depend on it.
#'
#' @noRd
predict_simulate <- function(object, rspecs, newdata, re_formula,
                             allow_new_levels, ndraws, ntrys,
                             param_uncertainty) {
  av <- list()
  nls <- list()
  for (nm in names(rspecs)) {
    rspec <- rspecs[[nm]]
    fam <- rspec$family
    if (!sim_can(fam)) {
      frm_stop("predict(): family '", fam[["family"]],
               "' has no simulator yet, so there is no predictive ",
               "distribution to summarize. fitted() reports the expected ",
               "response and frm_linpred() the linear predictor",
               sim_note(fam), call. = FALSE,
               package = frm_family_package(fam))
    }
    av[[nm]] <- if (is.null(newdata)) {
      object$frame[["aterm_values"]][[rspec$resp_name]]
    } else if (has_trunc(rspec)) {
      # truncation bounds must follow the newdata rows, or the draws
      # land outside the support the likelihood was normalized on
      aterms_for_newdata(rspec, newdata)
    } else {
      list()
    }
    if (!is.null(newdata) &&
        sim_is_structured(sim_context(object, rspec, list(),
                                      aterms = av[[nm]]))) {
      frm_stop("predict(newdata = ) is not supported for this model: its ",
               "draws are structured (a hidden state sequence, a ",
               "group-level latent class, or a correlated residual) and ",
               "that structure indexes the rows the model was fitted on. ",
               "Drop newdata to predict those rows", call. = FALSE)
    }
    nls[[nm]] <- predict_new_level_spec(object, rspec, newdata, re_formula,
                                        allow_new_levels)
  }
  # Every random quantity a replicate needs from the caller's stream is
  # taken UP FRONT: the parameter draws and one seed per replicate. A
  # replicate then simulates from its own seed. So the parameter draws
  # do not depend on how many rows newdata has, and a row's draws do not
  # depend on the rows AFTER it; they still shift with the rows before
  # it, which share the replicate's stream. The caller's stream is left
  # where the up-front draws put it.
  draw <- predict_par_drawer(object, param_uncertainty, ndraws)
  seeds <- sample.int(.Machine$integer.max, ndraws)
  saved <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", saved, envir = globalenv()), add = TRUE)
  joint <- isTRUE(object$spec[["rescor"]]) && length(rspecs) > 1L
  out <- stats::setNames(vector("list", length(rspecs)), names(rspecs))
  masked <- stats::setNames(vector("list", length(rspecs)), names(rspecs))
  for (s in seq_len(ndraws)) {
    set.seed(seeds[s])
    fs <- draw(s)
    dps <- lapply(stats::setNames(names(rspecs), names(rspecs)), function(nm) {
      predict_dpar_values(fs, rspecs[[nm]], newdata, re_formula,
                          allow_new_levels,
                          predict_new_level_draw(fs, nls[[nm]]))
    })
    # A draw from N(theta_hat, V) can land where a distributional
    # parameter is not finite: on the port's mixture fixture the
    # `(1 | patient)` log standard deviation has variance 7.97e6 in
    # vcov(full = TRUE), so a draw of 1365 makes an unseen level's
    # variance exp(2730) and mu1 is Inf. Such a CELL is NA and the rest
    # of the replicate stands. Dropping the whole replicate, as round 1
    # did, conditioned every other row on a selected subset of parameter
    # draws: a poisson row at x = 4 predicted 15.264 alone and 12.947
    # beside one overflowing row.
    oks <- lapply(dps, dpv_row_finite)
    if (joint) {
      ok <- Reduce(`&`, oks)
      oks <- lapply(oks, function(o) ok)
      ys_all <- predict_rescor_draw(fs, rspecs, dps, ok)
    }
    for (nm in names(rspecs)) {
      ok <- oks[[nm]]
      n <- length(ok)
      if (is.null(masked[[nm]])) masked[[nm]] <- integer(n)
      masked[[nm]] <- masked[[nm]] + !ok
      if (!any(ok)) next
      ys <- if (joint) {
        ys_all[[nm]]
      } else {
        predict_sim_rows(fs, rspecs[[nm]], dps[[nm]], av[[nm]], ok, ntrys)
      }
      if (is.null(ys)) {
        masked[[nm]] <- masked[[nm]] + ok
        next
      }
      if (is.null(out[[nm]])) {
        # a matrix-valued response (multinomial counts, the mvn
        # mixture's columns, lca item codes) gives a ROW per
        # observation, so the draws stack into a
        # draws x observations x columns array, the shape
        # posterior_predict() already uses on draws. The observation
        # margin is unnamed, as brms's posterior_predict() leaves it.
        out[[nm]] <- if (is.matrix(ys)) {
          array(NA_real_, c(ndraws, n, ncol(ys)),
                dimnames = list(NULL, NULL, colnames(ys)))
        } else {
          matrix(NA_real_, ndraws, n)
        }
      }
      if (length(dim(out[[nm]])) == 3L) {
        out[[nm]][s, ok, ] <- ys
      } else {
        out[[nm]][s, ok] <- as.numeric(ys)
      }
    }
  }
  predict_report_masked(out, masked, ndraws)
  out
}

#' Which rows of one replicate have every distributional parameter
#' finite.
#'
#' A dpar list holds vectors of one value per row, a matrix with a row
#' per observation (an ordinal fit's `.cs` offsets), a scalar that
#' applies to every row, and, for a structured family, something that
#' is not numeric at all. A non-numeric entry is not judged: erroring
#' inside a guard is worse than the thing the guard is for.
#'
#' @noRd
dpv_row_finite <- function(dpv) {
  n <- max(vapply(dpv, NROW, 1L), 1L)
  ok <- rep(TRUE, n)
  for (v in dpv) {
    if (!is.numeric(v)) next
    if (is.matrix(v) && nrow(v) == n) {
      ok <- ok & rowSums(!is.finite(v)) == 0
    } else if (length(v) == n) {
      ok <- ok & is.finite(v)
    } else if (!all(is.finite(v))) {
      ok[] <- FALSE
    }
  }
  ok
}

#' The rows of one entry of a per-row list: a vector of one value per
#' row or a matrix with a row per observation is subset, anything else
#' (a scalar, a shared parameter) is passed through.
#'
#' @noRd
subset_rows <- function(x, keep) {
  n <- length(keep)
  if (is.matrix(x) && nrow(x) == n) return(x[keep, , drop = FALSE])
  if (is.atomic(x) && length(x) == n) return(x[keep])
  x
}

#' One replicate's draws of one response, on the rows whose parameters
#' are finite.
#'
#' A structured family's rows are not independent (a hidden state
#' sequence, a group-level class), so they cannot be simulated apart;
#' there a replicate with any non-finite row draws nothing, and every
#' row of it is NA. `NULL` says so to the caller.
#'
#' @noRd
predict_sim_rows <- function(fs, rspec, dpv, av, ok, ntrys) {
  n <- length(ok)
  if (!all(ok)) {
    ctx0 <- sim_context(fs, rspec, list(), aterms = av)
    if (sim_is_structured(ctx0)) return(NULL)
    dpv <- lapply(dpv, subset_rows, keep = ok)
    av <- lapply(av %||% list(), subset_rows, keep = ok)
  }
  sim_draw(sim_context(fs, rspec, dpv, aterms = av, n = sum(ok),
                       extra = fit_extras(fs), max_iter = ntrys))
}

#' One replicate's JOINT draw of the responses of a `rescor` fit.
#'
#' brms draws a multivariate gaussian response from its joint law, so a
#' contrast or a sum across responses in `predict(summary = FALSE)`
#' carries the residual correlation. Drawing each response from its own
#' simulator gave draws that correlated at 0.009 on a fit whose rescor
#' was estimated at 0.7. The correlation matrix is read at THIS
#' replicate's parameters, so its own uncertainty is carried as well.
#' `rescor = TRUE` is refused for every family but gaussian, so the
#' joint law is always the multivariate normal.
#'
#' @noRd
predict_rescor_draw <- function(fs, rspecs, dps, ok) {
  if (!any(ok)) return(NULL)
  R <- rescor_matrix(fs)[names(rspecs), names(rspecs), drop = FALSE]
  L <- chol(R)
  m <- sum(ok)
  Z <- matrix(stats::rnorm(m * length(rspecs)), m) %*% L
  out <- list()
  for (k in seq_along(rspecs)) {
    dp <- dps[[k]]
    mu <- rep_len(as.numeric(dp[["mu"]]), length(ok))[ok]
    sg <- rep_len(as.numeric(dp[["sigma"]]), length(ok))[ok]
    out[[names(rspecs)[k]]] <- mu + sg * Z[, k]
  }
  out
}

#' Say which rows lost cells to a non-finite parameter draw, and stop
#' when there is nothing left to summarize.
#'
#' @noRd
predict_report_masked <- function(out, masked, ndraws) {
  none <- vapply(out, is.null, TRUE)
  if (any(none)) {
    frm_stop("predict(): every simulated replicate of response(s) ",
             paste0("'", names(out)[none], "'", collapse = ", "),
             " drew a distributional parameter that is not finite, so ",
             "there is no predictive distribution to summarize. A ",
             "parameter the draw law covers is barely identified; ",
             "vcov(object, full = TRUE) has its variance. ",
             "param_uncertainty = FALSE simulates at the estimates ",
             "alone", call. = FALSE)
  }
  msgs <- character(0)
  for (nm in names(masked)) {
    cnt <- masked[[nm]]
    hit <- which(cnt > 0L)
    if (!length(hit)) next
    what <- paste0(hit, " (", cnt[hit], " of ", ndraws, ")")
    msgs <- c(msgs, paste0(
      if (length(masked) > 1L) paste0("response '", nm, "': "),
      if (length(hit) == 1L) "row " else "rows ",
      paste(what, collapse = ", ")))
  }
  if (length(msgs)) {
    frm_warning("predict(): a simulated parameter draw was not finite for ",
                paste(msgs, collapse = "; "), ". Those cells are NA and ",
                "each row is summarized over the replicates it has; the ",
                "other rows are not affected. The draw law puts mass ",
                "where a parameter overflows: a row far outside the data, ",
                "or a parameter barely identified, whose variance ",
                "vcov(object, full = TRUE) shows. param_uncertainty = ",
                "FALSE simulates at the estimates alone", call. = FALSE)
  }
  invisible(NULL)
}

#' The unseen grouping levels one replicate has to draw an effect for,
#' and the design rows that load it.
#'
#' brms's `sample_new_levels = "gaussian"`: at a level the fit never
#' saw, the group effect is not zero, it is unknown, and its law is the
#' block's own estimated covariance. Holding it at zero made
#' `predict()`'s `Est.Error` at an unseen level bit-identical to a
#' known level's (0.84358 against 0.84842), and
#' its out-of-sample coverage 0.8618 against a nominal 0.95
#' (`dev/reviews/20260918-shapes.md`, BLOCKER 2).
#'
#' `frm_linpred(se.fit = TRUE)` already held the term, through
#' `lp_extra_var()` and `extra_var_blocks()`. This reads the SAME two
#' functions, so the variance `predict()` draws from and the variance
#' the standard error reports cannot drift apart. The design pieces do
#' not depend on the parameters, so `lp_eta_design()` runs once here
#' and only the block covariance is re-read per replicate.
#'
#' @noRd
predict_new_level_spec <- function(object, rspec, newdata, re_formula,
                                   allow_new_levels) {
  if (!isTRUE(allow_new_levels) || is.null(newdata)) return(NULL)
  use_re <- re_form_keeps(re_formula)
  if (!use_re) return(NULL)
  out <- list()
  hooked <- character(0)
  for (dnm in names(rspec$dpars)) {
    lp <- object$frame[["linpreds"]][[linpred_key(rspec$resp_name, dnm)]]
    if (is.null(lp) || !is.null(lp[["nl_body"]]) ||
          !is.null(lp[["constant"]]) || is.null(lp[["Z"]])) {
      next
    }
    ed <- tryCatch(suppressWarnings(
      lp_eta_design(object, lp, newdata, use_re, TRUE)),
      error = function(e) NULL)
    if (is.null(ed)) next
    if (!length(lp_extra_var(object, ed, use_re)$new_levels)) next
    if (!is.null(dpar_report_hook(rspec$family, dnm, rspec))) {
      # the response scale of such a dpar is not its own link inverse
      # (a mixture weight is the softmax over the component
      # predictors), so an offset added to eta cannot be mapped there
      hooked <- c(hooked, dnm)
      next
    }
    out[[dnm]] <- list(ed = ed, n = ed[["n"]])
  }
  if (length(hooked)) {
    frm_warning("predict(allow_new_levels = TRUE): the unseen level's ",
                "effect is NOT drawn for ",
                paste0("`", hooked, "`", collapse = ", "),
                ", whose response scale is a transform of several ",
                "predictors at once. Those rows are predicted at the ",
                "population level for it", call. = FALSE)
  }
  if (!length(out)) NULL else out
}

#' One replicate's unseen-level effects, as a per-dpar offset on the
#' LINK scale.
#'
#' The covariance is read at THIS replicate's parameters, so the
#' uncertainty in the variance component is carried as well as the
#' spread it implies. A block feeding several dpars shares one draw,
#' keyed as `extra_var_blocks()` keys it.
#'
#' @noRd
predict_new_level_draw <- function(fs, spec) {
  if (is.null(spec)) return(NULL)
  shared <- list()
  out <- list()
  for (dnm in names(spec)) {
    e <- spec[[dnm]]
    bl <- extra_var_blocks(lp_extra_var(fs, e$ed, TRUE)$new_levels, e$n)
    o <- numeric(e$n)
    for (key in names(bl)) {
      B <- bl[[key]]
      if (is.null(shared[[key]])) shared[[key]] <- mvn_draw_cov(B$S)
      o <- o + as.numeric(B$M %*% shared[[key]])
    }
    out[[dnm]] <- o
  }
  out
}

#' One draw from `N(0, S)`, with an eigen fallback for an `S` that is
#' singular, which an `|ID|`-merged block can be.
#'
#' @noRd
mvn_draw_cov <- function(S) {
  d <- nrow(S)
  L <- tryCatch(chol(S), error = function(e) NULL)
  if (!is.null(L)) return(drop(crossprod(L, stats::rnorm(d))))
  ev <- eigen(S, symmetric = TRUE)
  drop(ev$vectors %*% (sqrt(pmax(ev$values, 0)) * stats::rnorm(d)))
}

#' The distributional parameter values one replicate draws from, on
#' their natural scale.
#'
#' The in-sample, all-random-effects case is the one `eval_dpars()`
#' answers directly and is much the cheapest; anything else routes
#' through `frm_linpred()`, which is where `newdata` and `re_formula`
#' are interpreted.
#'
#' @noRd
predict_dpar_values <- function(fit, rspec, newdata, re_formula,
                                allow_new_levels, new_level_off = NULL) {
  resp <- rspec$resp_name
  if (is.null(newdata) && is.null(re_formula)) {
    return(with_cs_offsets(fit, rspec, eval_dpars(fit)[[resp]]))
  }
  dpv <- list()
  for (dnm in names(rspec$dpars)) {
    off <- new_level_off[[dnm]]
    # the NAMES are kept: they are the data's row names, which every
    # other method carries into its output and which as.vector() drops
    if (is.null(off)) {
      dpv[[dnm]] <- drop(frm_linpred(fit, newdata = newdata, dpar = dnm,
                                     resp = resp, re_formula = re_formula,
                                     type = "response",
                                     allow_new_levels = allow_new_levels))
      next
    }
    # an unseen level's drawn effect is an offset on the LINK scale,
    # where the model is additive, so the predictor is read there and
    # the link inverse applied after the offset is in
    lp <- fit$frame[["linpreds"]][[linpred_key(resp, dnm)]]
    eta <- drop(frm_linpred(fit, newdata = newdata, dpar = dnm,
                            resp = resp, re_formula = re_formula,
                            type = "link",
                            allow_new_levels = allow_new_levels))
    dpv[[dnm]] <- lp[["link"]]$linkinv(eta + off)
  }
  with_cs_offsets(fit, rspec, dpv)
}

#' A closure returning one fit-like object per call: the fit itself
#' under `param_uncertainty = FALSE`, and otherwise the fit with its
#' outer parameters replaced by a draw from
#' `N(theta_hat, vcov(full = TRUE))`.
#'
#' The Cholesky factor is taken once, outside the loop, because it is
#' the only expensive part of the draw.
#'
#' @noRd
predict_par_drawer <- function(object, param_uncertainty, ndraws) {
  if (!param_uncertainty) return(function(s) object)
  ds <- fit_draw_space(object)
  V <- ds$V
  L <- if (is.null(V) || !all(is.finite(V))) NULL else {
    tryCatch(chol(V + diag(0, nrow(V))), error = function(e) NULL)
  }
  if (is.null(L)) {
    # a fit whose covariance did not come back from the Hessian has no
    # law to draw from; the plug-in predictive distribution is still
    # defined and is what the caller gets, said out loud
    frm_warning("predict(): the covariance of the estimates could not be ",
                "recovered from the Hessian, so the draws are taken at the ",
                "estimates alone. The interval carries the observation ",
                "noise and NOT the uncertainty in the estimates. ",
                "param_uncertainty = FALSE asks for this on purpose",
                call. = FALSE)
    return(function(s) object)
  }
  map <- ds$map
  v0 <- fit_outer_vector(object, map)
  p <- length(v0)
  # all ndraws draws at once, column s for replicate s: the caller's
  # stream is consumed by the same amount whatever the rows are
  D <- v0 + crossprod(L, matrix(stats::rnorm(p * ndraws), p, ndraws))
  function(s) fit_set_outer(object, D[, s], map)
}
