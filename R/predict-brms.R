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
# and the uncertainty in the estimates. `propagate_error = FALSE`
# drops the second and simulates at the estimates alone, which is the
# plug-in predictive distribution; `dev/shapes-findings.md` reports what
# each one covers.
#
# The group effects at the levels the fit saw are DRAWN too, jointly,
# from their conditional law given the data and the replicate's
# parameters (predict_b_drawer()). brms's draws carry the posterior of
# each seen level's effect, and holding them at their modes made the
# interval at a known level narrow: 0.9445 and 0.9336 out of sample on
# two designs (dev/shapes-findings.md section 3). The user decided on
# 2026-09-22 to match brms; dev/reunc-findings.md has the construction
# and the coverage it claims.

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
#' All three sources.
#'
#' * Observation noise, from the family's simulator. This is almost
#'   always the largest term.
#' * Uncertainty in the estimates, from their asymptotic covariance.
#' * Uncertainty in the group effects of a level the fit SAW. Each
#'   replicate draws the whole vector of group effects from its
#'   conditional law given the data and that replicate's parameters,
#'   which the joint precision of the fit describes. For a linear mixed
#'   model this is the prediction error variance of the best linear
#'   unbiased predictor, jointly with the fixed effects. brms's draws
#'   carry the posterior of each seen level's effect, and this is the
#'   frequentist analogue. Two rows of one group share that group's
#'   draw, so `summary = FALSE` is right across rows as well as per
#'   row.
#'
#' `propagate_error = FALSE` drops the second and the third together
#' and simulates at the estimates alone, which leaves the observation
#' noise. Everything estimated is then treated as known.
#'
#' What the interval covers: a new observation at a known level, with
#' the nominal coverage AVERAGED over the groups. It is not a coverage
#' statement for one group's realized effect. The predictor shrinks a
#' group toward the population, so for a group far from the population
#' the interval covers less than nominal and for a group near it more.
#' `dev/reunc-findings.md` reports the measured coverage.
#'
#' A level the fit did NOT see is different: there is no conditional
#' law to draw from, so `allow_new_levels = TRUE` draws its effect from
#' the block's own estimated covariance, once per replicate, at that
#' replicate's parameters. This is brms's `sample_new_levels =
#' "gaussian"` and it is the whole between-group variance, so the
#' interval at an unseen level is wider than at a known one.
#'
#' @section re_formula and propagate_error are different questions:
#' The two are easy to confuse and neither can express the other.
#'
#' * `re_formula` chooses WHICH terms enter the prediction. Only
#'   `re_formula = NA` can say "predict for an average group", by
#'   leaving the group effects out of the linear predictor.
#' * `propagate_error` chooses whether the error in the ESTIMATES is
#'   propagated into the interval. Only `propagate_error = FALSE` can
#'   say "include this group's own effect but treat it as known".
#'
#' So a prediction for the group in front of you, with its effect held
#' at the fitted value, is `propagate_error = FALSE`, and a prediction
#' for a group you have not met is `re_formula = NA`. They combine:
#' `re_formula = NA` with `propagate_error = FALSE` is the population
#' curve with no estimation error at all.
#'
#' brms has no such argument, because its draws always carry both: a
#' posterior draw has its own parameters and its own group effects, so
#' there is nothing to switch off. A brms user looking for the missing
#' argument is looking for this one.
#'
#' Only the terms `re_formula` keeps are drawn. A population smooth and
#' a `gp()` curve stay at their modes under every setting.
#'
#' A fit whose covariance could not be recovered from the Hessian has no
#' second term to draw: the simulation falls back to the estimates and
#' warns, rather than drawing from a matrix of `NaN`. A fit made with
#' `quadrature = TRUE` has no group effects in its joint precision: its
#' group effects stay at their modes, and a warning says so.
#'
#' Under `REML = TRUE` (or `frmtmb_control(profile = TRUE)`) the fixed
#' effects are integrated out of the outer problem, so they are not in
#' `vcov(object, full = TRUE)`. They are drawn from their joint
#' covariance with the outer parameters, which the joint precision
#' gives, and the group effects are then drawn given both.
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
#' @param re_formula Which group-level terms enter the prediction.
#'   `NULL` (default) keeps all of them, `NA` keeps none. A one-sided
#'   formula keeps the terms it names and drops the others, as in brms:
#'   `~ (1 | g)` on a fit with `(1 + x | g) + (1 | h)` keeps the
#'   intercept of `g` alone. A formula that names no group-level term,
#'   such as `~0` or `~1`, is `NA`. A term the fit does not have is an
#'   error that names it; brms drops such a term silently. See
#'   [frm_linpred()] for the full rule.
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
#' @param propagate_error Whether the error in the estimates is
#'   propagated into the interval. `TRUE`, the default, draws the
#'   parameters and the group effects of a level the fit saw. `FALSE`
#'   holds both at their estimates, so the interval carries the
#'   observation noise alone: the plug-in predictive distribution. It
#'   is the only way to include a group's own effect and still treat it
#'   as known; `re_formula` cannot say that, because it chooses which
#'   terms are in the prediction and not whether their error is
#'   carried.
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
                               propagate_error = TRUE) {
  frm_check_dots(..., .unsupported = predict_type_retired)
  require_fitted(object, "predict()")
  check_flag(summary, "summary")
  check_flag(robust, "robust")
  check_flag(negative_rt, "negative_rt")
  check_flag(allow_new_levels, "allow_new_levels")
  check_flag(propagate_error, "propagate_error")
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
  # once, here: every replicate then predicts on the reduced design, and
  # the group-effect draw sees which blocks are left in it
  rr <- re_resolve(object, re_formula, "predict()")
  object <- rr$fit
  re_formula <- rr$re_formula
  ds <- predict_simulate(object, rspecs, newdata, re_formula,
                         allow_new_levels, ndraws, ntrys, propagate_error)
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
                             propagate_error) {
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
    ctx0 <- sim_context(object, rspec, list(), aterms = av[[nm]])
    # a cov = FALSE term is drawn row by row around the one-step mean,
    # which frm_linpred() forms on newdata from newdata's own response
    if (autocor_is_cond(ctx0[["autocor"]])) ctx0[["autocor"]] <- NULL
    if (!is.null(newdata) && sim_is_structured(ctx0)) {
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
  draw <- predict_par_drawer(object, propagate_error, ndraws)
  seeds <- sample.int(.Machine$integer.max, ndraws)
  # the stream is captured HERE, before the group-effect draw takes its
  # own seeds, so that draw costs the caller nothing: a conditional
  # predict() leaves .Random.seed exactly where the parameter draws and
  # the simulation seeds left it, which is where 0.61.0 left it. Taking
  # them after the capture was worth 18 differing quantities in a
  # reviewer's rebuild, on rows whose PREDICTIONS were identical.
  saved <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", saved, envir = globalenv()), add = TRUE)
  bdraw <- predict_b_drawer(object, re_formula, attr(draw, "delta"),
                            ndraws, propagate_error)
  joint <- isTRUE(object$spec[["rescor"]]) && length(rspecs) > 1L
  out <- stats::setNames(vector("list", length(rspecs)), names(rspecs))
  masked <- stats::setNames(vector("list", length(rspecs)), names(rspecs))
  for (s in seq_len(ndraws)) {
    fs <- draw(s)
    if (!is.null(bdraw)) fs <- bdraw(fs, s)
    # after the group-effect draw, which has a stream of its own; draw()
    # takes nothing from the stream, so this is where it always was
    set.seed(seeds[s])
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
  # a cov = FALSE term is already in dpv's mu (predict_dpar_values()),
  # conditional on the observed past, so the rows are drawn one by one
  # rather than by the recursion simulate() runs over its own draws
  rowwise <- autocor_is_cond(fs$frame[["autocor"]][[rspec$resp_name]])
  if (!all(ok)) {
    ctx0 <- sim_context(fs, rspec, list(), aterms = av)
    if (rowwise) ctx0[["autocor"]] <- NULL
    if (sim_is_structured(ctx0)) return(NULL)
    dpv <- lapply(dpv, subset_rows, keep = ok)
    av <- lapply(av %||% list(), subset_rows, keep = ok)
  }
  ctx <- sim_context(fs, rspec, dpv, aterms = av, n = sum(ok),
                     extra = fit_extras(fs), max_iter = ntrys)
  if (rowwise) ctx[["autocor"]] <- NULL
  sim_draw(ctx)
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
             "propagate_error = FALSE simulates at the estimates ",
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
                "vcov(object, full = TRUE) shows. propagate_error = ",
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
    # brms's posterior_predict() under cov = FALSE draws around the
    # one-step mean, which conditions on the observed earlier residuals;
    # the frm_linpred() route below already carries it
    dpv <- autocor_cond_dpars(fit, resp, eval_dpars(fit)[[resp]])
    return(cs_offsets_add(fit, resp, NULL, dpv))
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
  # with_cs_offsets() takes the list of EVERY response and reads the
  # cs() values of the training rows; handed this one response's list
  # it wrote the offsets one level down, where no simulator reads them,
  # and predict() drew every cs() model as if the term were absent
  cs_offsets_add(fit, resp, newdata, dpv)
}

#' A closure returning one fit-like object per call: the fit itself
#' under `propagate_error = FALSE`, and otherwise the fit with its
#' outer parameters replaced by a draw from
#' `N(theta_hat, vcov(full = TRUE))`.
#'
#' The Cholesky factor is taken once, outside the loop, because it is
#' the only expensive part of the draw.
#'
#' @noRd
predict_par_drawer <- function(object, propagate_error, ndraws) {
  if (!propagate_error) return(function(s) object)
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
                "propagate_error = FALSE asks for this on purpose",
                call. = FALSE)
    return(function(s) object)
  }
  map <- ds$map
  v0 <- fit_outer_vector(object, map)
  p <- length(v0)
  # all ndraws draws at once, column s for replicate s: the caller's
  # stream is consumed by the same amount whatever the rows are
  D <- v0 + crossprod(L, matrix(stats::rnorm(p * ndraws), p, ndraws))
  out <- function(s) fit_set_outer(object, D[, s], map)
  # the group-effect draw conditions on this replicate's parameters, so
  # it needs how far they moved (predict_b_drawer())
  attr(out, "delta") <- function(s) D[, s] - v0
  out
}

#' The replicate's group effects at the levels the fit SAW, drawn from
#' their conditional law given the data and the replicate's parameters.
#'
#' brms's draws carry the posterior of each seen level's effect `r_g`,
#' so its predictive interval at a known level includes the uncertainty
#' in that effect. The frequentist analogue is the conditional law of
#' `b` given the data, which is what the joint precision `Q` of
#' `sdreport()` describes: its `b` block is the Hessian of the inner
#' problem, and its off-diagonal block says how the modes move with the
#' parameters. Partitioning `Q` into the drawn parameters `d` (the
#' outer vector, plus `beta` under REML or profile) and the rest `r`,
#'
#'     r | d  ~  N(r_hat - Q_rr^-1 Q_rd (d - d_hat),  Q_rr^-1),
#'
#' which for a linear mixed model is Henderson's prediction error
#' variance of the BLUP, jointly with the fixed effects. `Q` is read
#' ONCE, at the estimates: the replicate's own parameters shift the
#' conditional MEAN, through `d - d_hat`, and not the conditional
#' variance, which would need the joint precision re-evaluated at every
#' draw. The draw is of the WHOLE vector, so
#' two rows of one group share that group's draw and rows of different
#' groups carry exactly the correlation `Q` implies. A per-row added
#' variance would get the summary right and `summary = FALSE` wrong.
#'
#' Only the blocks `re_formula = NA` would drop are written back: a
#' population smooth or a `gp()` curve is held at its mode as before.
#'
#' `propagate_error = FALSE` turns this draw OFF along with the
#' parameter draw. The two are one axis and the argument names it:
#' with it `FALSE` every estimated quantity is held at its estimate,
#' the group effects included, and the interval carries the
#' observation noise alone. That is the only way to ask for a group's
#' own effect while treating it as KNOWN, which `re_formula` cannot
#' express, because `re_formula` chooses WHICH terms are in the
#' prediction and not whether their error is propagated.
#'
#' Each replicate draws from its own seed, taken up front AFTER the
#' parameter draws and the simulation seeds, so those two are exactly
#' what they were before this draw existed and the rows at an unseen
#' level come out bitwise as they did.
#'
#' Returns `NULL`, and takes nothing from the caller's stream, when
#' `propagate_error` is `FALSE` or no kept block has a seen level to
#' draw.
#'
#' @noRd
predict_b_drawer <- function(object, re_formula, delta, ndraws,
                             propagate_error) {
  if (!propagate_error) return(NULL)
  if (!re_form_keeps(re_formula)) return(NULL)
  gov <- re_governed_b(object)
  if (!length(gov)) return(NULL)
  Q <- joint_precision(object)
  rn <- rownames(Q)
  if (is.null(Q) || !any(rn == "b")) {
    frm_warning("predict(): the fitted objective carries no random-effect ",
                "block in its joint precision (quadrature = TRUE ",
                "marginalizes it), so the group effects at the levels the ",
                "fit saw are held at their modes and their uncertainty is ",
                "NOT in the interval", call. = FALSE)
    return(NULL)
  }
  map <- fit_draw_space(object)$map
  pos_d <- unlist(lapply(unique(map$comp), function(cp) which(rn == cp)))
  pos_r <- setdiff(seq_len(nrow(Q)), pos_d)
  b_r <- which(rn[pos_r] == "b")
  nb <- length(object$estimates[["b"]])
  if (length(pos_d) != length(map$names) || length(b_r) != nb) {
    frm_warning("predict(): the joint precision does not line up with ",
                "the fit's parameters (", length(pos_d), " drawn rows for ",
                length(map$names), " parameters, ", length(b_r),
                " b rows for ", nb, " effects), so the group effects at the ",
                "levels the fit saw are held at their modes and their ",
                "uncertainty is NOT in the interval. Please report this ",
                "model", call. = FALSE)
    return(NULL)
  }
  Qrr <- methods::as(Matrix::forceSymmetric(Q[pos_r, pos_r, drop = FALSE]),
                     "CsparseMatrix")
  ch <- tryCatch(Matrix::Cholesky(Qrr, perm = TRUE, LDL = FALSE),
                 error = function(e) NULL)
  if (is.null(ch)) {
    frm_warning("predict(): the conditional precision of the random ",
                "effects is not positive definite, so the group effects at ",
                "the levels the fit saw are held at their modes and their ",
                "uncertainty is NOT in the interval. diagnose() names the ",
                "likely cause", call. = FALSE)
    return(NULL)
  }
  Qrd <- Q[pos_r, pos_d, drop = FALSE]
  sel <- b_r[gov]
  b0 <- object$estimates[["b"]][gov]
  nr <- length(pos_r)
  bseeds <- sample.int(.Machine$integer.max, ndraws)
  function(fs, s) {
    set.seed(bseeds[s])
    x <- Matrix::solve(ch, Matrix::solve(ch, stats::rnorm(nr), system = "Lt"),
                       system = "Pt")
    x <- as.numeric(x)
    if (!is.null(delta)) {
      x <- x - as.numeric(Matrix::solve(ch, Qrd %*% delta(s),
                                        system = "A"))
    }
    fs$estimates[["b"]][gov] <- b0 + x[sel]
    fs
  }
}

#' Positions in `b` of the group effects the ROWS being predicted
#' actually load, which is what a finite-difference delta method has to
#' perturb.
#'
#' `re_governed_b()` names every level of every kept block, and a
#' prediction of a few rows loads a few of them; the rest have a
#' derivative of exactly zero and cost one model evaluation each. On a
#' cumulative fit with 200 grouping levels, differencing all of them
#' took 10.92 s against 0.31 s before the group effects joined the
#' difference. Bounding the set by the rows makes the cost grow with the
#' rows rather than with the number of levels.
#'
#' The design that multiplies the coefficient vector is the source, so
#' this is the same set `re_eta()` and `lp_delta_A()` read: in sample
#' each predictor's `Z`, and on new data the rebuilt
#' `re_design_matrix()` plus the smooth parts. A column with a nonzero
#' entry is loaded; every other column contributes exactly nothing.
#'
#' Returns `NULL` when the design cannot be rebuilt, which means "no
#' bound": the caller then perturbs every kept level, as it did before.
#'
#' @noRd
re_used_b <- function(fit, newdata, resp, allow_new_levels) {
  blocks <- fit$frame[["re_blocks"]]
  if (!length(blocks)) return(integer(0))
  cols <- integer(0)
  n_c <- fit$frame[["n_c"]] %||% length(fit$estimates[["b"]])
  for (lp in fit$frame[["linpreds"]]) {
    if (!is.null(resp) && !identical(lp[["resp"]], resp)) next
    if (is.null(lp[["Z"]])) next
    if (is.null(newdata)) {
      cols <- c(cols, which(Matrix::colSums(abs(lp[["Z"]])) > 0))
      next
    }
    ed <- tryCatch(suppressWarnings(
      lp_eta_design(fit, lp, newdata, TRUE, allow_new_levels)),
      error = function(e) NULL)
    if (is.null(ed)) return(NULL)
    if (length(ed[["re_parts"]])) {
      Zn <- re_design_matrix(ed[["re_parts"]], ed[["n"]], n_c)
      cols <- c(cols, which(Matrix::colSums(abs(Zn)) > 0))
    }
    for (sp in ed[["sm_parts"]]) {
      used <- which(Matrix::colSums(abs(as.matrix(sp$Xr))) > 0)
      cols <- c(cols, sp$bk[["c_idx"]][used])
    }
  }
  cols <- unique(cols)
  out <- integer(0)
  for (bk in blocks) {
    ci <- bk[["c_idx"]]
    bi <- bk[["b_idx"]]
    hit <- which(ci %in% cols)
    if (!length(hit)) next
    # A loaded COLUMN names one b position only where expand_b() is
    # positionwise. An rr block's coefficient is that level's factors
    # through the loadings, and an esicar block's is b minus its
    # component's mean, so there one b entry reaches every coefficient
    # of its component and a bound would drop Jacobian columns that are
    # not zero. The whole block joins instead.
    out <- c(out, if (block_b_positionwise(bk)) bi[hit] else bi)
  }
  sort(unique(out))
}

#' The rows each candidate group effect reaches, as one sparse
#' indicator per linear predictor of `resp`, summed.
#'
#' It must cover EVERY column `re_used_b()` can return, because a
#' column with no support here looks to `re_b_batches()` like a
#' coefficient no row loads, and the batch then writes it a derivative
#' of exactly zero. In sample the predictor's `Z` carries every block,
#' including a factor smooth's basis. On new data the design is rebuilt
#' in two pieces, and reading only `re_parts` left every factor-smooth
#' column unsupported: `fitted(newdata = )` on a cumulative fit with
#' `s(x, g, bs = "fs")` AND an ordinary group term put `Est.Error` 68
#' percent wrong, silently, because the ordinary term gave the batch
#' something to attribute while the smooth gave it nothing
#' (`dev/reunc-log/fdsmooth-prefix.txt`). The error is not one-signed:
#' the variance carries `2 Jd' Vdb Jb`, so zeroing a `Jb` column can
#' move a cell either way. The smooth parts join it here.
#'
#' @noRd
re_row_support <- function(fit, newdata, resp, allow_new_levels) {
  n_c <- fit$frame[["n_c"]] %||% length(fit$estimates[["b"]])
  out <- NULL
  for (lp in fit$frame[["linpreds"]]) {
    if (!is.null(resp) && !identical(lp[["resp"]], resp)) next
    if (is.null(lp[["Z"]])) next
    M <- if (is.null(newdata)) {
      lp[["Z"]]
    } else {
      ed <- tryCatch(suppressWarnings(
        lp_eta_design(fit, lp, newdata, TRUE, allow_new_levels)),
        error = function(e) NULL)
      if (is.null(ed)) return(NULL)
      if (!length(ed[["re_parts"]]) && !length(ed[["sm_parts"]])) next
      Zn <- if (length(ed[["re_parts"]])) {
        re_design_matrix(ed[["re_parts"]], ed[["n"]], n_c)
      } else {
        Matrix::sparseMatrix(i = integer(0), j = integer(0),
                             x = numeric(0), dims = c(ed[["n"]], n_c))
      }
      for (sp in ed[["sm_parts"]]) {
        nz <- which(abs(as.matrix(sp$Xr)) > 0, arr.ind = TRUE)
        if (!nrow(nz)) next
        Zn <- Zn + Matrix::sparseMatrix(
          i = nz[, 1L], j = sp$bk[["c_idx"]][nz[, 2L]], x = 1,
          dims = c(ed[["n"]], n_c))
      }
      Zn
    }
    M <- methods::as(abs(M) > 0, "dMatrix")
    out <- if (is.null(out)) M else out + M
  }
  out
}

#' Group the perturbed effects into batches that can be differenced
#' TOGETHER, exactly.
#'
#' A finite difference over `b` costs one pair of model evaluations per
#' coefficient, and in sample every level of every block is loaded, so
#' the bound `re_used_b()` gives cannot help there: a cumulative fit
#' with 1000 levels took 142 s.
#'
#' It does not have to cost that. Within ONE block a row loads at most
#' one level, so perturbing every level of that block at once changes
#' each row by exactly its own level's perturbation, and the difference
#' read off that row is exactly the difference the single-coefficient
#' perturbation would have produced. The other entries of that column
#' are exactly zero either way. So a block costs ONE pair of
#' evaluations rather than one per level, and the answer is the same
#' number, not an approximation of it.
#'
#' The condition is checked, not assumed: a row that loads two levels of
#' one block (a multi-membership term) breaks it, and so does a block
#' whose coefficients are a function of several of its `b` entries (an
#' `rr` block). Those fall back to one coefficient at a time.
#'
#' Returns a list of batches, each `list(idx, owner)`, where `owner[r]`
#' is the position in `idx` of the effect that row `r` loads, or `NA`;
#' or `NULL` when the design could not be built, which means "no
#' batching".
#'
#' @noRd
re_b_batches <- function(fit, newdata, resp, allow_new_levels, b_idx) {
  if (!length(b_idx)) return(list())
  S <- re_row_support(fit, newdata, resp, allow_new_levels)
  if (is.null(S)) return(NULL)
  out <- list()
  for (bk in fit$frame[["re_blocks"]]) {
    ci <- bk[["c_idx"]]
    bi <- bk[["b_idx"]]
    if (!any(bi %in% b_idx)) next
    # The batch assumes expand_b() carries b to the coefficients
    # POSITIONWISE, so that perturbing every level at once moves each
    # row by its own level's step and nothing else. That is the whole
    # premise and it is asked as itself, not through a proxy: the
    # earlier test was `length(c_idx) != length(b_idx)`, which catches
    # rr below full rank and misses rr AT full rank and esicar, where
    # the lengths agree and `car_center()` subtracts the component
    # mean, so every row moved by its own step MINUS the mean of all of
    # them. No batching at all for such a block.
    if (!block_b_positionwise(bk)) return(NULL)
    # a positionwise block has one b entry per coefficient by
    # definition, so this cannot fire; it is kept as an indexing guard
    # because everything below indexes b_idx through c_idx's length
    if (length(ci) != length(bi)) return(NULL)
    D <- max(1L, bk[["dim"]])
    # one batch per COLUMN POSITION of the block, not per block: a row
    # of a `(1 + x | g)` term loads its level's intercept AND its slope,
    # so a batch over the whole block could not attribute the
    # difference, while a batch over one position holds one nonzero per
    # row. An |ID| block spanning two predictors splits the same way,
    # because each predictor is its own position.
    for (k in seq_len(D)) {
      at <- seq.int(k, length(ci), by = D)
      keep <- at[bi[at] %in% b_idx]
      if (!length(keep)) next
      Sb <- S[, ci[keep], drop = FALSE]
      # a row that loads two levels at this position cannot be
      # attributed (a multi-membership term does exactly that)
      if (max(Matrix::rowSums(Sb)) > 1) return(NULL)
      owner <- rep(NA_integer_, nrow(Sb))
      nz <- Matrix::which(Sb != 0, arr.ind = TRUE)
      owner[nz[, 1L]] <- nz[, 2L]
      out[[length(out) + 1L]] <- list(idx = bi[keep], owner = owner)
    }
  }
  out
}

#' Positions in `b` of the blocks `re_formula` governs: every block with
#' a component left in the (possibly reduced) design, and a smooth whose
#' basis is indexed by a grouping factor. A population smooth, `gp()`
#' and `hsgp()` stay out, because `re_formula = NA` keeps them too.
#'
#' @noRd
re_governed_b <- function(fit) {
  blocks <- fit$frame[["re_blocks"]]
  gs <- unlist(lapply(fit$frame[["linpreds"]], smooth_group_block_ids))
  ids <- integer(0)
  for (i in seq_along(blocks)) {
    bk <- blocks[[i]]
    if (bk[["covstruct"]] %in% c("smooth", "gp", "hsgp")) {
      if (i %in% gs) ids <- c(ids, bk[["b_idx"]])
      next
    }
    if (length(bk[["components"]])) ids <- c(ids, bk[["b_idx"]])
  }
  sort(unique(ids))
}
