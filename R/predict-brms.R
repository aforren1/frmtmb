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
#' A new row of a `me()` term also has no latent value, so each replicate
#' draws one from `normal(x, sdx)`, per row, or per level under `me(gr = )`, as
#' brms's `posterior_predict()` does on new data. The latent value of a new row
#' is
#' measurement noise in the grid, not an estimate, so `propagate_error = FALSE`
#' keeps it. `corme` is not read, as in
#' brms. Two responses or two dpars that read one `me()` term read ONE latent
#' value, a stated divergence from brms, which draws one per predictor against
#' its
#' own Stan program: at the training rows of a two-response fit the residuals
#' correlate at 0.734, one latent predicts 0.739 and brms's replicates give
#' -0.001, and the 95 percent width of `y1 - y2` at `sdx = 2` is 3.18 with one
#' latent and 15.08 with two. `newdata` must hold `sdx`. On the data the model
#' was fitted on the latent values are the fit's own.
#'
#' @section When the normal law does not describe the likelihood:
#' The normal law reads the curvature of the log-likelihood at the
#' estimates. Three checks read the log-likelihood itself.
#'
#' * Each side of each parameter, moved alone. Where the log-likelihood
#'   is flat toward one end of a parameter's range (a zero inflation
#'   with no excess zeros in the data), the draws do not move it that
#'   way. Where it falls much faster than the normal law says, the draws
#'   on that side are narrowed. A message of class `frmtmb_draw_held`
#'   says so, once per fit.
#' * The draws together. The log-likelihood is evaluated at 20 test
#'   draws (40 where the result is close to the threshold) and compared,
#'   draw by draw, with what the normal law says it is there. This is
#'   done once per fit. The check fires where the median of the excess
#'   drop is above 1 and above 0.15 times the median drop the normal law
#'   gives at those draws. A weakly identified mixture can give that
#'   with no parameter at an end. A model with many well-identified
#'   parameters does not: a small skew in each adds up to an excess that
#'   grows with their number, and so does the law's own drop.
#' * Each direction of the normal law. The median above is a sum over
#'   all directions, so one weak direction among many well-identified
#'   parameters, such as a mixture with many covariates, does not show
#'   in it. Where the draws together do not fire the check, each
#'   parameter is moved 2 of its standard deviations to each side, and
#'   the parameters correlated with it are moved by their regression on
#'   it. The excess of the drop there over the normal law's is read, and
#'   the excess of the parameter moved alone, by 2 of its conditional
#'   standard deviations (nearer the estimate), is taken off, so that a
#'   group standard deviation that the first check judges does not count
#'   twice. The check fires where what is left is above 6. That
#'   threshold is a heuristic: for a quadratic misfit it is a normal law
#'   twice as wide as the likelihood along that direction. A skew in a
#'   strongly correlated parameter alone can still count. Where the
#'   log-likelihood cannot be evaluated at a direction's end, the
#'   direction is read nearer the estimate, or not at all. This is done
#'   once per fit and costs two evaluations of the log-likelihood per
#'   parameter; on a model with random effects each is an inner
#'   optimization. A random-slope model whose slope variance and
#'   correlation are weakly identified together can fire it: there the
#'   normal law's far tail for a new level is too wide. Where every
#'   direction that fires is a covariance-structure parameter's, the
#'   draws are NOT reweighted. These are the parameters of every
#'   group-level covariance (standard deviations and correlations, and
#'   the ranges, Matern smoothness, decay rates, ar correlations and
#'   `rr()` loadings of structured covariances, spatial ones included),
#'   the `gp()` and smooth scales, the latent-residual standard
#'   deviation and its autocorrelation, and the autocorrelation,
#'   residual-correlation and `me()` parameters. The coefficients of
#'   every linear predictor, a sigma's among them, are not. On the
#'   random-slope fits that fire, the likelihood is flat toward a slope
#'   variance of 0, and weights by it made predictions worse on 10 of
#'   12. A warning of class `frmtmb_draw_wide_tail` names the parameter
#'   and says that the normal law's far tail along it is too wide. Where
#'   another direction fires too, the draws are reweighted. A weakness
#'   spread over many directions, each only a little off the normal law,
#'   can still go unflagged and keep the normal law's draws.
#'
#' Where the second check fires, or the third on any other direction,
#' each draw is weighted by the likelihood over the density it was
#' drawn from. The weights are Pareto smoothed with `loo::psis()` and
#' the draws are resampled by them. A parameter that the first check
#' holds at its estimate in part of the draws keeps that part. A
#' message of class
#' `frmtmb_draw_reweighted` gives the Pareto k. Where k is above 0.7 the
#' weights are not reliable and a warning of class `frmtmb_draw_pareto`
#' says so: `frmtmb.sample::frm_sample()` draws from the posterior
#' itself, though its `posterior_predict()` cannot yet predict a
#' grouping level the fit did not see, and the warning says so where the
#' prediction is for one. Whether k passes 0.7 depends on `ndraws`.
#' Neither check uses the random number stream, and where neither fires
#' the draws are not changed.
#'
#' @section predict() and fitted() can disagree on the mean:
#' [fitted()] is the expected response at the estimates. The mean of
#' `predict()`'s draws is the expected response AVERAGED over the
#' parameter draws, and where the response is not linear in the
#' parameters the two differ by more than Monte Carlo error. On a
#' hurdle_negbinomial fit of 300 rows the difference was 5 Monte Carlo
#' standard errors at 1000 draws, under 1 percent of the mean.
#' `propagate_error = FALSE` gives draws whose mean is `fitted()`'s. brms
#' has no such difference, because its `fitted()` is itself a mean over
#' the draws.
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
#' Only the terms `re_formula` keeps are drawn. A `gp()` term is drawn
#' under every setting, as brms's draws carry it: its coefficients given
#' the drawn parameters, and at a position the fit did not see, the field
#' there given the drawn field. A population smooth stays at its mode. An
#' exact field takes each draw's own inner mode, by Newton steps on numbers
#' and one factor of its kernel per draw: on a gaussian field at 500
#' positions and a poisson one at 300, a draw took 0.9 to about 2 times as
#' long as one that held the field at its mode. Beside another random effect
#' or in a nonlinear body the inner problem is solved on the tape instead, 3
#' to 21 times as long, and so is a draw whose Newton steps do not settle or
#' meet a curvature that is not positive: on binomial fields with many trials
#' 20 to 38 percent of the draws do, and
#' those cost what they did before.
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
#' @param sample_new_levels brms's argument. `"gaussian"` (or `NULL`)
#'   draws an unseen level's effect from its block's estimated
#'   covariance. `"old_levels"` gives each unseen level the effects of
#'   one seen level of its grouping factor, chosen at random once per
#'   call as brms chooses it, and every term of that factor reads the
#'   same seen level (within the row's by-level for `gr(by = )`),
#'   from each replicate's draw of the group effects. One choice serves
#'   every response of a multivariate call. Two different unseen values
#'   of one `mm()` term read two seen levels; brms 2.23.0 gives them
#'   one, because it numbers unseen values per member rather than by
#'   value. `"uncertainty"`
#'   mixes the two per posterior draw, and a maximum-likelihood fit has
#'   none, so it is refused by name.
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
  require_structure_packages(object)
  # its draws are taken at sigma's value
  autocor_check_fit(object)
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
  sample_new_levels <- check_sample_new_levels(
    sample_new_levels, "predict()", allow = c("gaussian", "old_levels"))
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
  resp <- resp_read(object, resp)
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
  # subset(): one response at a time, on newdata's rows where its subset
  # is TRUE, as brms asks (R/subset.R)
  subset_resp_check(object, resp, "predict()")
  newdata <- subset_newdata(object, resp, newdata)
  if (length(rspecs) > 1L) predict_mv_refuse(object, rspecs)
  # once, here: every replicate then predicts on the reduced design, and
  # the group-effect draw sees which blocks are left in it
  rr <- re_resolve(object, re_formula, "predict()")
  object <- rr$fit
  re_formula <- rr$re_formula
  na_said <- if (summary) {
    names(rspecs)[vapply(rspecs, function(rs) {
      ord_response_has_cs(object, rs)
    }, TRUE)]
  }
  ds <- predict_simulate(object, rspecs, newdata, re_formula,
                         allow_new_levels, ndraws, ntrys, propagate_error,
                         sample_new_levels, na_said = na_said %||% character(0))
  if (!is.null(transform)) ds <- lapply(ds, match.fun(transform))
  if (!summary) {
    return(if (length(ds) == 1L) pp_dimnames(ds[[1L]]) else
      predict_stack(ds))
  }
  per <- lapply(seq_along(ds), function(i) {
    cross <- predict_cross_rows(object, rspecs[[i]], newdata, re_formula,
                                allow_new_levels)
    predict_summarize_one(object, rspecs[[i]], ds[[i]], probs, robust,
                          cross)
  })
  if (length(per) == 1L) return(per[[1L]])
  out <- array(unlist(per), c(nrow(per[[1L]]), ncol(per[[1L]]), length(per)),
               dimnames = list(rownames(per[[1L]]), colnames(per[[1L]]),
                               names(ds)))
  out
}

#' brms's `sample_new_levels`, as far as a maximum-likelihood fit can
#' honor it, returned as the value in force. `fitted()` and `predict()`
#' both take the unseen level's effect from its block's estimated
#' covariance, which is brms's `"gaussian"`. Both also take
#' `"old_levels"`: `predict()` reads the chosen seen level's effect from
#' each replicate's draw, and `fitted()` reads that level's fitted
#' effect and conditional variance through the design. `"uncertainty"`
#' chooses per posterior draw, and is refused by name rather than
#' ignored.
#'
#' @noRd
check_sample_new_levels <- function(x, what, allow = "gaussian") {
  if (is.null(x)) return("gaussian")
  if (is.character(x) && length(x) == 1L && x %in% allow) return(x)
  frm_stop(what, " honors sample_new_levels = ",
           paste0("\"", allow, "\"", collapse = " and "), " only. ",
           "An unseen level's effect is taken from its block's own ",
           "estimated covariance, which is what brms's \"gaussian\" ",
           "does",
           if ("old_levels" %in% allow) {
             ", or from one seen level's draw, which is \"old_levels\""
           },
           ". ",
           if (!"old_levels" %in% allow) {
             paste0("\"old_levels\" resamples the draws of the levels ",
                    "that were seen, and this summary has none. ")
           },
           "\"uncertainty\" mixes the two per POSTERIOR draw, and a ",
           "maximum-likelihood fit has none", call. = FALSE)
}

#' The rows of an ordinal `cs()` response whose thresholds cross at the
#' plug-in estimate, where `fitted()` is `NaN`, or `NULL` where the
#' response cannot cross (`ord_response_has_cs()`). A row that cannot
#' be estimated is `NA` there, not `NaN`, and is not counted.
#'
#' @noRd
predict_cross_rows <- function(object, rspec, newdata, re_formula,
                               allow_new_levels) {
  if (!ord_response_has_cs(object, rspec)) return(NULL)
  P <- tryCatch(frm_linpred(object, newdata = newdata,
                            resp = rspec$resp_name, type = "response",
                            re_formula = re_formula,
                            allow_new_levels = allow_new_levels),
                error = keep_missing(NULL))
  if (!is.matrix(P)) return(NULL)
  rowSums(is.nan(P)) > 0L
}

#' brms's summary of one response's simulated draws.
#'
#' @noRd
predict_summarize_one <- function(object, rspec, d, probs, robust,
                                  cross = NULL) {
  if (fam_is_category_valued(rspec$family)) {
    return(predict_category_props(object, rspec, d, cross))
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

#' A one-response draws matrix with brms's dimnames, `list(NULL, NULL)`
#' rather than none, as brms's `posterior_predict()`, `predict(summary =
#' FALSE)`, `predictive_error()` and `log_lik()` return it (brms 2.23.0,
#' dev/names-findings.md). Any other shape is left alone.
#'
#' @noRd
pp_dimnames <- function(m) {
  if (is.matrix(m) && is.null(dimnames(m))) dimnames(m) <- list(NULL, NULL)
  m
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

#' Whether a response has a `cs()` term on a predictor whose family
#' keeps its thresholds ordered (`cumulative()`, `hurdle_cumulative()`,
#' or such a component of an ordinal mixture): the case where a row's
#' thresholds can cross.
#'
#' @noRd
ord_response_has_cs <- function(object, rspec) {
  if (!identical(rspec$family[["type"]], "ordinal")) return(FALSE)
  for (dp in names(rspec$dpars)) {
    lp <- object$frame[["linpreds"]][[linpred_key(rspec$resp_name, dp)]]
    if (!length(lp[["cs"]])) next
    tf <- cs_target_family(rspec$family, dp)
    if (!is.null(tf) && tf[["family"]] %in% ord_ordered_families) {
      return(TRUE)
    }
  }
  FALSE
}

#' brms's predictive summary of a category-valued response: the
#' simulated proportion of each category, one column per category.
#'
#' `cross`, when given, marks the rows whose thresholds cross at the
#' plug-in estimate, where `fitted()` is `NaN`: those rows are `NA`.
#' Under `cs()` on a family that keeps its thresholds ordered, a
#' replicate is `NA` where its parameter draw crosses two of the row's
#' thresholds; elsewhere the proportions are over the replicates that do
#' not cross, and the call warns once with how many were dropped (the
#' user's decision of 2026-10-06).
#'
#' @noRd
predict_category_props <- function(object, rspec, d, cross = NULL) {
  lv <- object$frame[["y_levels"]][[rspec$resp_name]]
  K <- if (!is.null(lv)) {
    length(lv)
  } else if (identical(rspec$family[["type"]], "ordinal")) {
    ordinal_ncat(object, rspec$resp_name)
  } else {
    max(c(as.integer(d), object$frame[["y"]][[rspec$resp_name]]))
  }
  # the drawn codes start at 0 for a hurdle family, whose category 0 is
  # the first column
  code0 <- ord_code0(rspec$family)
  out <- t(apply(d, 2L, function(col) {
    # over the replicates that are THERE: a replicate that drew a
    # non-finite dpar, or whose cs() thresholds cross, is NA, and
    # dividing by the full count would report proportions that do not
    # sum to one
    ok <- col[!is.na(col)]
    if (!length(ok)) return(rep(NA_real_, K))
    tabulate(as.integer(ok) + 1L - code0, K) / length(ok)
  }))
  if (ord_response_has_cs(object, rspec)) {
    if (is.null(cross) || length(cross) != ncol(d)) {
      cross <- rep(FALSE, ncol(d))
    }
    nas <- colSums(is.na(d))
    out[cross, ] <- NA_real_
    dropped <- which(nas > 0L & !cross)
    if (length(dropped) || any(cross)) {
      per <- if (length(dropped) <= 6L) {
        paste0("row ", dropped, " ", nas[dropped], collapse = ", ")
      } else {
        paste0(min(nas[dropped]), " to ", max(nas[dropped]), " per row")
      }
      frm_warning(
        if (length(dropped)) {
          paste0("On ", length(dropped), " of ", ncol(d), " row(s) the cs() ",
                 "offsets make two thresholds cross in some simulated ",
                 "parameter draws, where the row has no category ",
                 "distribution; their proportions are over the other ",
                 "replicates (replicates dropped, of ", nrow(d), ": ", per,
                 ").")
        },
        if (any(cross)) {
          paste0(if (length(dropped)) " ", sum(cross), " row(s) cross at ",
                 "the estimates themselves, as fitted() is NaN there, ",
                 "and their proportions are NA.")
        },
        call. = FALSE)
    }
  }
  colnames(out) <- brms_category_labels(
    lv %||% as.character(code0 + seq_len(K) - 1L), K)
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
#' `na_said` names the responses whose NA draws the caller reports
#' itself (the `cs()` crossing count of `predict_category_props()`), so
#' that they are not reported twice.
#'
#' @noRd
predict_simulate <- function(object, rspecs, newdata, re_formula,
                             allow_new_levels, ndraws, ntrys,
                             propagate_error,
                             sample_new_levels = "gaussian",
                             na_said = character(0)) {
  av <- list()
  nls <- list()
  # brms's "old_levels" choice is one per call and grouping factor,
  # whichever response a term belongs to
  pick_all <- list()
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
    # every term, as simulate() and fitted() read them: a subset here
    # (trunc(), rate() and thres(gr = ) were the only ones) left the
    # simulator without trials(), se() or dec(), and binomial drew one
    # trial per row where newdata had 1000
    av[[nm]] <- if (is.null(newdata)) {
      object$frame[["aterm_values"]][[rspec$resp_name]]
    } else {
      aterms_for_newdata(rspec, newdata)
    }
    ctx0 <- sim_context(object, rspec, list(), aterms = av[[nm]])
    # a cov = FALSE term is drawn row by row around the one-step mean,
    # which frm_linpred() forms on newdata from newdata's own response
    if (autocor_is_cond(ctx0[["autocor"]])) ctx0[["autocor"]] <- NULL
    # a registry residual term (fcor(), sar()) takes newdata over the
    # rows of its matrix, as brms does, and refuses any other
    if (!is.null(newdata) && autocor_is_seam(ctx0[["autocor"]]) &&
          is.null(fam_sim_ctx(ctx0[["family"]]))) {
      # the replicates draw through the block over newdata's rows,
      # which the fit's own frame then carries into every draw
      object$frame[["autocor"]][[rspec$resp_name]] <-
        autocor_seam_for_newdata(object, ctx0[["autocor"]], rspec, newdata,
                                 "predict()")
      ctx0[["autocor"]] <- NULL
    }
    if (!is.null(newdata) && !is.null(ctx0[["autocor"]]) &&
          is.null(fam_sim_ctx(ctx0[["family"]]))) {
      # the cov = TRUE terms: newdata's own groups and times, the route
      # simulate(newdata = ) takes, which agrees with brms's
      # posterior_predict() (dev/review-predread-explain.R)
      object$frame[["autocor"]][[rspec$resp_name]] <-
        autocor_for_newdata(object, ctx0[["autocor"]], rspec, newdata)
      ctx0[["autocor"]] <- NULL
    }
    if (!is.null(newdata) && sim_is_structured(ctx0)) {
      # said for the model at hand: a hidden state sequence, a latent class
      # and a mixture_mvn() fit each draw through their own simulator
      frm_stop("predict(newdata = ) is not supported for family '",
               fam[["family"]], "': it draws the whole response at once, ",
               "through its own simulator, and a prediction at new rows ",
               "draws one row at a time. Drop newdata to predict the rows ",
               "the model was fitted on", call. = FALSE,
               package = frm_family_package(fam))
    }
    nls[[nm]] <- predict_new_level_spec(object, rspec, newdata, re_formula,
                                        allow_new_levels, sample_new_levels,
                                        pick = pick_all)
    pick_all <- attr(nls[[nm]], "old_pick") %||% pick_all
  }
  # Every random quantity a replicate needs from the caller's stream is
  # taken UP FRONT: the parameter draws and one seed per replicate. A
  # replicate then simulates from its own seed. So the parameter draws
  # do not depend on how many rows newdata has, and a row's draws do not
  # depend on the rows AFTER it; they still shift with the rows before
  # it, which share the replicate's stream. The caller's stream is left
  # where the up-front draws put it.
  # the draw check's warning sends a user to frm_sample(), whose
  # posterior_predict() cannot yet predict a level the fit did not see
  new_lv <- any(vapply(nls, function(x) length(x) > 0L, TRUE))
  draw <- predict_par_drawer(object, propagate_error, ndraws,
                             new_levels = new_lv)
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
  simna <- list()
  for (s in seq_len(ndraws)) {
    fs <- draw(s)
    if (!is.null(bdraw)) fs <- bdraw(fs, s)
    # a new row's me() value is drawn per replicate, at the estimates too:
    # it is measurement noise in the grid, not an estimate (me_newdata_value())
    if (!is.null(newdata)) fs[["me_draw"]] <- new.env(parent = emptyenv())
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
        bad <- rowSums(!is.finite(as.matrix(ys))) > 0L
      } else {
        out[[nm]][s, ok] <- as.numeric(ys)
        bad <- !is.finite(as.numeric(ys))
      }
      # finite parameters the simulator gave no finite draw for
      # (thresholds that overflow, cs() offsets that cross, a process
      # with no end in its window): counted, never silent
      if (is.null(simna[[nm]])) simna[[nm]] <- integer(n)
      simna[[nm]][ok] <- simna[[nm]][ok] + bad
    }
  }
  predict_report_masked(out, masked, ndraws)
  predict_report_sim_na(simna[setdiff(names(simna), na_said)], ndraws)
  out
}

#' Say which rows the simulator left NA or infinite at finite parameters.
#'
#' @noRd
predict_report_sim_na <- function(simna, ndraws) {
  msgs <- character(0)
  for (nm in names(simna)) {
    cnt <- simna[[nm]]
    hit <- which(cnt > 0L)
    if (!length(hit)) next
    what <- paste0(hit, " (", cnt[hit], " of ", ndraws, ")")
    if (length(what) > 6L) {
      what <- c(what[1:6], paste0("and ", length(what) - 6L, " more"))
    }
    msgs <- c(msgs, paste0(
      if (length(simna) > 1L) paste0("response '", nm, "': "),
      if (length(hit) == 1L) "row " else "rows ",
      paste(what, collapse = ", ")))
  }
  if (length(msgs)) {
    frm_warning("predict(): the family's simulator returned no finite ",
                "response from a replicate's parameters for ",
                paste(msgs, collapse = "; "), ". Those cells are NA or ",
                "infinite. The parameters are finite there, but the ",
                "family has no finite draw at them. propagate_error = ",
                "FALSE simulates at the estimates alone", call. = FALSE,
                class = "frmtmb_predict_na")
  }
  invisible(NULL)
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
                     extra = fit_extras(fs, rspec$resp_name),
                     max_iter = ntrys)
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
#' The joint law is the multivariate normal, or for a Student-t model
#' the multivariate t with the replicate's one `nu`: a normal draw
#' divided by `sqrt(W / nu)`, `W ~ chi^2(nu)`, with one `W` per row
#' shared by the responses, which is what makes them jointly t rather
#' than K independent t's.
#'
#' `dps` is paired with `rspecs` by NAME when it has names, so the whole
#' list `eval_dpars()` returns serves any subset of the responses; an
#' unnamed `dps` is paired by position. Pairing by position alone gave a
#' subset `rspecs[c("y3", "y1")]` with the full list a draw labeled y3
#' that read y1's mean and sigma.
#'
#' @noRd
predict_rescor_draw <- function(fs, rspecs, dps, ok) {
  if (!is.null(names(dps))) {
    miss <- setdiff(names(rspecs), names(dps))
    if (length(miss)) {
      frm_stop("predict_rescor_draw(): `dps` has no entry for response '",
               miss[1L], "'. Pass the dpar values of every response in ",
               "`rspecs`, named by response", call. = FALSE)
    }
    dps <- dps[names(rspecs)]
  } else if (length(dps) != length(rspecs)) {
    frm_stop("predict_rescor_draw(): `dps` has ", length(dps), " entries ",
             "for ", length(rspecs), " responses. Name them by response, ",
             "or give one per response in the order of `rspecs`",
             call. = FALSE)
  }
  if (!any(ok)) return(NULL)
  R <- rescor_matrix(fs)[names(rspecs), names(rspecs), drop = FALSE]
  L <- chol(R)
  m <- sum(ok)
  Z <- matrix(stats::rnorm(m * length(rspecs)), m) %*% L
  nu <- rescor_shared_nu(fs)
  if (!is.null(nu)) Z <- Z / sqrt(stats::rchisq(m, nu) / nu)
  out <- list()
  for (k in seq_along(rspecs)) {
    dp <- dps[[k]]
    mu <- rep_len(as.numeric(dp[["mu"]]), length(ok))[ok]
    sg <- rep_len(as.numeric(dp[["sigma"]]), length(ok))[ok]
    out[[names(rspecs)[k]]] <- mu + sg * Z[, k]
  }
  out
}

#' The joint log-density of each row of a `rescor` fit, at the dpar
#' values `dpv` (every response's, as `eval_dpars()` returns them).
#'
#' Plain numeric R, written apart from the taped objective on purpose:
#' the Student-t branch is brms's `multi_student_t_lpdf()` term by term,
#' `lgamma((nu + K) / 2) - lgamma(nu / 2) - K / 2 log(nu pi)
#' - 1/2 log|Sigma| - (nu + K) / 2 log(1 + q / nu)`, with
#' `Sigma = D C D` entering as `log|C| + 2 sum(log sigma)`.
#'
#' The first two terms are NOT written as they read, and neither is
#' `log(nu pi)`. `lgamma_shift_diff()` and `log(nu) + log(pi)` are the
#' same numbers at a moderate `nu` and the only ones that survive a
#' large one, which a sampler reaches whenever `nu` is weakly
#' identified: against `mvt_std_loglik()`, the objective's own form, the
#' as-written version drifts by 5.6e-6 at `nu = 1e10`, by 227 log units
#' at `nu = 1e20`, and returns `NaN` from `Inf - Inf` at `nu = 1e306`,
#' where the density is still the gaussian limit and still wanted
#' (`dev/arcovsample-log/tnu.txt`). See `lgamma_shift_diff()` for why
#' the blend and the branch are both unavailable here.
#'
#' `resp` names a subset of the responses: brms's `log_lik(resp =)` of a
#' `rescor` model is the joint density of THOSE responses, the same law
#' restricted to the rows and columns of the correlation matrix they
#' name (a margin of a multivariate normal or t is one of the same
#' kind, with the same `nu`).
#'
#' @noRd
rescor_row_loglik <- function(fit, dpv, resp = NULL) {
  frame <- fit$frame
  all_rs <- names(fit$spec$responses)
  rs <- resp %||% all_rs
  K <- length(rs)
  n <- frame[["n_obs"]]
  Z <- matrix(vapply(rs, function(r) {
    as.numeric((frame[["y"]][[r]] - dpv[[r]]$mu) / dpv[[r]]$sigma)
  }, numeric(n)), n, K)
  lsig <- rowSums(matrix(vapply(rs, function(r) {
    rep(log(as.numeric(dpv[[r]]$sigma)), length.out = n)
  }, numeric(n)), n, K))
  C <- us_chol_cor(fit$estimates[["thetar"]], length(all_rs))
  # in the order `resp` gives, which is the order of the columns of Z
  if (!identical(rs, all_rs)) {
    k <- match(rs, all_rs)
    C <- C[k, k, drop = FALSE]
  }
  nu <- rescor_shared_nu(fit)
  if (is.null(nu)) {
    return(as.numeric(RTMB::dmvnorm(Z, 0, C, log = TRUE)) - lsig)
  }
  q <- rowSums((Z %*% solve(C)) * Z)
  ldet <- as.numeric(determinant(C, logarithm = TRUE)$modulus)
  lgamma_shift_diff(nu / 2, K / 2) - K / 2 * (log(nu) + log(pi)) -
    ldet / 2 - (nu + K) / 2 * log1p(q / nu) - lsig
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
#' brms's `sample_new_levels = "old_levels"` gives each unseen level the
#' effect of ONE seen level of its grouping factor, chosen at random once
#' per call and then read from every draw (`brms:::get_new_rdraws()`),
#' every term of that factor reading the same seen level. The
#' same holds here: the seen level is chosen here, before any replicate,
#' and each replicate reads that level's effect from its own draw of the
#' group effects, so the interval carries that level's conditional
#' uncertainty. The choice is recorded as the `old_pick` attribute.
#'
#' @noRd
predict_new_level_spec <- function(object, rspec, newdata, re_formula,
                                   allow_new_levels,
                                   sample_new_levels = "gaussian",
                                   pick = list()) {
  # brms draws a new latent residual for new rows, whatever re_formula
  # and allow_new_levels say: it is not a group-level term
  lat <- latres_new_spec(object, rspec, newdata)
  out <- predict_new_group_spec(object, rspec, newdata, re_formula,
                                allow_new_levels, sample_new_levels, pick)
  if (is.null(lat)) return(out)
  if (is.null(out)) out <- list()
  attr(out, "latres") <- lat
  out
}

#' `predict_new_level_spec()` for the group-level terms.
#'
#' @noRd
predict_new_group_spec <- function(object, rspec, newdata, re_formula,
                                   allow_new_levels, sample_new_levels,
                                   pick) {
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
      error = keep_missing(NULL))
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
  if (!length(out)) return(NULL)
  if (identical(sample_new_levels, "old_levels")) {
    # `pick` carries the choices another response of the same call made,
    # so a grouping factor two responses share is chosen once
    for (e in out) pick <- old_level_pick_add(pick, object, e$ed)
    attr(out, "old_pick") <- pick
  }
  out
}

#' brms's `"old_levels"` choice for the unseen levels of one design: one
#' seen level per grouping factor and unseen label, drawn once per call
#' as brms's `get_new_rdraws()` draws it, `sample(levels, 1)` in the
#' order the unseen levels first appear. brms runs that once per GROUP,
#' so every term of one grouping factor reads the same seen group; a key
#' already in `pick` keeps its level, which is what makes the second
#' block of `g` (a `sigma ~ (1 | g)` beside `mu`'s, or the slope half of
#' `(1 + x || g)`) read the level the first one chose. For
#' `gr(g, by = f)` the candidates are the levels of the row's own
#' by-level, as in brms, and a by-level's block that does not read the
#' row draws nothing. The value is the chosen level's label, which
#' `new_level_pick_apply()` and `predict_new_level_draw()` place in each
#' block. The keys depend on the design alone, so the estimates will
#' do. Blocks whose levels are their structure (`gr_cov`, `gr_prec`,
#' `car`, `spde`) have no unseen-level law, as in `lp_extra_var()`, and
#' are left out.
#'
#' @noRd
old_level_pick_add <- function(pick, object, ed) {
  for (rp in ed[["re_parts"]]) {
    bk <- rp$bk
    if (cs_levels_coupled(bk)) next
    rows <- which(rp$is_new %||% is.na(rp$j))
    if (!length(rows)) next
    keys <- new_level_key(bk, rp$new_key)
    lv <- as.character(bk[["levels"]])
    for (key in setdiff(unique(keys[rows]), names(pick))) {
      pick[[key]] <- lv[sample.int(length(lv), 1L)]
    }
  }
  pick
}

#' `fitted(sample_new_levels = "old_levels")`: the fit with the choice of
#' a seen level for each unseen one attached, where the design reads it
#' (`new_level_pick_apply()`). Every linear predictor of the responses
#' in `resp` is walked, nonlinear parameters and mixture weights
#' included: the design substitution is exact whatever transform follows
#' it, so none needs the exclusions `predict()`'s link-scale offset has.
#' A multivariate call passes all its responses at once, so a grouping
#' factor they share is chosen once, as brms chooses it.
#'
#' @noRd
fitted_old_levels <- function(object, resp, newdata, re_formula) {
  use_re <- re_form_keeps(re_formula)
  if (is.null(newdata) || !use_re) return(object)
  pick <- list()
  for (lp in object$frame[["linpreds"]] %||% list()) {
    if (!isTRUE(lp[["resp"]] %in% resp) || is.null(lp[["Z"]])) next
    ed <- tryCatch(suppressWarnings(
      lp_eta_design(object, lp, newdata, use_re, TRUE)),
      error = keep_missing(NULL))
    if (!is.null(ed)) pick <- old_level_pick_add(pick, object, ed)
  }
  if (length(pick)) object[["new_level_pick"]] <- pick
  object
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
  pick <- attr(spec, "old_pick")
  cvec <- if (!is.null(pick)) coef_b(fs)
  shared <- list()
  out <- list()
  for (dnm in names(spec)) {
    e <- spec[[dnm]]
    bl <- extra_var_blocks(lp_extra_var(fs, e$ed, TRUE)$new_levels, e$n)
    o <- numeric(e$n)
    for (key in names(bl)) {
      B <- bl[[key]]
      if (is.null(shared[[key]])) {
        shared[[key]] <- if (is.null(pick)) mvn_draw_cov(B$S) else {
          # the chosen seen level's effects in this replicate's draw;
          # a block's coefficients are laid out level by level. The
          # choice is keyed by grouping factor and label, the label
          # being what follows the block key in `key`
          d <- B$bk[["dim"]]
          lab <- pick[[new_level_key(B$bk, sub("^[^\r]*\r", "", key))]]
          k <- match(lab, as.character(B$bk[["levels"]]))
          # a by-level's block that does not read the row has no
          # choice, and its design rows are zero there
          if (is.na(k)) numeric(d) else {
            cvec[B$bk[["c_idx"]]][(k - 1L) * d + seq_len(d)]
          }
        }
      }
      o <- o + as.numeric(B$M %*% shared[[key]])
    }
    out[[dnm]] <- o
  }
  lat <- attr(spec, "latres")
  if (!is.null(lat)) {
    out[["mu"]] <- (out[["mu"]] %||% 0) + latres_new_draw(fs, lat)
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
  dpars_fn <- function(fit) {
    dpv <- list()
    for (dnm in names(rspec$dpars)) {
      off <- new_level_off[[dnm]]
      lp <- fit$frame[["linpreds"]][[linpred_key(resp, dnm)]]
      # the NAMES are kept: they are the data's row names, which every
      # other method carries into its output and which as.vector() drops
      if (is.null(off) && is.null(lp)) {
        dpv[[dnm]] <- drop(frm_linpred(fit, newdata = newdata, dpar = dnm,
                                       resp = resp, re_formula = re_formula,
                                       type = "response",
                                       allow_new_levels = allow_new_levels))
        next
      }
      # the link inverse of the predictor is the scale the simulator
      # reads (dpars_natural()). type = "response" is the REPORTING
      # scale, which for a mixture's theta is already the softmax
      # probability, and the simulator applied the softmax again: on
      # newdata a binomial mixture drew 443 where fitted() says 484. An
      # unseen level's drawn effect is an offset on the same scale,
      # where the model is additive
      off <- off %||% 0
      eta <- drop(frm_linpred(fit, newdata = newdata, dpar = dnm,
                              resp = resp, re_formula = re_formula,
                              type = "link",
                              allow_new_levels = allow_new_levels))
      dpv[[dnm]] <- lp[["link"]]$linkinv(eta + off)
    }
    dpv
  }
  # a cov = FALSE row whose response newdata leaves NA is filled with a
  # draw of this replicate before the rows after it read its residual
  dpv <- arma_cond_fill_dpars(fit, rspec, newdata, dpars_fn)
  # with_cs_offsets() takes the list of EVERY response and reads the
  # cs() values of the training rows; handed this one response's list
  # it wrote the offsets one level down, where no simulator reads them,
  # and predict() drew every cs() model as if the term were absent
  cs_offsets_add(fit, resp, newdata, with_shared_nu(fit, resp, dpv))
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
predict_par_drawer <- function(object, propagate_error, ndraws,
                               new_levels = FALSE) {
  if (!propagate_error) return(function(s) object)
  ds <- fit_draw_space(object)
  V <- ds$V
  # a draw cannot tell which predictions move along a lost direction, so
  # only a fit whose lost parameters are all variance components at
  # their boundary (held at zero, as the boundary message said) draws
  # from the covariance that holds them; any other loss keeps the
  # plug-in fallback and its warning
  if (!is.null(ds$jc) && !isTRUE(ds$boundary_only) &&
        !isTRUE(ds$mo_const_only)) {
    V <- NULL
  }
  # a "limit" that is not the family's own limiting case is a mean or
  # covariance parameter left undetermined by a fit with no maximum: its
  # raw variance (2.4e137 on an all-zero poisson) is no law to draw from
  if (any(ds$lost %in% "limit") && !draw_limit_quiet(object, ds$lost)) {
    V <- NULL
  }
  if (is.null(ds$root)) V <- draw_limit_restore(V, ds, object)
  # a covariance held as a factor (se_reparam_cov()) is drawn through
  # it: chol() of the whole matrix refuses one that ill-conditioned, or
  # returns a factor with none of its digits
  L <- if (!is.null(ds$root) && !is.null(V)) {
    t(ds$root)
  } else if (is.null(V) || !all(is.finite(V))) {
    NULL
  } else {
    tryCatch(chol(V + diag(0, nrow(V))), error = function(e) {
      if (is.null(ds$jc)) return(NULL)
      # held directions make the covariance singular: its symmetric
      # square root still draws from it
      e <- eigen(V, symmetric = TRUE)
      sqrt(pmax(e$values, 0)) * t(e$vectors)
    })
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
  Z <- matrix(stats::rnorm(p * ndraws), p, ndraws)
  law <- draw_edge_law(object, map, V, colSums(L^2), !is.null(ds$root),
                       rownames(V) %||% map$names)
  D <- v0 + if (is.null(law)) crossprod(L, Z) else draw_edge_apply(law, L, Z)
  held <- latres_held_pos(object, map)
  if (length(held)) D[held, ] <- v0[held]
  D <- draw_reweight(object, list(D = D, Z = Z, L = L, law = law, v0 = v0,
                                  held = held,
                                  labels = rownames(V) %||% map$names,
                                  comps = map$comp,
                                  new_levels = new_levels))
  out <- function(s) fit_set_outer(object, D[, s], map)
  # the group-effect draw conditions on this replicate's parameters, so
  # it needs how far they moved (predict_b_drawer())
  attr(out, "delta") <- function(s) D[, s] - v0
  out
}

#' The draw law, side by side, where the likelihood does not support the
#' normal one, or `NULL` where it does.
#'
#' The normal law N(theta_hat, V) reads the likelihood's curvature at the
#' estimates. Where an estimate sits at or near the end of its range, the
#' likelihood is flat on one side and falls far faster than that
#' curvature says on the other, and V's standard error measures only how
#' little it curves at the estimate. A zero-inflation intercept at -22.2
#' with no excess zeros in the data has a standard error of 4927: half
#' its normal draws put the zero inflation near 1, and predict() averaged
#' 3.04 where fitted() gave 5.74 (backlog 3269). One just inside its
#' range (-7.0, standard error 18.5) did the same (2.03 against 3.08).
#' A threshold gap of an ordinal mixture at exp(-13.1), standard error
#' 3840, made the thresholds after it infinite in half the draws (backlog
#' 3277). A group log sd at -3.0 with standard error 6.8 gave a new
#' level's prediction a standard deviation of 2.2e7.
#'
#' So each side of each parameter the objective reads is asked of the
#' likelihood, with the other parameters at the estimates
#' (draw_edge_side()):
#'
#' - a side where the log-likelihood changes by at most `draw_edge_flat`
#'   at `draw_edge_k` standard deviations is flat: the draws do not move
#'   the parameter that way (scale 0);
#' - a side where it falls by `draw_edge_k^2 / 2` (what the normal law
#'   gives at `draw_edge_k` conditional standard deviations) within
#'   `draw_edge_c` of that distance is steep: its draws are scaled by
#'   how much sooner the likelihood gets there, and a draw `z` sds out
#'   goes, in its own sds, no further than the one-parameter slice goes,
#'   in conditional sds, before the log-likelihood falls by `z^2 / 2`
#'   (draw_edge_map(); the slice's own point where the two sds agree);
#' - any other side keeps the normal law, so a fit whose likelihood the
#'   normal law describes keeps its draws bit for bit.
#'
#' A draw is never moved further than the normal law moves it: a side is
#' only ever narrowed. The conditional standard deviation `csd`
#' (draw_edge_csd(), the curve a one-parameter move sees) is the
#' reference, not the marginal one: an intercept beside an
#' uncentered slope falls about 100 times faster along its own axis than
#' its marginal standard deviation says, and is interior.
#'
#' A parameter with a flat side is drawn on its own, with its other side
#' scaled, and the rest are drawn from their normal law conditional on
#' it at its estimate (the Schur complement of V), then asked again under
#' that law. A slope on the zero inflation beside such an intercept
#' (`zi ~ x`) has a standard error of 4710 that its own two steep sides,
#' not the correlation (0.03), narrow.
#'
#' The verdict is kept in the fit's cache, and said once per fit (class
#' `frmtmb_draw_held`). It costs one objective evaluation per side of
#' each parameter where nothing changes, one more on a side that may be
#' flat, and up to about 50 on a steep side; on a fit with random
#' effects each is an inner solve.
#'
#' @noRd
draw_edge_law <- function(object, map, V, sd2, reparam, labels) {
  cache <- object$cache
  law <- if (is.environment(cache)) cache[["draw_edge_law"]]
  if (is.null(law)) {
    law <- draw_edge_build(object, map, V, sd2, reparam)
    if (is.environment(cache)) cache[["draw_edge_law"]] <- law
  }
  if (isFALSE(law)) return(NULL)
  if (!is.environment(cache) || !isTRUE(cache[["draw_edge_said"]])) {
    if (is.environment(cache)) cache[["draw_edge_said"]] <- TRUE
    frm_message(draw_edge_text(law, labels), class = "frmtmb_draw_held")
  }
  law
}

#' @noRd
draw_edge_text <- function(law, labels) {
  q <- function(i) paste0("`", labels[i], "`", collapse = ", ")
  paste0(
    "predict(): ",
    if (length(law$edge)) {
      paste0("the log-likelihood is flat toward one end of the range of ",
             q(law$edge), ", so the draws do not move ",
             if (length(law$edge) > 1L) "them" else "it",
             " that way and draw the other parameters given ",
             if (length(law$edge) > 1L) "their estimates" else
               "its estimate", ". ")
    },
    if (length(law$narrowed)) {
      paste0("The draws of ", q(law$narrowed), " are narrowed on a side ",
             "where the log-likelihood falls faster than the normal law ",
             "of vcov() says, so that they do not reach values the data ",
             "rule out. ")
    },
    "Said once per fit")
}

#' @noRd
draw_edge_build <- function(object, map, V, sd2, reparam = FALSE) {
  obj <- object$obj
  p <- object$opt$par
  np <- length(p)
  P <- length(map$names)
  off <- P - np
  if (is.null(obj) || !is.function(obj$fn) || off < 0L ||
        length(sd2) != P || NROW(V) != P) {
    return(FALSE)
  }
  sdv <- sqrt(pmax(sd2, 0))
  # before any evaluation moves the point the cached Hessian is checked at
  csd <- draw_edge_csd(object, V, sdv, off, reparam)
  if (is.null(csd)) return(FALSE)
  saved <- obj_state_save(obj)
  on.exit(obj_state_restore(obj, saved), add = TRUE)
  f0 <- tryCatch(obj$fn(p), error = function(e) NA_real_)
  if (!is.finite(f0)) return(FALSE)
  pos <- off + seq_len(np)
  lo <- hi <- as.list(rep(1, P))
  sd_law <- sdv
  sides <- draw_edge_sides(obj, p, off, pos, sdv, csd, f0)
  lo[pos] <- sides$lo
  hi[pos] <- sides$hi
  flat <- function(v) length(v) == 1L && v == 0
  edge <- pos[vapply(lo[pos], flat, TRUE) | vapply(hi[pos], flat, TRUE)]
  Lc <- NULL
  if (length(edge)) {
    Vc <- draw_edge_condition(V, edge)
    Lc <- draw_root(Vc)
    rest <- setdiff(pos, edge)
    sd_law[rest] <- sqrt(pmax(diag(Vc)[rest], 0))
    sides <- draw_edge_sides(obj, p, off, rest, sd_law, csd, f0)
    lo[rest] <- sides$lo
    hi[rest] <- sides$hi
  }
  steep <- function(v) length(v) > 1L
  narrowed <- setdiff(pos[vapply(lo[pos], steep, TRUE) |
                            vapply(hi[pos], steep, TRUE)], edge)
  if (!length(edge) && !length(narrowed)) return(FALSE)
  list(lo = lo, hi = hi, s = sd_law, Lc = Lc, edge = edge,
       narrowed = narrowed)
}

#' The law of each side, `list(lo, hi)` of the draw positions `js`, each
#' read at its own standard deviation `sdv[j]` (draw_edge_side()).
#'
#' @noRd
draw_edge_sides <- function(obj, p, off, js, sdv, csd, f0) {
  lo <- hi <- as.list(rep(1, length(js)))
  for (i in seq_along(js)) {
    j <- js[i]
    s <- sdv[j]
    if (!is.finite(s) || s <= 0) next
    cs <- min(csd[j], s)
    lo[[i]] <- draw_edge_side(obj, p, j - off, -1, s, cs, f0)
    hi[[i]] <- draw_edge_side(obj, p, j - off, 1, s, cs, f0)
  }
  list(lo = lo, hi = hi)
}

#' One side's law: 0 where it is flat; 1 where the normal law stands;
#' on a steep side, the distances (in conditional standard deviations
#' `cs`) at which the log-likelihood falls by `z^2 / 2`, z = 1 to 4, each
#' at most `z`. A move the objective cannot evaluate counts as a steep
#' one: the data rule it out.
#'
#' @noRd
draw_edge_side <- function(obj, p, k, dir, s, cs, f0) {
  at <- function(t) {
    q <- p
    q[k] <- p[k] + dir * t
    f1 <- tryCatch(suppressWarnings(obj$fn(q)), error = function(e) NA_real_)
    if (is.finite(f1)) f1 - f0 else Inf
  }
  K <- draw_edge_k
  tc <- min(draw_edge_c * K * cs, K * s)
  dc <- at(tc)
  if (dc >= K^2 / 2) {
    tk <- vapply(1:4, function(z) draw_edge_reach(at, z^2 / 2, z * cs), 0)
    return(cummax(tk) / cs)
  }
  if (abs(dc) <= draw_edge_flat &&
        (tc >= K * s || abs(at(K * s)) <= draw_edge_flat)) {
    return(0)
  }
  1
}

#' The distance at which `at()` first reaches `target`, at most `cap`: a
#' geometric search down from `cap`, then six halvings of the bracket in
#' log scale, which leaves it to a few percent.
#'
#' @noRd
draw_edge_reach <- function(at, target, cap) {
  if (at(cap) < target) return(cap)
  hi <- cap
  lo <- cap
  repeat {
    lo <- lo / 4
    if (lo < cap * 1e-8) return(lo)
    if (at(lo) < target) break
    hi <- lo
  }
  for (i in seq_len(6L)) {
    mid <- sqrt(lo * hi)
    if (at(mid) < target) lo <- mid else hi <- mid
  }
  sqrt(lo * hi)
}

#' Deviations `x` on one side of a position whose draw law has standard
#' deviation `s`, mapped through that side's law `v` (draw_edge_side()).
#' `v` holds, in conditional standard deviations `cs`, the distances
#' along the one-parameter slice at which the log-likelihood falls by
#' `z^2 / 2`, z = 1 to 4 (linear between, and past 4 at the last knots'
#' slope). A draw `z = |x| / s` out goes to `s * g`, `g` the least of
#' `z` (the normal law), `z * v[2] / 2` (the side's one scale) and `v`
#' at `z`: the draw's distance in its own sds is at most the slice's in
#' conditional sds. That is the point where the slice falls by `z^2 / 2`
#' only where `s = cs`; where the parameter is correlated with others
#' (`s > cs`) the draw goes `s / cs` times as far along this axis, as the
#' normal law's own draws do, and the slice is not where it lands.
#'
#' Either alone failed a case. The one scale lets the normal tail past
#' z = 2 run through a wall: a log sd with the residual sd at 0 under
#' se(sigma = TRUE) drew a row of standard error 0.05 with a standard
#' deviation of 0.61, where brms 2.23.0 has 0.108. The knots alone
#' follow an exponential wall too closely inside it: on 3269's zero
#' inflation they gave half brms's share of zeros (0.0055 against
#' 0.0097), and on `zi ~ x` a mean 8 Monte Carlo standard errors from
#' fitted() against 3.4 (dev/ppdraw-findings.md, "Fix round").
#'
#' A side only just steep enough is blended back toward the normal law
#' (draw_edge_blend()), so that the map does not jump at `draw_edge_c`.
#'
#' @noRd
draw_edge_map <- function(x, s, v) {
  if (length(v) == 1L) return(x * v)
  z <- abs(x) / s
  g <- stats::approx(0:4, c(0, v), xout = pmin(z, 4))$y
  far <- z > 4
  g[far] <- v[4] + (z[far] - 4) * max(0, v[4] - v[3])
  h <- sign(x) * s * pmin(z, z * v[2] / 2, g)
  a <- draw_edge_blend(v)
  if (a < 1) h <- (1 - a) * x + a * h
  h
}

#' The weight of a steep side's narrowing against the normal law: 1 where
#' the side's ratio r, `v[2] / 2` (the distance at which the
#' log-likelihood falls by 2, over the normal law's), is at most
#' `draw_edge_r0`, falling linearly to 0 at `draw_edge_c`, where a side
#' stops being steep. Without it the 98% width of a new level's
#' prediction stepped from 5.36 to 5.76 between neighboring data sets
#' (backlog 4932); with it the steep sides with r in (0.4, 0.5] covered
#' 0.902 and 0.981 at 90% and 98%, where the full narrowing covered 0.886
#' and 0.974 (dev/predraw2-findings.md).
#'
#' @noRd
draw_edge_blend <- function(v) {
  r <- v[2] / 2
  min(1, max(0, (draw_edge_c - r) / (draw_edge_c - draw_edge_r0)))
}

#' The standard deviation of each draw position along its own axis, the
#' others held: `1 / sqrt(H[j, j])` from the objective's Hessian in the
#' outer parameters where the fit has it (`fit$cache$hessian_fixed`,
#' taken at the estimates), which is exactly the curve a one-parameter
#' move of the objective sees, under REML too; else `1 / sqrt(diag(V^-1))`.
#' `NULL` where neither can be trusted: a covariance held as a factor
#' (se_reparam_cov()) is one solve() takes every digit from. On a raw
#' degree-8 polynomial it narrowed draws 20-fold that the normal law
#' describes.
#'
#' @noRd
draw_edge_csd <- function(object, V, sdv, off, reparam) {
  np <- NROW(V) - off
  hc <- if (is.environment(object$cache)) object$cache$hessian_fixed
  H <- hc$H
  out <- sdv
  if (is.null(hc$Hq) && is.matrix(H) && NROW(H) == np &&
        identical(hc$at, sdr_outer_point(object$obj))) {
    # a curvature that rounds below zero (a correlation held at its limit
    # reads +-1e-9) has no conditional sd: Inf here, and the marginal one
    # stands, with no "NaNs produced" from sqrt()
    cs <- 1 / sqrt(pmax(diag(H), 0))
    ok <- is.finite(cs) & cs > 0
    j <- off + which(ok)
    out[j] <- pmin(cs[ok], sdv[j])
    return(out)
  }
  if (reparam) return(NULL)
  s <- which(is.finite(sdv) & sdv > 0)
  if (!length(s)) return(out)
  Q <- tryCatch(solve(V[s, s, drop = FALSE]), error = function(e) NULL)
  if (is.null(Q)) return(out)
  cs <- 1 / sqrt(pmax(diag(Q), 0))
  ok <- is.finite(cs) & cs > 0
  out[s[ok]] <- pmin(cs[ok], sdv[s[ok]])
  out
}

#' The draw covariance with the positions `edge` drawn on their own (their
#' variances, no covariance) and the rest conditional on them at their
#' estimates.
#'
#' @noRd
draw_edge_condition <- function(V, edge) {
  rest <- setdiff(seq_len(NROW(V)), edge)
  Vee <- V[edge, edge, drop = FALSE]
  W <- tryCatch(solve(Vee, V[edge, rest, drop = FALSE]),
                error = function(e) {
                  ev <- eigen(Vee, symmetric = TRUE)
                  inv <- ifelse(ev$values > max(ev$values) * 1e-12,
                                1 / ev$values, 0)
                  ev$vectors %*% (inv * t(ev$vectors)) %*%
                    V[edge, rest, drop = FALSE]
                })
  Vc <- matrix(0, NROW(V), NROW(V), dimnames = dimnames(V))
  Vc[rest, rest] <- V[rest, rest, drop = FALSE] -
    crossprod(V[edge, rest, drop = FALSE], W)
  Vc[rest, rest] <- (Vc[rest, rest] + t(Vc[rest, rest])) / 2
  diag(Vc)[edge] <- diag(V)[edge]
  Vc
}

#' A square root `L` of a covariance, `crossprod(L) == V`, through its
#' symmetric root where held positions make it singular.
#'
#' @noRd
draw_root <- function(V) {
  tryCatch(chol(V), error = function(e) {
    ev <- eigen(V, symmetric = TRUE)
    sqrt(pmax(ev$values, 0)) * t(ev$vectors)
  })
}

#' The draws' deviations from the estimates under a side-by-side law:
#' the normal draw through `law$Lc` (the conditional law) or `L`, each
#' coordinate then mapped through the law of the side it fell on.
#'
#' @noRd
draw_edge_apply <- function(law, L, Z) {
  delta <- crossprod(law$Lc %||% L, Z)
  one <- function(v) length(v) == 1L && v == 1
  for (j in which(!vapply(law$lo, one, TRUE) | !vapply(law$hi, one, TRUE))) {
    x <- delta[j, ]
    neg <- x < 0
    x[neg] <- draw_edge_map(x[neg], law$s[j], law$lo[[j]])
    x[!neg] <- draw_edge_map(x[!neg], law$s[j], law$hi[[j]])
    delta[j, ] <- x
  }
  delta
}

#' The draws, reweighted where the log-likelihood says the draw law is
#' wrong inside the range, and as they came everywhere else.
#'
#' The side law (draw_edge_law()) reads one parameter at a time, with
#' the others at their estimates. A weakly identified mixture has no
#' parameter at an end, yet its log-likelihood falls far faster than the
#' normal law says along directions no single parameter shows: on
#' one-component data a two-component gaussian mixture's draws fell by a
#' median of 5.7 where the law says 3.2, and predict()'s mean left
#' fitted()'s by 13 Monte Carlo standard errors while the simulator was
#' right (backlog 4882).
#'
#' The check (draw_check()) is paired, draw by draw: at each test draw
#' the objective's drop from its minimum is compared with the drop the
#' draw law implies at that same point, half the draw's Mahalanobis
#' distance less the side law's log Jacobian (draw_logratio()). Where
#' the objective is the law's quadratic the two are equal at every draw,
#' so the comparison has no sampling noise of its own, and a chi-square
#' reference for the drops alone would add some. It fires where the
#' median excess drop over `draw_check_n` test draws (twice as many
#' near the thresholds) is above `draw_check_tol` and above
#' `draw_check_rel` times the median drop the law implies
#' (draw_check_verdict()). That median is a sum over every direction, so
#' a block of well-identified parameters dilutes one weak direction
#' until it no longer shows; where it does not fire, each parameter's
#' direction of the law is read on its own (draw_check_dirs()).
#' The test draws come from a private stream, so the
#' verdict belongs to the fit: it is cached, and where it does not fire
#' the draws are returned untouched.
#'
#' Where it fires, each draw gets the importance ratio exp(-objective) /
#' proposal density, Pareto smoothed (loo::psis()), and the draws are
#' resampled by those weights (draw_resample()). Above a Pareto k of
#' `draw_pareto_tol` no normal-based proposal reaches the likelihood, and
#' a warning says so and names frmtmb.sample::frm_sample().
#'
#' Where only directions of variance components fire the per-direction
#' reading (`warn_only`), the draws are returned as they came and a
#' warning says so (draw_tail_say()): on 12 random-slope fires the
#' reweighted draws made the in-sample predictive means worse against
#' the truth on 10, and narrowed a new level's 90% interval to 0.89 of a
#' sampled posterior's, while the normal law's 90% interval matched it;
#' the likelihood's flat side toward a variance of 0 leaves the
#' reweighting a target with no proper mass (dev/drawdir-findings.md,
#' section 12).
#'
#' @noRd
draw_reweight <- function(object, dw) {
  dw$off <- length(dw$v0) - length(object$opt$par)
  chk <- draw_check(object, dw)
  if (!isTRUE(chk$fire)) return(dw$D)
  if (isTRUE(chk$warn_only)) {
    draw_tail_say(chk, isTRUE(dw$new_levels))
    return(dw$D)
  }
  if (!requireNamespace("loo", quietly = TRUE)) {
    frm_warning(draw_check_said_what(chk), ", and the draws need the loo ",
                "package to be reweighted. They are NOT corrected; install ",
                "loo", class = "frmtmb_draw_pareto", call. = FALSE)
    return(dw$D)
  }
  r <- draw_logratio(object, dw, dw$D, dw$Z)
  if (is.null(r)) return(dw$D)
  rs <- draw_resample(r$lw, r$stratum)
  draw_reweight_say(object, chk, rs$k, ncol(dw$D), isTRUE(dw$new_levels))
  dw$D[, rs$idx, drop = FALSE]
}

#' The fit's verdict on its draw law, cached: `list(fire, excess, rel,
#' n, dir, dir_par, dir_n, warn_only)`, `dir` the per-direction reading
#' (draw_check_dirs(); NA where the summed check fires first) and
#' `warn_only` TRUE where every direction that fired is a variance
#' component's.
#' A median within a factor of 2 of the threshold is read again over
#' `draw_check_n` more test draws, deciding on all of them: on a random
#' slope fit whose excess over 1000 draws has median 0.11 and a long
#' right tail, 12% of random 20-draw medians were above 1
#' (dev/predraw2-findings.md).
#'
#' @noRd
draw_check <- function(object, dw) {
  cache <- object$cache
  chk <- if (is.environment(cache)) cache[["draw_check"]]
  if (is.null(chk)) {
    chk <- list(fire = FALSE, excess = NA_real_, rel = NA_real_, n = 0L)
    Z0 <- draw_check_z(length(dw$v0), 2L * draw_check_n)
    D0 <- dw$v0 + if (is.null(dw$law)) crossprod(dw$L, Z0) else
      draw_edge_apply(dw$law, dw$L, Z0)
    if (length(dw$held)) D0[dw$held, ] <- dw$v0[dw$held]
    first <- seq_len(draw_check_n)
    r <- draw_logratio(object, dw, D0[, first, drop = FALSE],
                       Z0[, first, drop = FALSE])
    if (!is.null(r)) {
      e <- r$e
      m <- r$m
      v <- draw_check_verdict(e, m)
      if (v$near) {
        r2 <- draw_logratio(object, dw, D0[, -first, drop = FALSE],
                            Z0[, -first, drop = FALSE])
        if (!is.null(r2)) {
          e <- c(e, r2$e)
          m <- c(m, r2$m)
          v <- draw_check_verdict(e, m)
        }
      }
      chk <- list(fire = v$fire, excess = v$excess, rel = v$rel,
                  n = length(e), dir = NA_real_, dir_par = NA_character_,
                  dir_n = 0L, warn_only = FALSE)
      if (!isTRUE(v$fire)) {
        dd <- draw_check_dirs(object, dw)
        if (!is.null(dd)) {
          chk$dir <- dd$stat
          chk$dir_par <- dd$par
          chk$dir_n <- dd$n
          chk$fire <- isTRUE(dd$stat > draw_dir_tol)
          chk$warn_only <- chk$fire && isTRUE(dd$vc_only)
        }
      }
    }
    if (is.environment(cache)) cache[["draw_check"]] <- chk
  }
  chk
}

#' The per-direction reading of the draw law: `list(stat, par, n,
#' vc_only)`, or `NULL` where the objective cannot be read. `vc_only` is
#' TRUE where every direction above `draw_dir_tol` belongs to a variance
#' component (`draw_vc_comps`).
#'
#' The summed check's median is over draws that move every parameter at
#' once, so a block of K well-identified coordinates adds about K to the
#' law's drop and nothing to the excess, and one weak direction stops
#' showing: on 18 weakly identified mixtures the check stopped firing at
#' K = 8 to 128 added coordinates (backlog 5584). Here each parameter j
#' of the objective is read along its own direction of the law: j moved
#' `draw_edge_k` of its standard deviations, the others by their
#' regression on it (column j of the law's factor), through the side
#' law, on each side. The law's drop there is `draw_edge_k^2 / 2`
#' whatever the other parameters do, so a separable block adds nothing.
#'
#' The side law already reads each parameter's own slice (the others at
#' the estimates). A group log sd near 0 fell 15.2 further than the law
#' says at 2 standard deviations along its direction, 10.8 of it along
#' its own slice, a side the side law leaves to the normal law. So the
#' slice's positive excess eB, read at `draw_edge_k` of j's CONDITIONAL
#' standard deviations, is taken off the direction's eA: `stat = eA -
#' max(eB, 0)`. The slice is read at a shorter distance than the
#' direction (the ratio of j's standard deviation to its conditional one,
#' `sratio`, is 1 or more), so this is not "the part no single parameter
#' shows": for a skew in j alone that grows as the cube of the distance,
#' `stat = eB (sratio^3 - 1)`, 0 where j is uncorrelated and most of the
#' error where it is strongly correlated. It is a heuristic that keeps
#' the GLMMs' group sds, which the side law judges, out of the reading
#' (dev/drawdir-findings.md, and the review of 2026-11-17).
#'
#' A direction whose excess eA is not above `draw_dir_tol` needs no
#' slice (stat is at most eA), so `stat` is eA there, an upper bound. A
#' direction whose end the objective cannot evaluate is read again at
#' `draw_dir_back` of the distance, its excess put on the full distance's
#' scale by `1 / t^2` (exact for a quadratic misfit), and left unread if
#' none is finite: a random-slope fit's log slope sd moved 76 log units
#' out failed in the inner solve, along a direction flatter than the
#' law. A slice the objective cannot evaluate leaves its direction
#' unread too. `stat` is never infinite.
#'
#' Cost: two objective evaluations per parameter, one more for each
#' direction above the threshold and each end read again; once per fit.
#'
#' @noRd
draw_check_dirs <- function(object, dw) {
  obj <- object$obj
  p <- object$opt$par
  np <- length(p)
  off <- dw$off
  if (is.null(obj) || !is.function(obj$fn) || off < 0L || np == 0L) {
    return(NULL)
  }
  pos <- off + seq_len(np)
  Lf <- dw$law$Lc %||% dw$L
  sv <- colSums(Lf^2)
  R <- setdiff(pos[sv[pos] > 0], c(dw$held, dw$law$edge))
  k <- length(R)
  if (!k) return(NULL)
  B <- Lf[, R, drop = FALSE]
  U <- B / rep(sqrt(sv[R]), each = nrow(B))
  ck <- draw_edge_k
  Z <- cbind(ck * U, -ck * U)
  D <- dw$v0 + if (is.null(dw$law)) crossprod(dw$L, Z) else
    draw_edge_apply(dw$law, dw$L, Z)
  if (length(dw$held)) D[dw$held, ] <- dw$v0[dw$held]
  r <- draw_logratio(object, dw, D, Z)
  if (is.null(r)) return(NULL)
  eA <- r$e
  n <- 2L * k
  # the share of draw_edge_k each direction was read at; the excess is
  # put on draw_edge_k's scale by 1 / t^2, as a quadratic misfit scales
  at <- rep(1, 2L * k)
  for (t in draw_dir_back) {
    bad <- which(!is.finite(eA))
    if (!length(bad)) break
    Zt <- t * Z[, bad, drop = FALSE]
    Dt <- dw$v0 + if (is.null(dw$law)) crossprod(dw$L, Zt) else
      draw_edge_apply(dw$law, dw$L, Zt)
    if (length(dw$held)) Dt[dw$held, ] <- dw$v0[dw$held]
    rt <- draw_logratio(object, dw, Dt, Zt)
    if (is.null(rt)) break
    n <- n + length(bad)
    eA[bad] <- rt$e / t^2
    at[bad] <- t
  }
  # a direction the objective cannot read even near the estimate is left
  # unread: an end the inner solve fails at is no evidence of a wall
  eA[!is.finite(eA)] <- -Inf
  stat <- eA
  cand <- which(eA > draw_dir_tol)
  if (length(cand)) {
    csd <- draw_dirs_csd(B)
    saved <- obj_state_save(obj)
    on.exit(obj_state_restore(obj, saved), add = TRUE)
    f0 <- tryCatch(obj$fn(p), error = function(e) NA_real_)
    if (!is.finite(f0)) return(NULL)
    for (i in cand) {
      j <- (i - 1L) %% k + 1L
      s <- at[i] * ck
      q <- p
      q[R[j] - off] <- p[R[j] - off] + (if (i > k) -1 else 1) * s * csd[j]
      obj_state_restore(obj, saved)
      f1 <- tryCatch(suppressWarnings(obj$fn(q)),
                     error = function(e) NA_real_)
      n <- n + 1L
      stat[i] <- if (!is.finite(csd[j]) || !is.finite(f1)) {
        -Inf
      } else {
        eA[i] - max(f1 - f0 - s^2 / 2, 0) / at[i]^2
      }
    }
  }
  i <- which.max(stat)
  if (!length(i) || !is.finite(stat[i])) return(NULL)
  lab <- dw$labels %||% character(0)
  par <- if (length(lab) >= R[(i - 1L) %% k + 1L]) {
    lab[R[(i - 1L) %% k + 1L]]
  } else {
    NA_character_
  }
  fired <- R[(which(stat > draw_dir_tol) - 1L) %% k + 1L]
  comp <- dw$comps %||% character(0)
  vc <- length(fired) > 0L && length(comp) >= max(fired) &&
    all(comp[fired] %in% draw_vc_comps)
  list(stat = stat[i], par = par, n = n, vc_only = vc)
}

#' Each column's conditional standard deviation given the others, under
#' the law whose factor holds those columns, `1 / sqrt(diag(solve(B'B)))`
#' through a QR of `B`, so the covariance itself is never inverted.
#'
#' @noRd
draw_dirs_csd <- function(B) {
  k <- ncol(B)
  qf <- qr(B)
  if (qf$rank < k) return(rep(NA_real_, k))
  Ri <- backsolve(qr.R(qf), diag(k))
  out <- numeric(k)
  out[qf$pivot] <- 1 / sqrt(rowSums(Ri^2))
  out
}

#' What the check found, for its message and warning.
#'
#' @noRd
draw_check_said_what <- function(chk) {
  if (isTRUE(chk$dir > draw_dir_tol)) {
    return(paste0(
      "predict(): along one direction of the normal law of vcov(), the one ",
      "that moves ", if (is.na(chk$dir_par)) "a parameter" else {
        paste0("`", chk$dir_par, "`")
      }, " and the parameters correlated with it, the log-likelihood ",
      "falls further than that law says at ", draw_edge_k, " standard ",
      "deviations (an excess of ", format(chk$dir, digits = 3), " after ",
      "the excess of that parameter moved alone, nearer the estimate, is ",
      "taken off)"))
  }
  paste0("predict(): the log-likelihood at draws from the normal law of ",
         "vcov() falls further than that law says (median excess ",
         format(chk$excess, digits = 3), " over ", chk$n, " test draws)")
}

#' The warning where only variance components' directions fired: the
#' draws are the normal law's, said on every call, as the Pareto warning
#' is.
#'
#' @noRd
draw_tail_say <- function(chk, new_levels = FALSE) {
  who <- if (is.na(chk$dir_par)) "a variance component" else {
    paste0("`", chk$dir_par, "`, a variance component,")
  }
  frm_warning(draw_check_said_what(chk), ". So the normal law's far ",
              "tail along ", who, " is too wide, and the draws' extreme ",
              "quantiles (a new level's 98% interval, for one) are too ",
              "wide with it. The draws are NOT reweighted: toward a ",
              "variance of 0 the likelihood is flat, and weights by it ",
              "made the predictions worse. ",
              "frmtmb.sample::frm_sample() draws from the posterior itself",
              if (new_levels) {
                paste0(", but its posterior_predict() cannot yet predict ",
                       "a grouping level the fit did not see")
              },
              class = "frmtmb_draw_wide_tail", call. = FALSE)
  invisible(NULL)
}

#' The check's verdict from the test draws' excesses `e` and Mahalanobis
#' distances `m`. Two conditions, both needed:
#'
#' - the median excess is above `draw_check_tol`: half the draws land
#'   where the likelihood is more than e times lower than the law says;
#' - the median excess over the median drop the law implies, `median(m)
#'   / 2`, is above `draw_check_rel`. A small skew in each of p
#'   coordinates adds up to an excess that grows like p, while that
#'   ratio stays near a quarter of each coordinate's standardized
#'   fourth-order term whatever p is. A weakly identified mixture has a
#'   few directions far from the law and a ratio of 0.3 or more.
#'
#' `near` asks for more test draws: either quantity within a factor of 2
#' of its threshold while the other is at least half of its own.
#'
#' @noRd
draw_check_verdict <- function(e, m) {
  ex <- stats::median(e)
  half <- stats::median(m) / 2
  rel <- if (is.finite(ex) && half > 0) ex / half else if (is.na(ex)) {
    NA_real_
  } else {
    Inf
  }
  a <- ex / draw_check_tol
  b <- rel / draw_check_rel
  ok <- !is.na(a) && !is.na(b)
  list(fire = ok && a > 1 && b > 1, excess = ex, rel = rel,
       near = ok && a >= 0.5 && b >= 0.5 && (a <= 2 || b <= 2))
}

#' Standard normal deviates for draw_check()'s test draws, from a fixed
#' seed, leaving the caller's stream (and its generator kind) as it was.
#'
#' @noRd
draw_check_z <- function(p, n) {
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = genv)
  on.exit(if (had) {
    assign(".Random.seed", old, envir = genv)
  } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
    rm(".Random.seed", envir = genv)
  }, add = TRUE)
  set.seed(draw_check_seed, kind = "Mersenne-Twister",
           normal.kind = "Inversion", sample.kind = "Rejection")
  matrix(stats::rnorm(p * n), p, n)
}

#' The log importance ratio of each draw (column of `D`, made from the
#' deviates `Z` by the law in `dw`): `lw` = -(objective drop) + m / 2 +
#' log Jacobian, and `e` = -lw anchored to 0 at the estimates, the
#' excess of the drop over what the law implies. `m` is the draw's
#' Mahalanobis distance under the normal law the deviates went through,
#' taken through the factor (a QR of its columns) so that no inverse of
#' the covariance is formed. A side the law narrows maps its deviate
#' through draw_edge_map(), whose slope is the Jacobian. A flat side
#' holds its parameter at the estimate: those draws form their own
#' stratum, as the likelihood along a flat side has no mass to set
#' against the other side's. `NULL` where the objective cannot be read.
#'
#' Under REML the drawn fixed effects are not arguments of the
#' objective; the ratio is that of the outer parameters' marginal, and
#' the fixed effects keep their conditional law given them.
#'
#' @noRd
draw_logratio <- function(object, dw, D, Z) {
  obj <- object$obj
  p <- object$opt$par
  np <- length(p)
  off <- dw$off
  if (is.null(obj) || !is.function(obj$fn) || off < 0L || np == 0L) {
    return(NULL)
  }
  pos <- off + seq_len(np)
  saved <- obj_state_save(obj)
  on.exit(obj_state_restore(obj, saved), add = TRUE)
  f0 <- tryCatch(obj$fn(p), error = function(e) NA_real_)
  if (!is.finite(f0)) return(NULL)
  drop <- vapply(seq_len(ncol(D)), function(s) {
    # every inner solve starts from the estimates' modes: a far draw left
    # as the start of the next one can send that one astray
    obj_state_restore(obj, saved)
    f1 <- tryCatch(suppressWarnings(obj$fn(D[pos, s])),
                   error = function(e) NA_real_)
    if (is.finite(f1)) f1 - f0 else Inf
  }, 0)
  law <- dw$law
  L <- law$Lc %||% dw$L
  X <- crossprod(L, Z)
  sv <- colSums(L^2)
  edge <- law$edge %||% integer(0)
  R <- setdiff(pos[sv[pos] > 0], c(dw$held, edge))
  m <- numeric(ncol(Z))
  if (length(R)) {
    qf <- qr(L[, R, drop = FALSE])
    if (qf$rank > 0L) {
      m <- colSums(qr.qty(qf, Z)[seq_len(qf$rank), , drop = FALSE]^2)
    }
  }
  # la, the log Jacobian itself, is the proposal's; lj, anchored to 0 at
  # the estimate on each side, is the check's, which a side law that
  # follows its slice exactly would otherwise read as an excess
  lj <- la <- numeric(ncol(Z))
  stratum <- character(ncol(Z))
  if (!is.null(law)) {
    one <- function(v) length(v) == 1L && v == 1
    mapped <- which(!vapply(law$lo, one, TRUE) | !vapply(law$hi, one, TRUE))
    for (j in intersect(mapped, c(R, edge))) {
      x <- X[j, ]
      neg <- x < 0
      sl <- s0 <- numeric(length(x))
      sl[neg] <- draw_edge_slope(abs(x[neg]) / law$s[j], law$lo[[j]])
      sl[!neg] <- draw_edge_slope(abs(x[!neg]) / law$s[j], law$hi[[j]])
      s0[neg] <- draw_edge_slope(0, law$lo[[j]])
      s0[!neg] <- draw_edge_slope(0, law$hi[[j]])
      mv <- sl > 0
      if (j %in% edge) {
        stratum <- paste0(stratum, ifelse(mv, "s", "f"))
        m[mv] <- m[mv] + x[mv]^2 / sv[j]
      }
      la[mv] <- la[mv] + log(sl[mv])
      lj[mv] <- lj[mv] + log(sl[mv]) - log(s0[mv])
    }
  }
  list(e = drop - m / 2 - lj, lw = -drop + m / 2 + la, m = m,
       stratum = stratum)
}

#' The slope of draw_edge_map()'s map at `z` (a deviate in its own
#' standard deviations) on a side with law `v`: of whichever of its
#' three pieces is lowest there, and the first knot's at `z = 0`, blended
#' with the normal law's 1 as the map is.
#'
#' @noRd
draw_edge_slope <- function(z, v) {
  if (length(v) == 1L) return(rep(v, length(z)))
  if (!length(z)) return(numeric(0))
  vv <- c(0, v)
  tail <- max(0, v[4] - v[3])
  k <- pmin(floor(z), 3)
  sg <- ifelse(z >= 4, tail, vv[k + 2] - vv[k + 1])
  g <- stats::approx(0:4, vv, xout = pmin(z, 4))$y
  g[z > 4] <- v[4] + (z[z > 4] - 4) * tail
  # at z = 0 every piece is 0; the lowest slope is the one that leads
  a <- cbind(z, z * v[2] / 2, g)
  s <- cbind(1, v[2] / 2, sg)
  w <- ifelse(z > 0, max.col(-a, ties.method = "first"),
              max.col(-s, ties.method = "first"))
  sl <- s[cbind(seq_along(z), w)]
  b <- draw_edge_blend(v)
  if (b < 1) sl <- (1 - b) + b * sl
  sl
}

#' Draw indices resampled by Pareto-smoothed importance weights, each
#' stratum on its own and keeping its count, by systematic resampling (one
#' uniform from the caller's stream per stratum). `k` is the largest
#' Pareto k of the strata, `Inf` where one has too few draws to fit it.
#'
#' @noRd
draw_resample <- function(lw, stratum) {
  idx <- seq_along(lw)
  k <- -Inf
  for (st in unique(stratum)) {
    j <- which(stratum == st)
    if (!any(is.finite(lw[j]))) {
      k <- Inf
      next
    }
    ps <- suppressWarnings(loo::psis(lw[j], r_eff = 1))
    kj <- ps$diagnostics$pareto_k
    # loo returns Inf where every weight of its tail is the same; such
    # weights are bounded, which is k at most 0, not too few draws
    tl <- utils::tail(sort(lw[j]), ceiling(min(0.2 * length(j),
                                               3 * sqrt(length(j)))))
    if (!is.finite(kj) && length(j) >= 25L && all(is.finite(tl)) &&
          diff(range(tl)) <= 1e-12 * max(1, abs(tl[1]))) {
      kj <- 0
    }
    k <- max(k, kj)
    lwj <- as.numeric(ps$log_weights)
    w <- exp(lwj - max(lwj))
    cw <- cumsum(w / sum(w))
    cw[length(cw)] <- 1
    u <- (stats::runif(1) + seq_along(j) - 1) / length(j)
    idx[j] <- j[findInterval(u, cw) + 1L]
  }
  list(idx = idx, k = k)
}

#' @noRd
draw_reweight_say <- function(object, chk, k, ndraws, new_levels = FALSE) {
  head <- paste0(draw_check_said_what(chk), ", so the draws are reweighted ",
                 "by likelihood over proposal density (Pareto-smoothed ",
                 "importance sampling)")
  if (!is.finite(k) || k > draw_pareto_tol) {
    frm_warning(head, ", but ",
                if (is.finite(k)) {
                  paste0("the weights' Pareto k is ", format(k, digits = 3),
                         ", above ", draw_pareto_tol, ": no normal law ",
                         "reaches this likelihood, and")
                } else {
                  paste0(ndraws, " draws are too few to judge the weights,",
                         " and")
                },
                " the reweighted draws are not reliable. ",
                "frmtmb.sample::frm_sample() draws from the posterior itself",
                if (new_levels) {
                  paste0(", but its posterior_predict() cannot yet ",
                         "predict a grouping level the fit did not see")
                },
                class = "frmtmb_draw_pareto", call. = FALSE)
    return(invisible(NULL))
  }
  cache <- object$cache
  if (is.environment(cache) && isTRUE(cache[["draw_check_said"]])) {
    return(invisible(NULL))
  }
  if (is.environment(cache)) cache[["draw_check_said"]] <- TRUE
  frm_message(head, ", Pareto k ", format(k, digits = 2),
              ". Said once per fit", class = "frmtmb_draw_reweighted")
  invisible(NULL)
}

#' Are the parameters the standard-error check calls "limit" the family's
#' own limiting case (a zero inflation at 0, a precision at the binomial
#' limit), which keeps the other standard errors? The check's own
#' predicate decides (se_limit_quiet(), from the boundary rule for family
#' parameters).
#'
#' @noRd
draw_limit_quiet <- function(object, lost) {
  isTRUE(se_limit_quiet(object, lost))
}

#' Which "limit" parameters are a family parameter at the END of its
#' range in a family limiting case (se_limit_draw_held(), the check's own
#' predicate), named by `names(lost)`. A family parameter the limit leaves
#' undefined (a `nu ~ x` slope at the gaussian limit) is not one: it keeps
#' its variance 0, as the check holds it.
#'
#' @noRd
draw_limit_held <- function(object, lost) {
  out <- stats::setNames(logical(length(lost)), names(lost))
  v <- se_limit_draw_held(object, lost)
  out[names(v)[v %in% TRUE]] <- TRUE
  out
}

#' A family parameter at the end of its range (draw_limit_held()) is held
#' in `cov_fixed_prop`, variance 0, which would drop the side the data do
#' identify. Its row of sdreport()'s own inverse is put back, where that
#' is finite, so that draw_edge_law() decides each side of it as it does
#' without that check. The caller has already sent any "limit" that is not
#' the family's limiting case to the fallback.
#'
#' @noRd
draw_limit_restore <- function(V, ds, object) {
  lim <- names(ds$lost)[ds$lost == "limit"]
  raw <- ds$raw
  if (is.null(V) || !length(lim) || is.null(raw) ||
        !all(dim(raw) == dim(V))) {
    return(V)
  }
  held <- draw_limit_held(object, ds$lost)
  lim <- intersect(lim, names(held)[held])
  j <- match(lim, ds$map$names)
  j <- j[!is.na(j)]
  if (!length(j) || !all(is.finite(raw[j, ])) ||
        !all(diag(raw)[j] > 0)) {
    return(V)
  }
  V[j, ] <- raw[j, ]
  V[, j] <- raw[, j]
  V
}

#' How many standard deviations draw_edge_side() moves a parameter for
#' its flat test, and the normal law's fall there, `k^2 / 2`: the
#' distance about one draw in twenty passes.
#'
#' @noRd
draw_edge_k <- 2

#' How far, as a share of `draw_edge_k` conditional standard deviations,
#' a side has to reach the normal law's fall before it is narrowed. A
#' side the normal law describes reaches it at 1, a skewed one somewhat
#' sooner; 0.5 narrows only a side that falls at least four times
#' faster than the law says (dev/ppdraw-findings.md, "Fix round").
#'
#' @noRd
draw_edge_c <- 0.5

#' Below this ratio a steep side is narrowed in full (draw_edge_blend()).
#'
#' @noRd
draw_edge_r0 <- 0.4

#' The largest change in log-likelihood draw_edge_side() calls flat.
#' Fourteen estimates at the end of their range change by 1.3e-8 to
#' 1.1e-7 on their flat side (dev/ppdraw-probe6.R); interior estimates
#' near an end, from about 0.003 (a group log sd lme4 does not call
#' singular). Between the two, the gain over the limit is what decides,
#' and below 1e-4 holding that side is the honest reading of it.
#'
#' @noRd
draw_edge_flat <- 1e-4

#' How many test draws draw_check() reads, and from which seed.
#' (dev/predraw2-findings.md has the measurement behind 20.)
#'
#' @noRd
draw_check_n <- 20L

#' @noRd
draw_check_seed <- 20261113L

#' The median excess drop (draw_logratio()) above which draw_check()
#' fires, the first of its two conditions (the second is
#' `draw_check_rel`; draw_check_verdict()). On 430 fits (310 GLMs and
#' GLMMs, and zero-inflated, hurdle and mixture fits at and away from an
#' end) the median over 20 test draws was at most 0.29; on the three
#' weakly identified mixtures of backlog 4882, 1.07 to 9.0. Across the
#' fix round's 718 verdicts the largest ratio of an unfired fit outside
#' the mixtures was 0.116 (dev/predraw2-findings.md, section 11).
#'
#' @noRd
draw_check_tol <- 1

#' The median excess over the median drop the law implies above which
#' draw_check() fires (with `draw_check_tol`). On GLMs with one
#' coefficient per level of a factor and 2 to 20 rows per level
#' (poisson, bernoulli, gamma, gaussian; 10 to 500 levels), the excess
#' reached 7.4 at 500 levels while this ratio stayed at most 0.12 (0.066
#' where the excess passed 1); on the fired weakly identified mixtures
#' it was 0.39 or more (dev/predraw2-findings.md, fix round).
#'
#' @noRd
draw_check_rel <- 0.15

#' The excess (draw_check_dirs()) above which one direction of the draw
#' law fires the check, a heuristic borrowed from the side law's
#' `draw_edge_c`: for a quadratic misfit it is a law whose standard
#' deviation along the direction is twice the likelihood's (a drop of 8
#' at `draw_edge_k` standard deviations where the law says 2); for any
#' other shape it is a density error of exp(6) there. It is not a
#' calibrated false-alarm rate. It was set between the largest reading
#' of the 472 gaussian GLMs and GLMMs it was chosen on (5.6, leaving
#' out four random-slope fits at 7.0 to 46 whose normal law gave a new
#' level's 98% interval 1.36 to 2.04 times as wide as the law with the
#' parameters known) and the smallest of the 18 weakly identified
#' mixtures the summed check fires on (13.3). On 160 fits of other
#' designs without a random slope it fired on one (an identity-link
#' sigma that the normal law drew below 0), on 126 random-slope fits on
#' 8 (dev/drawdir-findings.md).
#'
#' @noRd
draw_dir_tol <- draw_edge_k^2 / 2 * (1 / draw_edge_c^2 - 1)

#' The shares of the distance at which draw_check_dirs() reads a
#' direction again where the objective cannot be evaluated at its end.
#'
#' @noRd
draw_dir_back <- c(0.75, 0.5)

#' The template components whose outer parameters are covariance-
#' structure parameters, on their unconstrained scale. `theta` holds the
#' parameters of every random-effect block: log sds and correlations,
#' and for structured covariances their ranges, Matern smoothness,
#' decay rates, ar correlations and `rr()` loadings (frmtmb.spatial's
#' `car` and `spde` included), the `gp()` and smooth scales, and the
#' latent-residual sd and its autocorrelation; `thetaac` the
#' autocorrelation parameters; `thetar` the residual correlations;
#' `logsdme` and `thetame` the `me()` latent sds and correlations
#' (`R/frame.R`, where the template is built). A variance in them can
#' have a side toward 0 along which the likelihood flattens, which the
#' reweighting cannot weigh; the other parameters of the set ride with
#' the same covariance. `beta`,
#' `betad` (the coefficients of every linear predictor, a sigma's and a
#' mixture weight's among them) and `meanme` are not: they are
#' regression coefficients, and on slpw:12 of the review the
#' reweighting along a sigma direction mended a new level's interval
#' the normal law had put at 162 times the posterior's width.
#'
#' @noRd
draw_vc_comps <- c("theta", "thetaac", "thetar", "logsdme", "thetame")

#' The Pareto k above which reweighted draws are not trusted (loo's
#' threshold for importance sampling).
#'
#' @noRd
draw_pareto_tol <- 0.7

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
#' Only the blocks `re_formula = NA` would drop are written back, and a
#' `gp()` field's coefficients at every `re_formula`; a population smooth
#' is held at its mode as before. With a `gp()` field the replicate also
#' draws the field at a position the fit did not see given its own field
#' (`krige_draw`).
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
  # a gp() field is a population term, which re_formula does not drop: its
  # coefficients are drawn at every re_formula, as brms's draws carry them
  gov <- sort(unique(c(if (re_form_keeps(re_formula)) re_governed_b(object),
                       gp_drawn_b(object))))
  if (!length(gov)) return(NULL)
  Q <- joint_precision(object)
  rn <- rownames(Q)
  if (is.null(Q) || !any(rn == "b")) {
    frm_warning("predict(): the fitted objective carries no random-effect ",
                "block in its joint precision (quadrature = TRUE or ",
                "importance = marginalizes it), so the group effects at the ",
                "levels the ",
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
  gpx <- gp_field_draw_plan(object, gov)
  bseeds <- sample.int(.Machine$integer.max, ndraws)
  function(fs, s) {
    set.seed(bseeds[s])
    x <- Matrix::solve(ch, Matrix::solve(ch, stats::rnorm(nr), system = "Lt"),
                       system = "Pt")
    x <- as.numeric(x)
    xn <- x[sel]
    if (!is.null(delta)) {
      x <- x - as.numeric(Matrix::solve(ch, Qrd %*% delta(s),
                                        system = "A"))
    }
    fs$estimates[["b"]][gov] <- b0 + x[sel]
    if (!is.null(gpx) && !is.null(delta)) fs <- gp_field_draw(fs, gpx, xn)
    # at a position the fit did not see, the field's value given the
    # replicate's own field, as a posterior draw of brms makes it
    fs[["krige_draw"]] <- TRUE
    fs
  }
}

#' Positions in `b` of every `gp()` block, exact and approximate: the
#' coefficients `predict()` draws whatever `re_formula` says.
#'
#' @noRd
gp_drawn_b <- function(fit) {
  sort(unlist(lapply(fit$frame[["re_blocks"]] %||% list(), function(bk) {
    if (bk[["covstruct"]] %in% c("gp", "hsgp")) bk[["b_idx"]]
  })))
}

#' What gp_field_draw() needs of the exact gp() blocks among the drawn
#' coefficients `gov`: where their z sit in the objective's full vector
#' and in `gov`. `NULL` where the replicate's parameters are not the
#' objective's outer vector: under `REML = TRUE` or `profile = TRUE` beta
#' is inner and drawn as well, and under quadrature or importance
#' sampling there is no inner mode. Those keep the first-order draw.
#'
#' @noRd
gp_field_draw_plan <- function(fit, gov) {
  if (isTRUE(fit$REML) || isTRUE(fit$control$profile) ||
        isTRUE(fit$quadrature) || !is.null(fit$importance) ||
        !is.null(fit$beta_whiten) || is.null(fit$obj)) {
    return(NULL)
  }
  bnames <- which(names(fit$obj$env$last.par) == "b")
  blocks <- list()
  for (bk in fit$frame[["re_blocks"]] %||% list()) {
    if (!identical(bk[["covstruct"]], "gp")) next
    at <- match(bk[["b_idx"]], gov)
    if (anyNA(at)) next
    blocks[[length(blocks) + 1L]] <- list(b_idx = bk[["b_idx"]], at = at,
                                          pos = bnames[bk[["b_idx"]]])
  }
  if (!length(blocks) || anyNA(unlist(lapply(blocks, `[[`, "pos")))) {
    return(NULL)
  }
  list(obj = fit$obj, blocks = blocks,
       newton = tryCatch(gp_field_newton_plan(fit), error = function(e) NULL))
}

#' What gp_field_newton() needs, or `NULL` where it does not apply: one exact
#' gp() block, the model's only random effect, in a linear predictor that is
#' no nonlinear body's, each row reading at most one of its positions. The
#' rows that read one (`ri`), the position (`pj`, sorted) and the row's
#' multiplier (`pm`), and the field at the estimates there (`c0`).
#'
#' @noRd
gp_field_newton_plan <- function(fit) {
  fr <- fit$frame
  bks <- fr[["re_blocks"]] %||% list()
  if (length(bks) != 1L || !identical(bks[[1L]][["covstruct"]], "gp") ||
        length(fit$spec$responses) != 1L || length(fr[["autocor"]]) ||
        length(fr[["me"]]$terms) || length(fr[["par_template"]][["miss"]])) {
    return(NULL)
  }
  bk <- bks[[1L]]
  key <- gp_block_lp(bk, fr)
  lp <- fr[["linpreds"]][[key]]
  if (is.null(lp) || !is.null(lp[["nl_body"]]) ||
        !identical(lp[["par"]], "beta")) {
    return(NULL)
  }
  Zb <- methods::as(lp_Z(fit, lp)[, bk[["c_idx"]], drop = FALSE],
                    "TsparseMatrix")
  nz <- Matrix::summary(Zb)
  nz <- nz[nz$x != 0, , drop = FALSE]
  if (!nrow(nz) || any(tabulate(nz$i, nrow(Zb)) > 1L)) return(NULL)
  est <- fit$estimates
  r <- lp[["resp"]]
  list(bk = bk, key = key, lp = lp, resp = r, ri = nz$i, pj = nz$j,
       # the log density quadratic in the field: one step is the mode
       quadratic = gp_laplace_exact(fit, list(bk)),
       pm = nz$x, m = ncol(Zb),
       c0 = gp_field(est[["b"]][bk[["b_idx"]]],
                     est[["theta"]][bk[["theta_idx"]]], bk),
       fam = fit$spec$responses[[r]]$family, y = fr[["y"]][[r]],
       av = fr[["aterm_values"]][[r]])
}

#' The field's inner mode at a replicate's own parameters by Newton steps
#' from the fit's field `c0`, on numbers: with `g` and `w` each position's
#' gradient and curvature of the log-likelihood in the field at `f`,
#' `f <- C (C + W^-1)^-1 (f + g / w)`, `C` the replicate's kernel, until the
#' step is under 1e-2 of the field. One step is the mode where the log
#' density is quadratic in the field (a gaussian mean with the identity
#' link). A factor of the kernel per step, where the tape's inner solve took
#' 3 to 21 times the replicate's whole cost
#' (dev/reviews/2026-11-13-gpfield.md, F2). `NULL` where a curvature is not
#' positive or the steps do not settle: the tape's solve decides then.
#'
#' @noRd
gp_field_newton <- function(fs, nt) {
  bk <- nt$bk
  th <- fs$estimates[["theta"]][bk[["theta_idx"]]]
  dp <- eval_dpars(fs, b = NULL)[[nt$resp]]
  dpar <- nt$lp[["dpar"]]
  lk <- nt$lp[["link"]]
  eta <- lk$linkfun(as.numeric(dp[[dpar]]))
  ex <- fit_extras(fs, nt$resp)
  wc <- nt$av[["weights"]] %||% 1
  K <- gp_corr(th, bk)
  sd2 <- sd_floored(th[1L])^2
  C <- sd2 * K
  dg <- cbind(seq_len(nt$m), seq_len(nt$m))
  f <- nt$c0
  done <- FALSE
  for (it in seq_len(if (nt$quadratic) 1L else gp_field_newton_steps)) {
    e <- eta
    e[nt$ri] <- e[nt$ri] + nt$pm * f[nt$pj]
    h <- 1e-3 * pmax(1, abs(e))
    ll <- function(k) {
      dp[[dpar]] <- lk$linkinv(e + k * h)
      wc * as.numeric(row_lpdf(nt$fam, nt$y, nt$y, dp, nt$av, ex))
    }
    l1 <- ll(-1)
    l2 <- ll(0)
    l3 <- ll(1)
    g <- ((l3 - l1) / (2 * h))[nt$ri]
    w <- (-(l3 - 2 * l2 + l1) / h^2)[nt$ri]
    gp <- as.numeric(rowsum(nt$pm * g, nt$pj, reorder = TRUE))
    wp <- as.numeric(rowsum(nt$pm^2 * w, nt$pj, reorder = TRUE))
    if (length(wp) != nt$m || !all(is.finite(c(gp, wp))) || any(wp <= 0)) {
      return(NULL)
    }
    M <- C
    M[dg] <- M[dg] + 1 / wp
    R <- tryCatch(chol(M), error = function(e) NULL)
    if (is.null(R)) return(NULL)
    f1 <- as.vector(C %*% backsolve(R, backsolve(R, f + gp / wp,
                                                  transpose = TRUE)))
    step <- max(abs(f1 - f))
    f <- f1
    # Newton converges quadratically: a step of 1e-2 of the field leaves
    # the field it reaches about 1e-4 of it off
    if (nt$quadratic || step <= 1e-2 * max(1, abs(f))) {
      done <- TRUE
      break
    }
  }
  if (!done) return(NULL)
  p <- bk[["gp_order"]]
  # gp_factor_ordered()'s own arithmetic, so the field the replicate reads
  # through the cached factor is the one it would have factored
  L <- t(gp_chol_num(K[p, p, drop = FALSE]))
  list(z = forwardsolve(L, f[p]) / sd_floored(th[1L]),
       cache = list(theta = as.numeric(th), L = L))
}

# How many Newton steps gp_field_newton() takes toward a replicate's mode
# before the tape's own solve decides. From the fit's field a poisson
# field at 300 positions settles in a few; one step alone left the field
# up to 1.7 posterior sds off (dev/gpfield-findings.md, fix round F2).
gp_field_newton_steps <- 8L

#' One replicate's exact gp() field: z at the inner mode of the replicate's
#' own parameters, plus the drawn deviation `xn` from it. The first-order
#' shift of the mode, which the group effects take, is not enough for a
#' field: z is the field over `sd L(theta)`, so where the data pin the
#' field the mode of z moves with theta along a curve. Read through each
#' replicate's sd and length scale, the first-order z spread the field at
#' a fitted position 1.38 where its noise is 0.2, and the field it gave was
#' not one that kernel draws: kriged 5 length scales past the positions it
#' spread 52 where brms reads 3.4 (dev/gpfield-findings.md, 1919).
#'
#' @noRd
gp_field_draw <- function(fs, plan, xn) {
  if (!is.null(plan$newton)) {
    zn <- tryCatch(gp_field_newton(fs, plan$newton), error = function(e) NULL)
    if (!is.null(zn) && all(is.finite(zn$z))) {
      g <- plan$blocks[[1L]]
      fs$estimates[["b"]][g$b_idx] <- zn$z + xn[g$at]
      # the factor at this replicate's theta, which its field reads next
      fs$frame[["re_blocks"]][[1L]][["gp_factor_cache"]] <- zn$cache
      return(fs)
    }
  }
  obj <- plan$obj
  env <- obj$env
  par <- outer_from_template(fs$estimates, obj, fs$frame)
  if (is.null(par)) return(fs)
  saved <- obj_state_save(obj)
  on.exit(obj_state_restore(obj, saved), add = TRUE)
  v <- tryCatch(suppressWarnings(obj$fn(par)), error = function(e) NA_real_)
  if (!isTRUE(is.finite(v))) return(fs)
  for (g in plan$blocks) {
    fs$estimates[["b"]][g$b_idx] <- env$last.par[g$pos] + xn[g$at]
  }
  fs
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
      error = keep_missing(NULL))
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
        error = keep_missing(NULL))
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
    # an exact gp()'s z reaches every position, so no batch can hold it, and
    # fd_eta_chain() takes it where a design exists; the rest of the fit keeps
    # its batches
    if (identical(bk[["covstruct"]], "gp")) next
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
    # The WHOLE kept set as ONE batch, when no row loads two of its
    # columns. The per-position split below is what a `(1 + x | g)` term
    # needs, but it costs one batch per POSITION, and for a smooth,
    # `gp()` or `hsgp()` block `dim` IS the coefficient count, so there
    # the split degenerates to one batch per coefficient and buys
    # nothing. An exact `gp()` in sample loads one column per row, its
    # design being an indicator, so the whole block is one exact batch:
    # 329 model evaluations become 11 and `fitted()` goes from 0.63 s to
    # 0.02 s on a 160-coefficient cumulative fit
    # (dev/resmooth-batchcost-before.txt against -after.txt).
    # A `(1 | g)` block has one position, so it reaches the same single
    # batch by a shorter route; a `(1 + x | g)` block FAILS this test,
    # because a row loads its level's intercept and its slope, and falls
    # through to the split unchanged.
    whole <- re_batch_try(S, ci, bi, which(bi %in% b_idx))
    if (!is.null(whole)) {
      out[[length(out) + 1L]] <- whole
      next
    }
    # one batch per COLUMN POSITION of the block: a row of a
    # `(1 + x | g)` term loads its level's intercept AND its slope, so a
    # batch over the whole block could not attribute the difference,
    # while a batch over one position holds one nonzero per row. An |ID|
    # block spanning two predictors splits the same way, because each
    # predictor is its own position.
    for (k in seq_len(D)) {
      at <- seq.int(k, length(ci), by = D)
      keep <- at[bi[at] %in% b_idx]
      if (!length(keep)) next
      one <- re_batch_try(S, ci, bi, keep)
      # a row that loads two levels at this position cannot be
      # attributed (a multi-membership term does exactly that)
      if (is.null(one)) return(NULL)
      out[[length(out) + 1L]] <- one
    }
  }
  out
}

#' One batch over `keep`, positions within a block's `c_idx`, or `NULL`
#' when a row loads two of those columns, so that the difference read off
#' that row could not be attributed to one of them.
#'
#' @noRd
re_batch_try <- function(S, ci, bi, keep) {
  if (!length(keep)) return(NULL)
  Sb <- S[, ci[keep], drop = FALSE]
  if (max(Matrix::rowSums(Sb)) > 1) return(NULL)
  owner <- rep(NA_integer_, nrow(Sb))
  nz <- Matrix::which(Sb != 0, arr.ind = TRUE)
  owner[nz[, 1L]] <- nz[, 2L]
  list(idx = bi[keep], owner = owner)
}

#' The chain-rule columns `fit_fd_se()` can take for the smooth, `gp()`
#' and `hsgp()` blocks among `b_idx`, and the variance that is not
#' coefficient uncertainty on each predictor of `resp`.
#'
#' A row of such a block's design loads every one of its columns once
#' the row is off the fitted positions (a kriging row, a smooth basis
#' row), so no batch can attribute a difference, and one pair of
#' evaluations per coefficient was the whole cost: 329 evaluations for a
#' 160-coefficient exact `gp()` at three new rows. Its effects reach the
#' prediction only through one predictor's `eta`, linearly through `Z`,
#' so the derivative is `(d out / d eta) Z`: one pair per predictor.
#' `Z` is the design `lp_eta_design()` builds, in sample the predictor's
#' own `Z`.
#'
#' Returns `list(chain, extra)`, each a list of entries keyed by the
#' predictor; empty where the design cannot be rebuilt, which leaves the
#' caller on the routes it had.
#'
#' @noRd
fd_eta_chain <- function(fit, newdata, resp, allow_new_levels, use_re,
                         b_idx) {
  chain <- list()
  extra <- list()
  for (lp in fit$frame[["linpreds"]]) {
    if (!is.null(resp) && !identical(lp[["resp"]], resp)) next
    if (!is.null(lp[["nl_body"]])) next
    key <- linpred_key(lp[["resp"]], lp[["dpar"]])
    ed <- tryCatch(suppressWarnings(
      lp_eta_design(fit, lp, newdata, use_re, allow_new_levels)),
      error = keep_missing(NULL))
    if (is.null(ed)) next
    for (bk in fit$frame[["re_blocks"]] %||% list()) {
      if (!bk[["covstruct"]] %in% c("smooth", "gp", "hsgp") ||
          (!block_b_positionwise(bk) && bk[["covstruct"]] != "gp")) {
        next
      }
      mine <- vapply(bk[["components"]], function(cp) {
        identical(cp$lp_key, key)
      }, NA)
      if (!any(mine)) next
      keep <- which(bk[["b_idx"]] %in% b_idx)
      if (!length(keep)) next
      Z <- if (is.null(newdata)) {
        if (is.null(lp[["Z"]])) next
        as.matrix(lp[["Z"]][, bk[["c_idx"]][keep], drop = FALSE])
      } else {
        sp <- Find(function(s) identical(s$bk[["c_idx"]], bk[["c_idx"]]),
                   ed[["sm_parts"]])
        if (is.null(sp)) next
        as.matrix(sp$Xr)[, keep, drop = FALSE]
      }
      # an exact gp()'s design is over its field, and b is z: linear in z
      # at a theta, so the chain still holds through d cvec / d b
      if (bk[["covstruct"]] == "gp") {
        Js <- rr_jacobians(fit)$Jb[bk[["c_idx"]], bk[["b_idx"]], drop = FALSE]
        Zc <- if (is.null(newdata)) {
          lp[["Z"]][, bk[["c_idx"]], drop = FALSE]
        } else {
          sp$Xr
        }
        Z <- as.matrix(Zc %*% Js[, keep, drop = FALSE])
      }
      chain[[length(chain) + 1L]] <- list(key = key,
                                          idx = bk[["b_idx"]][keep], Z = Z)
    }
    if (!is.null(newdata)) {
      ev <- lp_extra_var_vec(fit, ed, use_re)
      if (any(ev > 0)) {
        extra[[length(extra) + 1L]] <- list(key = key, var = ev)
      }
    }
  }
  list(chain = chain, extra = extra)
}

#' Positions in `b` of EVERY smooth, `gp()` and `hsgp()` block.
#'
#' `re_formula` does not govern these: `NA` keeps them, so the
#' finite-difference standard error has to difference them at every
#' `re_formula` or it reports the uncertainty of a smaller model than the
#' one the estimate came from. Measured before this was added: on a
#' cumulative fit with `s(x, g, bs = "fs")`, `Est.Error` at
#' `re_formula = NA` was short of a reference that differences the
#' smooth's 30 coefficients by 0.52 of itself at the median row and 0.68
#' at the worst (`dev/resmooth-fdse.txt`).
#'
#' @noRd
smooth_b_idx <- function(fit) {
  ids <- integer(0)
  for (bk in fit$frame[["re_blocks"]]) {
    if (bk[["covstruct"]] %in% c("smooth", "gp", "hsgp", "latres")) {
      ids <- c(ids, bk[["b_idx"]])
    }
  }
  sort(unique(ids))
}

#' Positions in `b` the finite-difference standard error has to
#' difference when `re_formula` KEEPS the group-level effects: every
#' block with a component left in the (possibly reduced) design, and a
#' smooth whose basis is indexed by a grouping factor.
#'
#' The group-indexed smooth is here for the standard error's sake, not
#' because `re_formula` drops it: it does not, any more than it drops
#' `s(x)` (`lp_eta_design()`). It stays because leaving its columns out
#' of the differenced set put `Est.Error` 68 percent wrong on a
#' cumulative fit (`dev/reunc-log/fdsmooth-prefix.txt`).
#'
#' A population smooth, `gp()` and `hsgp()` are deliberately NOT here,
#' and their absence is no longer a gap: `smooth_b_idx()` carries every
#' smooth block at every `re_formula`, and each caller unions the two.
#' What this function adds on top is the `(x | g)` group-level part,
#' which is the part `re_formula` governs. See dev/resmooth-findings.md.
#'
#' @noRd
re_governed_b <- function(fit) {
  blocks <- fit$frame[["re_blocks"]]
  gs <- unlist(lapply(fit$frame[["linpreds"]], smooth_group_block_ids))
  ids <- integer(0)
  for (i in seq_along(blocks)) {
    bk <- blocks[[i]]
    if (bk[["covstruct"]] %in% c("smooth", "gp", "hsgp", "latres")) {
      if (i %in% gs) ids <- c(ids, bk[["b_idx"]])
      next
    }
    if (length(bk[["components"]])) ids <- c(ids, bk[["b_idx"]])
  }
  sort(unique(ids))
}
