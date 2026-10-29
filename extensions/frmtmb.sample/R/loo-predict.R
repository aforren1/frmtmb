# brms's PSIS-weighted functions on frm_sample() draws: loo_epred(),
# loo_linpred(), loo_predict(), loo_predictive_interval(), loo_R2() and
# loo_model_weights(). Each one is a function of the log_lik() matrix and
# of one prediction matrix, so each takes ONE log_lik() and ONE
# loo::psis() per call (per response for loo_R2(), per model for
# loo_model_weights()) and hands them to the loo package. The generics
# and their refusals on a fit are frmtmb's (R/loo-generics.R there).

#' Leave-one-out expectations of draws
#'
#' brms's `loo_epred()`, `loo_linpred()`, `loo_predict()` and
#' `loo_predictive_interval()` on draws from [frm_sample()]. Each takes
#' the matching prediction ([posterior_epred()], [posterior_linpred()]
#' or [posterior_predict()]) and summarizes every observation's column
#' under that observation's leave-one-out importance weights, by
#' `loo::E_loo()`. The weights are the Pareto-smoothed ones of
#' [psis()]: one [frmtmb::log_lik()] matrix and one `loo::psis()` per
#' call, on the same draws as the prediction.
#'
#' @section Leave-one-out, and newdata:
#' The weight of a draw for observation `i` comes from that observation's
#' own column of [frmtmb::log_lik()], so a leave-one-out expectation exists
#' only for a row the model was fitted on. brms passes `newdata` to both
#' `loo()` and the prediction, and so computes the weights on the new
#' rows' own response. Here [frmtmb::loo()] and [frmtmb::log_lik()] refuse
#' `newdata`, and these functions refuse it in the same words, with
#' `re_formula` and the new-level arguments. For a summary on other
#' rows, use [posterior_predict()] or [posterior_epred()] with
#' `newdata`: on a row the model has not seen, there is no observation
#' to leave out.
#'
#' @section Draws:
#' `ndraws` and `draw_ids` choose the draws ONCE, and the same draws give
#' the weights and the prediction. brms draws a random subsample for
#' each of the two, so its weights then belong to other draws than the
#' prediction they weight. As in [frmtmb::loo()] and brms, the relative
#' efficiency is taken on the chain structure when every draw is used,
#' and on one chain for a subset (`ndraws` or `draw_ids`). Draws from
#' `frm_sample(thin =)` keep their chains.
#'
#' @section Multivariate models:
#' Without `resp`, the weights come from the joint [frmtmb::log_lik()]
#' of the responses, and the result has a third dimension with one layer
#' per response, named as brms names the responses (`ya` for a column
#' `y_a`). With one response in `resp`, the weights come from that
#' response's own log-likelihood, as brms's `loo(resp =)` computes
#' them, and the result is that response's matrix. With several, the
#' weights come from [frmtmb::log_lik()] of those responses (their
#' joint density under `set_rescor(TRUE)`, the sum of their densities
#' without it), and the result has one layer per response named, in
#' that order, as in brms.
#'
#' @param object A `frmtmb_draws` from [frm_sample()].
#' @param type The summary of each column: `"mean"` (the default),
#'   `"var"` or `"quantile"`, as in brms.
#' @param probs The quantiles for `type = "quantile"`.
#' @param prob For `loo_predictive_interval()`, the mass of the central
#'   interval, one number in `(0, 1)`.
#' @param psis_object A `psis` object from `loo::psis()`, [psis()] or
#'   `loo(save_psis = TRUE)$psis_object`, for the same draws and rows.
#'   Given, it is used as it is and no [frmtmb::log_lik()] is computed.
#' @param resp For a multivariate model, the responses to summarize, in
#'   brms's spelling (`ya`) or frmtmb's (`y_a`). `NULL` is every
#'   response.
#' @param ... brms's further arguments. `ndraws` and `draw_ids` choose
#'   the draws for both halves. The arguments of the prediction
#'   function (`dpar`, `nlpar`, `transform`, `incl_thres`, `sort`)
#'   reach it. `r_eff` and `cores` reach `loo::psis()`. `moment_match`
#'   and `reloo` are refused when `TRUE`, as [frmtmb::loo()] refuses
#'   them. `newdata`, `re_formula` and the new-level arguments are
#'   refused unless they hold brms's default value (`NULL`;
#'   `allow_new_levels = FALSE`; `sample_new_levels = "uncertainty"`),
#'   which asks for nothing (see *Leave-one-out, and newdata*). Any
#'   other name is an error that names it.
#' @return A matrix with one row per observation and one column per
#'   summary, named `mean`, `var` or `q<100 p>`. Where the prediction
#'   has a third dimension (categories, or the responses of a
#'   multivariate model without `resp` or with several), an array
#'   with that dimension third.
#' @seealso [psis()], [frmtmb::loo()], [loo_R2()]
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE) &&
#'     requireNamespace("loo", quietly = TRUE) &&
#'     !frmtmb.sample:::tmbstan_build_broken()) {
#'   set.seed(9)
#'   dd <- data.frame(x = rnorm(60))
#'   dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
#'   ds <- frm_sample(bf(y ~ x), family = gaussian(), data = dd,
#'                    chains = 2, iter = 500, refresh = 0)
#'   head(loo_epred(ds))
#'   head(loo_predictive_interval(ds, prob = 0.8))
#'
#'   # one set of weights, reused for several summaries
#'   ps <- psis(ds)
#'   head(loo_predict(ds, type = "var", psis_object = ps))
#' }
#' }
#' @name sample-loo_predict
NULL

#' @rdname sample-loo_predict
#' @exportS3Method brms::loo_epred
#' @rawNamespace S3method(rstantools::loo_epred,frmtmb_draws)
#' @export
loo_epred.frmtmb_draws <- function(object, type = c("mean", "var",
                                                    "quantile"),
                                   probs = 0.5, psis_object = NULL,
                                   resp = NULL, ...) {
  type <- frm_match_arg(type)
  loo_draws_expect(object, posterior_epred, "loo_epred()", type, probs,
                   psis_object, resp, list(...),
                   c("dpar", "nlpar", "sort"))
}

#' @rdname sample-loo_predict
#' @exportS3Method rstantools::loo_linpred
#' @export
loo_linpred.frmtmb_draws <- function(object, type = c("mean", "var",
                                                      "quantile"),
                                     probs = 0.5, psis_object = NULL,
                                     resp = NULL, ...) {
  type <- frm_match_arg(type)
  loo_draws_expect(object, posterior_linpred, "loo_linpred()", type,
                   probs, psis_object, resp, list(...),
                   c("transform", "dpar", "nlpar", "incl_thres", "sort"))
}

#' @rdname sample-loo_predict
#' @exportS3Method rstantools::loo_predict
#' @export
loo_predict.frmtmb_draws <- function(object, type = c("mean", "var",
                                                      "quantile"),
                                     probs = 0.5, psis_object = NULL,
                                     resp = NULL, ...) {
  type <- frm_match_arg(type)
  loo_draws_expect(object, posterior_predict, "loo_predict()", type,
                   probs, psis_object, resp, list(...),
                   c("transform", "negative_rt", "sort", "ntrys"))
}

#' @rdname sample-loo_predict
#' @exportS3Method rstantools::loo_predictive_interval
#' @export
loo_predictive_interval.frmtmb_draws <- function(object, prob = 0.9,
                                                 psis_object = NULL, ...) {
  if (length(prob) != 1L) {
    # brms's words
    frm_stop("Argument 'prob' should be of length 1.", call. = FALSE)
  }
  check_probability(prob, "prob")
  alpha <- (1 - prob) / 2
  draws_as_caller("loo_predictive_interval()", loo_predict(
    object, type = "quantile", probs = c(alpha, 1 - alpha),
    psis_object = psis_object, ...))
}

#' The arguments a `loo_*()` function on draws takes in `...`, sorted by
#' where each one goes.
#'
#' brms hands its dots to `loo()` AND to the prediction, and each of the
#' two ignores what it does not know. Here a name either reaches the
#' half that honors it or is refused, so a brms spelling that would
#' change nothing (a `scale` that brms 2.23.0's `posterior_linpred()`
#' reads nowhere) says so instead of changing nothing in silence.
#'
#' @noRd
loo_draws_opts <- function(dots, what, pred_args) {
  nms <- names(dots) %||% rep("", length(dots))
  if (any(!nzchar(nms))) {
    frm_stop(what, " takes its further arguments by name, and argument ",
             which(!nzchar(nms))[1L], " of `...` has none", call. = FALSE)
  }
  # newdata and the arguments that only mean something with it: the
  # weights are built from log_lik(), which refuses them, and the policy
  # is loo()'s, so its reason is given in its words. brms's own default
  # values ask for none of that and pass.
  given <- function(a) a %in% nms && !loo_default_value(a, dots[[a]])
  ll_refused <- c("newdata", "re_formula", "re.form", draws_new_level_args,
                  "nug")
  hit <- Filter(given, ll_refused)
  if (length(hit)) {
    frm_stop(what, " cannot honor `", hit[1L], "`: its weights are ",
             "leave-one-out importance ratios from log_lik() of the draws, ",
             "which takes the rows the model was fitted on and refuses `",
             hit[1L], "`. A row the model was not fitted on has no ",
             "observation to leave out. brms computes the weights on ",
             "newdata's own response; here, use posterior_predict(newdata ",
             "=) or posterior_epred(newdata =) for those rows",
             call. = FALSE)
  }
  for (a in intersect(nms, c("moment_match", "reloo"))) {
    check_flag(dots[[a]], a)
    if (dots[[a]]) {
      frm_stop(what, " cannot honor ", a, " = TRUE: brms refits the model ",
               "for the observations whose Pareto k is high, and draws ",
               "carry no Stan program to refit. Read the Pareto k values of ",
               "psis() instead", call. = FALSE)
    }
  }
  if ("pointwise" %in% nms) check_flag(dots$pointwise, "pointwise")
  psis_args <- c("r_eff", "cores")
  quiet <- c("ndraws", "draw_ids", "pointwise", "k_threshold",
             "moment_match", "reloo", "moment_match_args", "reloo_args",
             ll_refused)
  unknown <- setdiff(nms, c(psis_args, quiet, pred_args))
  if (length(unknown)) {
    frm_stop(what, " has no argument `", unknown[1L], "`. Its further ",
             "arguments are ndraws, draw_ids, ",
             paste(c(pred_args, psis_args), collapse = ", "),
             call. = FALSE)
  }
  list(ndraws = dots$ndraws, draw_ids = dots$draw_ids,
       pred = dots[intersect(nms, pred_args)],
       psis = dots[intersect(nms, psis_args)])
}

#' Whether `value` is brms's default for the argument `name`, so that
#' spelling it out asks for nothing: `FALSE` for the new-level switches,
#' `"uncertainty"` for `sample_new_levels`, `NULL` for the rest.
#'
#' @noRd
loo_default_value <- function(name, value) {
  switch(name,
    allow_new_levels = ,
    allow.new.levels = isFALSE(value),
    sample_new_levels = identical(value, "uncertainty"),
    is.null(value))
}

#' The response `resp` names, in frmtmb's spelling, from either
#' spelling: frmtmb's own (`y_a`) or brms's (`ya`). NULL stays NULL.
#'
#' @noRd
loo_resp_name <- function(fit, resp, what) {
  if (is.null(resp)) return(NULL)
  resps <- names(fit$spec$responses)
  stan <- brms_stan_name(resps)
  ok <- is.character(resp) && length(resp) == 1L && !is.na(resp)
  hit <- if (ok) {
    if (resp %in% resps) resp else resps[match(resp, stan)]
  }
  if (!ok || is.na(hit)) {
    frm_stop(what, ": `resp` must name one response of this model; it ",
             "has ", paste(unique(c(stan, resps)), collapse = ", "),
             call. = FALSE)
  }
  hit
}

#' The smoothed importance weights of draws for leaving out each column
#' of `ll`, one `loo::psis()` with the relative efficiency `loo()` uses.
#'
#' One function, so that [psis()], `loo(save_psis = TRUE)` and every
#' `loo_*()` function weight a draw the same way.
#'
#' @noRd
draws_psis_of <- function(ll, x, ...) {
  dots <- list(...)
  if (is.null(dots$r_eff)) dots$r_eff <- loo_r_eff(ll, x)
  # the ratios are the NEGATIVE log-likelihood, which is what
  # loo::loo.matrix() smooths and what brms's psis.brmsfit() passes on
  do.call(loo::psis, c(list(-ll), dots))
}

#' The PSIS object for the rows `rows` of draws `x`: log_lik() once, and
#' its smoothing once. `ll` comes back too, for loo_R2(), which needs the
#' ratios themselves.
#'
#' @noRd
loo_draws_weights <- function(x, rows, resp, what, psis_args = list()) {
  ll <- draws_as_caller(what, re_formula = FALSE,
                        loo_matrix(x, NULL, resp, what, draw_ids = rows))
  list(ll = ll, psis = do.call(draws_psis_of, c(list(ll, x), psis_args)))
}

#' The refusal when the loo package is missing, one sentence for every
#' function here that needs it.
#'
#' @noRd
loo_need_pkg <- function(what) {
  frm_stop(what, " needs the 'loo' package; install it", call. = FALSE)
}

#' Shared body of loo_epred(), loo_linpred() and loo_predict().
#'
#' @noRd
loo_draws_expect <- function(object, pred_fun, what, type, probs,
                             psis_object, resp, dots, pred_args) {
  # loo_predictive_interval() reaches here through loo_predict(), and a
  # refusal names the function the user called
  what <- draws_caller(what)
  if (!requireNamespace("loo", quietly = TRUE)) {
    loo_need_pkg(what)
  }
  op <- loo_draws_opts(dots, what, pred_args)
  if (identical(type, "quantile")) {
    ok <- is.numeric(probs) && length(probs) >= 1L && !anyNA(probs) &&
      all(probs > 0 & probs < 1)
    if (!ok) {
      frm_stop(what, ": `probs` must be numbers strictly between 0 and 1",
               call. = FALSE)
    }
  }
  fit <- draws_base_fit(object)
  resps <- names(fit$spec$responses)
  # several responses, as brms takes them: their layers, weighted by the
  # log_lik() of that set (the joint density under rescor). An empty
  # `resp` is every response, as brms's validate_resp() reads it
  if (!length(resp)) resp <- NULL
  resp <- if (length(resp) > 1L) {
    draws_resp_set(fit, resp)
  } else {
    loo_resp_name(fit, resp, what)
  }
  rows <- draws_subsample(object, op$ndraws, op$draw_ids)
  # the prediction first: its own refusals (a family with no simulator,
  # an unknown dpar) cost nothing, and the weights cost a log_lik()
  one <- function(r) {
    draws_as_caller(what, do.call(pred_fun, c(
      list(object, resp = r, draw_ids = rows), op$pred)))
  }
  pred <- if (length(resps) > 1L && length(resp) != 1L) {
    # brms's multivariate layout: one layer per response, named as brms
    # names them (y_a is ya). Each layer is drawn on its own: a
    # leave-one-out summary is per row and per response, so a joint
    # draw of a rescor model would give it the same expectation
    want <- resp %||% resps
    draws_resp_stack(fit, want, lapply(want, one), what)
  } else {
    one(resp)
  }
  if (is.null(psis_object)) {
    # brms's message
    frm_message("Running PSIS to compute weights")
    psis_object <- loo_draws_weights(object, rows, resp, what,
                                     op$psis)$psis
  }
  if (!inherits(psis_object, "psis") ||
        !identical(dim(psis_object), dim(pred)[1:2])) {
    frm_stop(what, ": `psis_object` must be a psis object for the ",
             dim(pred)[1L], " draws and ", dim(pred)[2L], " observations ",
             "this call predicts, as psis() of the same draws gives; it ",
             "has ", if (inherits(psis_object, "psis")) {
               paste(dim(psis_object), collapse = " x ")
             } else paste(class(psis_object), collapse = "/"),
             call. = FALSE)
  }
  loo_expect_value(pred, psis_object, type, probs)
}

#' brms's `E_loo_value()`: `loo::E_loo()` per column, observations as
#' rows and the summary as columns, a third dimension kept third.
#'
#' @noRd
loo_expect_value <- function(x, psis_object, type, probs) {
  labs <- if (identical(type, "quantile")) paste0("q", probs * 100) else type
  # loo's own name for the variance; brms's partial match reaches it
  ltype <- if (identical(type, "var")) "variance" else type
  one <- function(m) {
    v <- loo::E_loo(m, psis_object, type = ltype, probs = probs)$value
    # E_loo() gives a summaries by observations matrix for several
    # quantiles and a vector otherwise
    v <- if (is.matrix(v)) t(v) else matrix(v)
    colnames(v) <- labs
    v
  }
  if (length(dim(x)) != 3L) return(one(x))
  k <- dim(x)[3L]
  out <- array(NA_real_, c(dim(x)[2L], length(labs), k),
               dimnames = list(NULL, labs, dimnames(x)[[3L]]))
  for (j in seq_len(k)) out[, , j] <- one(x[, , j])
  out
}

#' LOO-adjusted R-squared of draws
#'
#' brms's `loo_R2()` on draws from [frm_sample()]: one minus the ratio
#' of the variance of the leave-one-out residuals to the variance of the
#' response, `1 - Var(y - E_loo(ypred)) / Var(y)`, with uncertainty from
#' the Bayesian bootstrap of Vehtari and Lampinen (2002). The
#' leave-one-out prediction of each row is its [posterior_epred()]
#' column under the row's Pareto-smoothed importance weights, as in
#' brms; each bootstrap draw reweights the observations by a flat
#' Dirichlet. Values are clipped to `[-1, 1]`, as brms clips them.
#'
#' The weights come from one [frmtmb::log_lik()] and one `loo::psis()` per
#' response, on the same draws as the prediction. With
#' `add_criterion(x, "loo_R2")` the bootstrap draws are stored, and a
#' later `loo_R2()` without further arguments reads them, as brms's
#' does.
#'
#' An ordinal or categorical family is refused, as [frmtmb::bayes_R2()]
#' refuses it: brms refuses the categorical families and treats ordinal
#' probabilities as a continuous prediction, with a warning. `newdata`
#' is refused for the reason [loo_epred()] gives: a row the model was
#' not fitted on has no observation to leave out.
#'
#' @param object A `frmtmb_draws` from [frm_sample()].
#' @param resp For a multivariate model, which responses, in brms's
#'   spelling (`ya` for a column `y_a`), as in brms and [frmtmb::bayes_R2()].
#'   `NULL`, the default, is every response, one `R2<resp>` column each.
#' @param summary If `TRUE` (the default), summarize the bootstrap draws
#'   into estimate, error and quantiles; if `FALSE`, return the
#'   `ndraws x responses` matrix of them.
#' @param robust If `TRUE`, the summary is the median and MAD.
#' @param probs Quantiles for the summary.
#' @param seed A seed for the Bayesian bootstrap. The random number
#'   state is restored afterwards, as in brms.
#' @param args_epred,args_loglik Named lists of arguments for
#'   [posterior_epred()] alone and for [frmtmb::log_lik()] alone.
#' @param ... `ndraws` or `draw_ids`, which choose the draws for both,
#'   and further arguments of [posterior_epred()].
#' @return A summary matrix with one row per response, or with
#'   `summary = FALSE` the matrix of bootstrap draws.
#' @references Vehtari, A. and Lampinen, J. (2002). Bayesian model
#'   assessment and comparison using cross-validation predictive
#'   densities. Neural Computation, 14(10), 2439-2468.
#' @seealso [frmtmb::bayes_R2()], [loo_epred()]
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE) &&
#'     requireNamespace("loo", quietly = TRUE) &&
#'     !frmtmb.sample:::tmbstan_build_broken()) {
#'   set.seed(9)
#'   dd <- data.frame(x = rnorm(60))
#'   dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
#'   ds <- frm_sample(bf(y ~ x), family = gaussian(), data = dd,
#'                    chains = 2, iter = 500, refresh = 0)
#'   loo_R2(ds, seed = 1)
#'   # in-sample R2 is larger: it scores rows the draws were fitted to
#'   bayes_R2(ds)
#' }
#' }
#' @name sample-loo_R2
NULL

#' @rdname sample-loo_R2
#' @exportS3Method rstantools::loo_R2
#' @export
loo_R2.frmtmb_draws <- function(object, resp = NULL, summary = TRUE,
                                robust = FALSE, probs = c(0.025, 0.975),
                                seed = NULL, args_epred = list(),
                                args_loglik = list(), ...) {
  check_flag(summary, "summary")
  check_flag(robust, "robust")
  if (!is.list(args_epred) || !is.list(args_loglik)) {
    frm_stop("loo_R2(): `args_epred` and `args_loglik` must be lists",
             call. = FALSE)
  }
  if (!requireNamespace("loo", quietly = TRUE)) {
    loo_need_pkg("loo_R2()")
  }
  fit <- draws_base_fit(object)
  resps <- names(fit$spec$responses)
  stan <- brms_stan_name(resps)
  sel <- brms_resp_index(resp, stan)
  lab <- if (length(resps) > 1L) paste0("R2", stan[sel]) else "R2"
  dots <- list(...)
  # brms reads a stored criterion whatever else is asked; here only when
  # nothing is asked that would change the draws behind it
  st <- object[["criteria"]][["loo_R2"]]
  if (is.matrix(st) && !length(dots) && !length(args_epred) &&
        !length(args_loglik) && all(lab %in% colnames(st))) {
    R2 <- st[, lab, drop = FALSE]
    return(if (summary) {
      posterior_summary(R2, probs = probs, robust = robust)
    } else R2)
  }
  op <- loo_draws_opts(dots, "loo_R2()", c("dpar", "nlpar"))
  for (a in c("ndraws", "draw_ids")) {
    if (a %in% c(names(args_epred), names(args_loglik))) {
      frm_stop("loo_R2(): give `", a, "` in `...`, where it chooses the ",
               "draws for both the prediction and the weights, not in ",
               "args_epred or args_loglik, where the two would weight ",
               "different draws", call. = FALSE)
    }
  }
  rows <- draws_subsample(object, op$ndraws, op$draw_ids)
  if (!is.null(seed)) {
    # brms's handling: seed the bootstrap, then put the caller's stream
    # back as it was
    if (exists(".Random.seed", envir = globalenv())) {
      old <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }
  R2 <- lapply(sel, function(r) {
    rn <- resps[r]
    y <- draws_observed_y(fit, rn)
    fam <- fit$spec$responses[[rn]]$family
    if (is.null(y) || is.matrix(y) ||
          isTRUE(fam[["type"]] %in% c("ordinal", "categorical")) ||
          !is.null(fam[["simplex"]])) {
      frm_stop("loo_R2() is not defined for the '", fam[["family"]],
               "' response '", rn, "': R-squared needs a numeric response ",
               "and a numeric prediction, and this family predicts a ",
               "category distribution. bayes_R2() refuses it for the same ",
               "reason. Use loo() to compare such models", call. = FALSE)
    }
    ep <- draws_as_caller("loo_R2()", do.call(posterior_epred, c(
      list(object, resp = rn, draw_ids = rows), op$pred, args_epred)))
    if (length(dim(ep)) > 2L) {
      frm_stop("loo_R2() is not defined for a prediction with several ",
               "values per row", call. = FALSE)
    }
    ll <- draws_as_caller("loo_R2()", re_formula = FALSE, do.call(
      loo_matrix, c(list(object, NULL, rn, "loo_R2()", draw_ids = rows),
                    args_loglik)))
    ps <- do.call(draws_psis_of, c(list(ll, object), op$psis))
    loo_R2_boot(as.numeric(y), ep, ll, ps)
  })
  R2 <- matrix(unlist(R2), ncol = length(sel), dimnames = list(NULL, lab))
  if (summary) posterior_summary(R2, probs = probs, robust = robust) else R2
}

#' brms's `.loo_R2()` on given matrices: the leave-one-out residuals'
#' variance against the response's, under `S` Bayesian bootstrap
#' weightings of the `N` observations.
#'
#' The `rexp()` draws are taken in brms's order, so one seed gives brms's
#' bootstrap. The weighted sums are matrix products rather than brms's
#' `rowSums(sweep())`, which builds three more `S x N` matrices; they
#' agree to rounding.
#'
#' @noRd
loo_R2_boot <- function(y, ypred, ll, psis_object) {
  ypredloo <- loo::E_loo(ypred, psis_object, log_ratios = -ll)$value
  err <- ypredloo - y
  S <- nrow(ypred)
  N <- ncol(ypred)
  w <- matrix(stats::rexp(S * N, rate = 1), nrow = S, ncol = N)
  w <- w / rowSums(w)
  m <- cbind(y, y^2, err, err^2)
  s <- w %*% m
  var_y <- (N / (N - 1)) * (s[, 2L] - s[, 1L]^2)
  var_err <- (N / (N - 1)) * (s[, 4L] - s[, 3L]^2)
  out <- unname(1 - var_err / var_y)
  out[out < -1] <- -1
  out[out > 1] <- 1
  out
}

#' Model weights by stacking or pseudo-BMA
#'
#' brms's `loo_model_weights()` on draws from [frm_sample()]: the
#' [frmtmb::log_lik()] matrix of each model and its relative efficiency, handed
#' to `loo::loo_model_weights()`, which stacks the leave-one-out
#' predictive densities (`method = "stacking"`, the default) or weights
#' the models by pseudo-BMA. [posterior_average()] takes its
#' `"stacking"` and `"pseudobma"` weights from the same function.
#'
#' @param x A `frmtmb_draws` from [frm_sample()].
#' @param ... Further draws objects, without names unless the name is
#'   the label to use; then `loo::loo_model_weights()`'s own options
#'   (`method`, `optim_method`, `optim_control`, `BB`, `BB_n`, `alpha`,
#'   `cores`), and `ndraws`, `draw_ids` and `resp` (either spelling),
#'   which reach every model's [frmtmb::log_lik()]. `newdata` and the
#'   new-level arguments are refused unless they hold brms's default.
#' @param model_names Labels for the models. The default deparses the
#'   arguments, as brms does.
#' @return loo's weights object, named by model.
#' @seealso [posterior_average()], [frmtmb::loo()]
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE) &&
#'     requireNamespace("loo", quietly = TRUE) &&
#'     !frmtmb.sample:::tmbstan_build_broken()) {
#'   set.seed(9)
#'   dd <- data.frame(x = rnorm(60))
#'   dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
#'   d1 <- frm_sample(bf(y ~ x), family = gaussian(), data = dd,
#'                    chains = 2, iter = 500, refresh = 0)
#'   d2 <- frm_sample(bf(y ~ 1), family = gaussian(), data = dd,
#'                    chains = 2, iter = 500, refresh = 0)
#'   loo_model_weights(d1, d2)
#' }
#' }
#' @name sample-loo_model_weights
NULL

#' @rdname sample-loo_model_weights
#' @exportS3Method loo::loo_model_weights
#' @export
loo_model_weights.frmtmb_draws <- function(x, ..., model_names = NULL) {
  dots <- list(...)
  nms <- names(dots) %||% rep("", length(dots))
  is_model <- vapply(dots, inherits, NA, "frmtmb_draws")
  bad <- which(!nzchar(nms) & !is_model)
  if (length(bad)) {
    frm_stop("loo_model_weights() takes models without a name, and ",
             "argument ", bad[1L] + 1L, " is a ",
             paste(class(dots[[bad[1L]]]), collapse = "/"), ", not draws. ",
             "Name the options, as method = \"pseudobma\"", call. = FALSE)
  }
  models <- c(list(x), dots[is_model])
  # brms's labels: the deparsed arguments, a name where one is given
  ex <- as.list(substitute(list(x, ...)))[-1L]
  labels <- vapply(ex[c(TRUE, is_model)],
                   function(e) paste(deparse(e), collapse = ""), "")
  given <- c("", nms[is_model])
  labels[nzchar(given)] <- given[nzchar(given)]
  if (!is.null(model_names)) {
    if (!is.character(model_names) ||
          length(model_names) != length(models)) {
      frm_stop("loo_model_weights(): `model_names` must name each of the ",
               length(models), " models", call. = FALSE)
    }
    labels <- model_names
  }
  opts <- dots[!is_model]
  ll_args <- intersect(names(opts), c("ndraws", "draw_ids", "resp"))
  out <- draws_model_weights(models, opts[setdiff(names(opts), ll_args)],
                             opts[ll_args], "loo_model_weights()")
  names(out) <- labels
  out
}

#' brms's `loo_model_weights.brmsfit()` on a list of draws: each model's
#' log_lik() and relative efficiency, then `loo::loo_model_weights()` on
#' the matrices. posterior_average() reaches it too, as brms's
#' `model_weights()` reaches brms's. `ll_opts` holds `ndraws`,
#' `draw_ids` and `resp`, which reach every model's log_lik().
#'
#' @noRd
draws_model_weights <- function(models, opts, ll_opts, what) {
  if (!requireNamespace("loo", quietly = TRUE)) {
    loo_need_pkg(what)
  }
  # as in the other loo_*() functions: brms's default values ask for
  # nothing log_lik() refuses, and pass
  ll_refused <- c("newdata", "re_formula", "re.form", draws_new_level_args,
                  "nug")
  hit <- Filter(function(a) !loo_default_value(a, opts[[a]]),
                intersect(names(opts), ll_refused))
  if (length(hit)) {
    frm_stop(what, " cannot honor `", hit[1L], "`: the weights are ",
             "computed from log_lik() of each model, which takes the rows ",
             "the model was fitted on and refuses it", call. = FALSE)
  }
  if ("r_eff_list" %in% names(opts)) {
    frm_stop(what, " cannot honor `r_eff_list`: the relative efficiencies ",
             "are each model's own, taken on its chains as loo() takes ",
             "them", call. = FALSE)
  }
  opts <- opts[setdiff(names(opts), ll_refused)]
  known <- setdiff(names(formals(loo::loo_model_weights.default)),
                   c("x", "...", "r_eff_list"))
  unknown <- setdiff(names(opts), known)
  if (length(unknown)) {
    frm_stop(what, " has no argument `", unknown[1L], "`. It takes ",
             paste(c(known, "ndraws", "draw_ids", "resp"), collapse = ", "),
             call. = FALSE)
  }
  ll <- lapply(models, function(m) {
    r <- loo_resp_name(draws_base_fit(m), ll_opts$resp, what)
    draws_as_caller(what, re_formula = FALSE,
                    loo_matrix(m, ll_opts$ndraws, r, what,
                               draw_ids = ll_opts$draw_ids))
  })
  n <- vapply(ll, ncol, 1L)
  if (length(unique(n)) != 1L) {
    frm_stop(what, " compares the models row by row, so they must be ",
             "fitted on the same rows; they have ",
             paste(n, collapse = ", "), " observations", call. = FALSE)
  }
  r_eff <- Map(loo_r_eff, ll, models)
  do.call(loo::loo_model_weights,
          c(list(unname(ll), r_eff_list = unname(r_eff)), opts))
}
