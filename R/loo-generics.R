# brms's leave-one-out functions that need posterior draws: the
# generics, and their refusals on a maximum-likelihood fit. The draws
# methods are frmtmb.sample's, registered on these generics.
#
# Why the generics live here and not in frmtmb.sample. A user who has
# fitted with frm() and calls brms's loo_R2() must be told to sample.
# With the generics in frmtmb.sample that user met "could not find
# function" with frmtmb alone attached, and "no applicable method for
# 'loo_R2' applied to an object of class "frmtmb_fit"" with rstantools
# or brms attached: neither names the reason or the route. log_lik(),
# loo() and bayes_R2() are here for the same reason. kfold(), reloo(),
# loo_moment_match() and loo_subsample() had the same gap (with loo
# attached, loo_moment_match(fit) died inside loo's default method on a
# missing `loo`), and frmtmb.sample's draws methods refuse them too, so
# their refusal on a fit names what to use instead of sampling.

#' Leave-one-out functions that need posterior draws
#'
#' brms's `loo_epred()`, `loo_linpred()`, `loo_predict()`,
#' `loo_predictive_interval()`, `loo_R2()`, `loo_model_weights()`,
#' `kfold()`, `reloo()`, `loo_moment_match()` and `loo_subsample()`. Each
#' weights or refits posterior draws, and
#' a `frmtmb_fit` has one parameter vector, so on a fit each refuses and
#' says what to do instead. The first six are computed on draws by the
#' `frmtmb.sample` package, which registers its methods on these
#' generics: sample with `frmtmb.sample::frm_sample()` first.
#' The last four refit the sampler or approximate `loo()` for data too
#' large to evaluate, which neither package
#' does; [frm_bootstrap()] resamples a fit, and [AIC()] compares fits.
#'
#' @param object,x A `frmtmb_fit`, or (with `frmtmb.sample` loaded)
#'   draws.
#' @param ... Passed to methods.
#' @return These methods signal an error on a maximum-likelihood fit.
#' @seealso [loo()], [log_lik()], [bayes_R2()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(40))
#' dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
#' fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
#' try(loo_predict(fit))
#' try(kfold(fit))
#' try(loo_subsample(fit))
#' @name loo_predict
NULL

#' @rdname loo_predict
#' @export
loo_epred <- function(object, ...) UseMethod("loo_epred")

#' @rdname loo_predict
#' @export
loo_linpred <- function(object, ...) UseMethod("loo_linpred")

#' @rdname loo_predict
#' @export
loo_predict <- function(object, ...) UseMethod("loo_predict")

#' @rdname loo_predict
#' @export
loo_predictive_interval <- function(object, ...) {
  UseMethod("loo_predictive_interval")
}

#' @rdname loo_predict
#' @export
loo_R2 <- function(object, ...) UseMethod("loo_R2")

#' @rdname loo_predict
#' @export
loo_model_weights <- function(x, ...) UseMethod("loo_model_weights")

#' @rdname loo_predict
#' @export
kfold <- function(x, ...) UseMethod("kfold")

#' @rdname loo_predict
#' @export
reloo <- function(x, ...) UseMethod("reloo")

#' @rdname loo_predict
#' @export
loo_moment_match <- function(x, ...) UseMethod("loo_moment_match")

#' @rdname loo_predict
#' @export
loo_subsample <- function(x, ...) UseMethod("loo_subsample")

#' The refusal of a draws-weighting function on a fit, naming the route
#' to draws.
#'
#' @noRd
loo_fit_refusal <- function(what, example = paste0(
                              sub("()", "", what, fixed = TRUE),
                              "(frmtmb.sample::frm_sample(fit))")) {
  frm_stop(what, " weights posterior draws by their leave-one-out ",
           "importance ratios, and this is a maximum-likelihood fit with ",
           "one parameter vector. Sample first: ", example,
           call. = FALSE)
}

#' The refusal of a refitting function on a fit. Sampling would not
#' help: frmtmb.sample refuses it on draws too. `loo_moment_match()`
#' is one: it moment-matches the draws toward each fold and refits the
#' folds that stay bad.
#'
#' @noRd
loo_refit_refusal <- function(what) {
  frm_stop(what, " refits the sampler on modified data, and neither ",
           "frmtmb nor frmtmb.sample does that, on a fit or on draws. ",
           "frm_bootstrap() resamples a maximum-likelihood fit, and AIC() ",
           "compares fits", call. = FALSE)
}

#' The refusal of `loo_subsample()` on a fit: an approximation of
#' `loo()` for data too large to evaluate, which frmtmb.sample does not
#' build either, because its `log_lik()` holds the whole matrix.
#'
#' @noRd
loo_unbuilt_refusal <- function(what) {
  frm_stop(what, " approximates loo() on a subsample of the ",
           "observations, and neither frmtmb nor frmtmb.sample builds it, ",
           "on a fit or on draws. loo(frmtmb.sample::frm_sample(fit)) ",
           "computes the full estimate; AIC() compares fits",
           call. = FALSE)
}

#' @rdname loo_predict
#' @exportS3Method brms::loo_epred
#' @rawNamespace S3method(rstantools::loo_epred,frmtmb_fit)
#' @export
loo_epred.frmtmb_fit <- function(object, ...) loo_fit_refusal("loo_epred()")

#' @rdname loo_predict
#' @exportS3Method rstantools::loo_linpred
#' @export
loo_linpred.frmtmb_fit <- function(object, ...) {
  loo_fit_refusal("loo_linpred()")
}

#' @rdname loo_predict
#' @exportS3Method rstantools::loo_predict
#' @export
loo_predict.frmtmb_fit <- function(object, ...) {
  loo_fit_refusal("loo_predict()")
}

#' @rdname loo_predict
#' @exportS3Method rstantools::loo_predictive_interval
#' @export
loo_predictive_interval.frmtmb_fit <- function(object, ...) {
  loo_fit_refusal("loo_predictive_interval()")
}

#' @rdname loo_predict
#' @exportS3Method rstantools::loo_R2
#' @export
loo_R2.frmtmb_fit <- function(object, ...) loo_fit_refusal("loo_R2()")

#' @rdname loo_predict
#' @exportS3Method loo::loo_model_weights
#' @export
loo_model_weights.frmtmb_fit <- function(x, ...) {
  loo_fit_refusal("loo_model_weights()", paste0(
    "loo_model_weights(frmtmb.sample::frm_sample(fit1), ",
    "frmtmb.sample::frm_sample(fit2))"))
}

#' @rdname loo_predict
#' @exportS3Method loo::kfold
#' @export
kfold.frmtmb_fit <- function(x, ...) loo_refit_refusal("kfold()")

#' @rdname loo_predict
#' @exportS3Method brms::reloo
#' @export
reloo.frmtmb_fit <- function(x, ...) loo_refit_refusal("reloo()")

#' @rdname loo_predict
#' @exportS3Method loo::loo_moment_match
#' @export
loo_moment_match.frmtmb_fit <- function(x, ...) {
  loo_refit_refusal("loo_moment_match()")
}

#' @rdname loo_predict
#' @exportS3Method loo::loo_subsample
#' @export
loo_subsample.frmtmb_fit <- function(x, ...) {
  loo_unbuilt_refusal("loo_subsample()")
}
