# The leave-one-out and information-criterion generics, and what a
# MAXIMUM-LIKELIHOOD fit answers them with.
#
# The estimators themselves are posterior quantities and live in the
# frmtmb.sample package with the sampler that produces the draws
# (dev/draws-extraction.md). What stays here is every generic core has a
# method for, so that frmtmb.sample registers its `frmtmb_draws` methods
# on ONE generic per name rather than masking a pair, and so that a
# script ported from brms that calls one of these on a fit gets the
# reason and the route rather than "could not find function".

#' Approximate leave-one-out cross-validation
#'
#' `loo()` and `waic()` estimate the expected log predictive density of
#' a model from posterior draws: `loo()` by Pareto-smoothed importance
#' sampling, `waic()` by the widely applicable information criterion.
#' `loo_compare()` ranks several of them. `LOO()` and `WAIC()` are
#' brms's deprecated capitalized spellings.
#'
#' All of these average a likelihood over posterior draws, so on a
#' `frmtmb_fit` - which is one maximum-likelihood parameter vector -
#' they refuse and name the two routes to an answer: [AIC()] and
#' [BIC()], which are the maximum-likelihood analogues already
#' available on the fit, or sampling the model first.
#'
#' @section Sampling:
#' The estimators for posterior draws are in the `frmtmb.sample`
#' package, which also provides `frm_sample()`:
#'
#' ```
#' # install.packages("remotes")
#' remotes::install_github("aforren1/frmtmb",
#'                         subdir = "extensions/frmtmb.sample")
#' library(frmtmb.sample)
#' loo(frm_sample(fit))
#' ```
#'
#' It registers `frmtmb_draws` methods on these generics, so the
#' spellings on this page are the ones that keep working once it is
#' loaded.
#'
#' @param x A `frmtmb_fit`, or (with `frmtmb.sample` loaded) draws.
#' @param ... Passed to methods.
#' @return These methods signal an error on a maximum-likelihood fit.
#' @seealso [AIC()] and [BIC()] for the maximum-likelihood comparison,
#'   [frm_bootstrap()] for a resampling one.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(40))
#' dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
#' fit <- frm(bf(y ~ x) + gaussian(), data = dd)
#'
#' # the maximum-likelihood comparison is available directly
#' AIC(fit)
#' # the predictive one needs draws, and says so
#' try(loo(fit))
#' @export
loo <- function(x, ...) UseMethod("loo")

#' @rdname loo
#' @exportS3Method loo::loo
#' @export
loo.frmtmb_fit <- function(x, ...) {
  stop("loo() is a posterior quantity and this is a maximum-likelihood ",
       "fit: an elpd averages the likelihood over draws. Sample first, ",
       "with frmtmb.sample::loo(frmtmb.sample::frm_sample(fit)) once ",
       "that package is installed, or compare maximum-likelihood fits ",
       "with AIC() or BIC()", call. = FALSE)
}

#' @rdname loo
#' @export
waic <- function(x, ...) UseMethod("waic")

#' @rdname loo
#' @exportS3Method loo::waic
#' @export
waic.frmtmb_fit <- function(x, ...) {
  stop("waic() averages the likelihood over posterior draws and a ",
       "maximum-likelihood fit has none. Sample first, with ",
       "frmtmb.sample::waic(frmtmb.sample::frm_sample(fit)) once that ",
       "package is installed; AIC() is the maximum-likelihood analogue ",
       "already on the fit", call. = FALSE)
}

#' @rdname loo
#' @export
loo_compare <- function(x, ...) UseMethod("loo_compare")

#' @rdname loo
#' @export
loo_compare.default <- function(x, ...) {
  # frmtmb's generic would otherwise mask loo's own function for anyone
  # who attaches both, and `loo_compare(loo(d1), loo(d2))`, the
  # spelling the sampling package's help page recommends, would stop at
  # "no applicable method", because loo does not export its default
  # method for the search path to find. Reaching for loo's METHOD rather
  # than calling loo::loo_compare() is deliberate: dispatch from inside
  # this namespace would find this function again and recurse forever.
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("loo_compare() on already-computed criteria is the loo ",
         "package's own function, and the package is not installed",
         call. = FALSE)
  }
  fn <- utils::getS3method("loo_compare", "default",
                           envir = asNamespace("loo"))
  fn(x, ...)
}

#' @rdname loo
#' @export
LOO <- function(x, ...) UseMethod("LOO")

#' @rdname loo
#' @exportS3Method brms::LOO
#' @export
LOO.frmtmb_fit <- function(x, ...) {
  stop("LOO() is the deprecated brms spelling, and on a ",
       "maximum-likelihood fit there are no draws to average anyway. ",
       "Compare fits with AIC(), or install frmtmb.sample and call ",
       "loo(frm_sample(fit)) - the lowercase spelling is the current ",
       "one there too", call. = FALSE)
}

#' @rdname loo
#' @export
WAIC <- function(x, ...) UseMethod("WAIC")

#' @rdname loo
#' @exportS3Method brms::WAIC
#' @export
WAIC.frmtmb_fit <- function(x, ...) {
  stop("WAIC() is the deprecated brms spelling, and it averages over ",
       "posterior draws a maximum-likelihood fit does not have. AIC() ",
       "is already on the fit; install frmtmb.sample for the sampled ",
       "version, waic(frm_sample(fit))", call. = FALSE)
}

#' Bayesian R-squared
#'
#' The proportion of the outcome's variance a model explains, computed
#' per posterior draw. A maximum-likelihood fit has one parameter
#' vector rather than a posterior, so this refuses on a `frmtmb_fit`
#' and names the route to draws; the estimator itself is in the
#' `frmtmb.sample` package, on the same generic.
#'
#' @param object A `frmtmb_fit`, or (with `frmtmb.sample` loaded) draws.
#' @param ... Passed to methods.
#' @return This method signals an error on a maximum-likelihood fit.
#' @seealso [loo()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(40))
#' dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
#' fit <- frm(bf(y ~ x) + gaussian(), data = dd)
#' try(bayes_R2(fit))
#' @export
bayes_R2 <- function(object, ...) UseMethod("bayes_R2")

#' @rdname bayes_R2
#' @exportS3Method rstantools::bayes_R2
#' @export
bayes_R2.frmtmb_fit <- function(object, ...) {
  stop("bayes_R2() is computed per posterior draw and this is a ",
       "maximum-likelihood fit. Install frmtmb.sample and sample ",
       "first: bayes_R2(frm_sample(fit))", call. = FALSE)
}
