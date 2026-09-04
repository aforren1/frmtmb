# Leave-one-out cross-validation and the pointwise log-likelihood behind
# it, for frm_sample() draws. The statistical content of the brmsfit
# post-processing surface: everything here is a function of the
# ndraws x nobs matrix log_lik() returns.

#' Chain membership of each row of the draws matrix.
#'
#' `frm_sample()` stacks the chains with `rbind` in chain order, so the
#' chain identifier `loo::relative_eff()` wants is a `rep(each =)` over
#' the per-chain draw count. It is only meaningful for the WHOLE matrix:
#' a thinned subsample no longer has the autocorrelation structure the
#' relative efficiency estimates.
#'
#' @noRd
draws_chain_id <- function(x) {
  nc <- x$stanfit@sim$chains %||% 1L
  n <- nrow(x$draws)
  if (nc <= 1L || n %% nc != 0L) return(rep(1L, n))
  rep(seq_len(nc), each = n %/% nc)
}

#' Refuse draws that do not carry the random effects.
#'
#' `frm_sample(laplace = TRUE)` marginalizes the inner parameters, so
#' the draws hold only the outer ones. Every quantity in this file is
#' CONDITIONAL on a draw's own `b` (brms's convention, and the only one
#' available when the sampler visits `b` itself), so there is nothing to
#' condition on.
#'
#' @noRd
draws_require_b <- function(x, what) {
  fit <- x$fit
  if (!length(fit$frame[["re_blocks"]])) return(invisible(NULL))
  if (any(startsWith(colnames(x$draws), "b["))) return(invisible(NULL))
  stop(what, " needs draws of the random effects, and these draws come ",
       "from frm_sample(laplace = TRUE), which integrates them out ",
       "instead of sampling them. The pointwise log-density is the one ",
       "CONDITIONAL on each draw's own group-level values, so there is ",
       "nothing left to condition on. Resample without laplace = TRUE",
       call. = FALSE)
}

#' Refuse the likelihoods that do not factor into one term per row.
#'
#' The columns of a `log_lik()` matrix are observations, and every
#' consumer of it (`loo()`, `waic()`, PSIS) leaves ONE of them out. A
#' likelihood whose smallest independent unit is a group (an R-side
#' autocorrelation block, a hidden-Markov sequence, a group-level latent
#' class) has no such column: dropping one row of a sequence is not a
#' model refit anyone asked for. brms has no families in this position,
#' so there is no convention to follow and the honest answer is to say
#' so rather than to hand back a matrix whose leave-one-out meaning is
#' undefined.
#'
#' @noRd
draws_loglik_factors <- function(fit, what) {
  # [[ ]] throughout: `$` partial-matches, and a frame field that gains
  # a longer sibling name would silently change which structure is read
  frame <- fit$frame
  for (r in names(fit$spec$responses)) {
    # a structure that supplies its own `loglik` is exactly the one
    # whose density replaces the rowwise product, and it names its own
    # smallest independent unit; a structure that declares capabilities
    # and nothing else leaves the columns of the matrix well defined
    # `structure` and its `unit` field are documented parts of the
    # frmtmb_structure() protocol, so they are read directly rather
    # than through core's one-line accessors
    st <- fit$spec$responses[[r]]$family[["structure"]]
    unit <- if (!is.null((frame[["autocor"]] %||% list())[[r]])) {
      "an R-side residual correlation (ar/ma/arma/cosy/unstr) block"
    } else if (!is.null(st[["loglik"]])) {
      st[["unit"]] %||% "a group the likelihood does not factor within"
    } else {
      next
    }
    stop(what, " needs a likelihood that factors into one term per ",
         "observation, and the response '", r, "' does not: its ",
         "smallest independent unit is ", unit, ", so a column of the ",
         "matrix would be a GROUP and leaving one out would drop a ",
         "whole sequence. brms has no family in this position, so ",
         "there is no leave-one-out convention to follow. Compare ",
         "these models with AIC() on the ML fits, or with ",
         "frm_bootstrap()", call. = FALSE)
  }
  if (length(frame[["mi_map"]] %||% list())) {
    stop(what, " is not defined for a model with in-model imputation ",
         "(mi() / me()): a row whose response or predictor is latent ",
         "contributes the density of a PARAMETER as well as of an ",
         "observation, so its column would not be the observation's ",
         "own likelihood and leaving it out would not leave out the ",
         "latent value. Fit the completed data with frm_multiple() and ",
         "compare with AIC(), or use frm_bootstrap()", call. = FALSE)
  }
  invisible(NULL)
}

#' The per-observation conditional log-density at ONE parameter vector.
#'
#' Runs the same composition the objective runs (`row_lpdf()` for the
#' family density with `cens()` and `trunc()` folded in, then the case
#' weights), on numeric dpar values instead of on the tape.
#'
#' @noRd
draws_row_loglik <- function(fit, resp) {
  frame <- fit$frame
  rspecs <- fit$spec$responses
  dpv <- with_cs_offsets(fit, NULL, eval_dpars(fit))
  extra <- fit_extras(fit)
  n <- frame[["n_obs"]]
  if (isTRUE(fit$spec$rescor)) {
    # the joint gaussian likelihood contributes ONE K-variate density
    # per row, so the columns still index observations; this is brms's
    # log_lik for a set_rescor(TRUE) model
    rs <- names(rspecs)
    K <- length(rs)
    Z <- vapply(rs, function(r) {
      as.numeric((frame[["y"]][[r]] - dpv[[r]]$mu) / dpv[[r]]$sigma)
    }, numeric(n))
    lsig <- rowSums(vapply(rs, function(r) {
      rep(log(as.numeric(dpv[[r]]$sigma)), length.out = n)
    }, numeric(n)))
    C <- us_chol_cor(fit$estimates[["thetar"]], K)
    return(as.numeric(RTMB::dmvnorm(Z, 0, C, log = TRUE)) - lsig)
  }
  use <- resp %||% names(rspecs)
  out <- numeric(n)
  for (r in use) {
    av <- frame[["aterm_values"]][[r]]
    fam <- rspecs[[r]]$family
    yv <- frame[["y"]][[r]]
    ll <- row_lpdf(fam, yv, yv, dpv[[r]], av, extra)
    out <- out + (av[["weights"]] %||% 1) * as.numeric(ll)
  }
  out
}

#' Pointwise log-likelihood of posterior draws
#'
#' The `ndraws x nobs` matrix of per-observation log-densities, each row
#' evaluated at one draw's own parameter vector. Every leave-one-out and
#' information-criterion quantity in this package is a function of this
#' matrix, and [frmtmb::loo()], [frmtmb::waic()] and `loo::psis()` take it directly.
#'
#' @section What is conditioned on:
#' The density is CONDITIONAL on the draw's own group-level values:
#' `frm_sample()` samples `b` alongside everything else, so each row of
#' the draws matrix is a complete parameter vector and no integration is
#' left to do. This is exactly brms's convention, where the Stan model
#' also samples the group-level parameters. The consequence is worth
#' stating plainly: `loo()` on such a matrix is leave-one-OBSERVATION-out
#' with the groups held fixed, not leave-one-group-out, and for a model
#' with few observations per group the two differ.
#'
#' Addition terms enter exactly as they enter the fitted objective,
#' because they enter through the same code: `cens()` replaces a row's
#' density by the matching CDF difference, `trunc()` divides by the
#' window mass, `weights()` multiplies the row's contribution, and
#' `trials()` is the family's own argument. brms composes them in the
#' same order (`log_lik_censor()`, `log_lik_truncate()`,
#' `log_lik_weight()`).
#'
#' @section Multivariate models:
#' With `set_rescor(TRUE)` a column is the joint density of the row's
#' response VECTOR, so the matrix keeps one column per observation.
#' Without `rescor`, the responses are independent given the predictors
#' and the default sums their log-densities per row; pass `resp` to get
#' one response's contribution alone.
#'
#' @section Likelihoods with no per-observation column:
#' A model whose smallest independent unit is a group has no
#' per-observation column to leave out, and this refuses rather than
#' inventing one: R-side residual correlation ([frmtmb::frmtmb-autocor]), a
#' `frmtmb.latent::hmm()` sequence, and a group-level mixture (`mixture(groups = )`).
#' An `frmtmb.latent::lca()` subject is one row, so its column is well defined and is
#' not refused. In-model imputation (`mi()`, `me()`) is refused for the
#' same kind of reason: a latent value is a parameter, not an
#' observation. Use `AIC()` on the maximum-likelihood fits or
#' [frmtmb::frm_bootstrap()] for those.
#'
#' @param object A `frmtmb_draws` from [frm_sample()].
#' @param ndraws Number of draws to use, evenly spaced through the
#'   matrix (default: all of them).
#' @param resp For a multivariate model without `rescor`, the response
#'   whose contribution to report; the default sums over responses.
#' @param ... Unused.
#' @return A numeric matrix with one row per draw and one column per
#'   observation (the rows the model was fitted on).
#' @seealso [frmtmb::loo()], [frmtmb::waic()], [frmtmb::bayes_R2()]
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE)) {
#'   set.seed(9)
#'   dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
#'   dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
#'   fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
#'   ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0)
#'
#'   ll <- log_lik(ds)
#'   dim(ll)
#'   # the column means are the per-observation expected log-densities
#'   head(colMeans(ll))
#' }
#' }
#' @export
log_lik <- function(object, ...) {
  # own generic, for the reason pp_check() and as_draws() have one: the
  # packages that define this name (rstantools, brms) stay out of the
  # dependency list, and log_lik(ds) has to work without them
  UseMethod("log_lik")
}

#' @rdname log_lik
#' @exportS3Method rstantools::log_lik
#' @export
log_lik.frmtmb_draws <- function(object, ndraws = NULL, resp = NULL,
                                 ...) {
  fit <- draws_base_fit(object)
  draws_require_b(object, "log_lik()")
  draws_loglik_factors(fit, "log_lik()")
  if (!is.null(resp) && !resp %in% names(fit$spec$responses)) {
    stop("log_lik(resp = \"", resp, "\") names no response of this ",
         "model; it has ",
         paste(names(fit$spec$responses), collapse = ", "),
         call. = FALSE)
  }
  idx <- draws_par_index(fit)
  rows <- draws_subsample(object, ndraws)
  out <- NULL
  for (k in seq_along(rows)) {
    v <- draws_row_loglik(draws_fit_at(object, rows[k], idx), resp)
    if (is.null(out)) out <- matrix(NA_real_, length(rows), length(v))
    out[k, ] <- v
  }
  attr(out, "chain_id") <- if (length(rows) == nrow(object$draws)) {
    draws_chain_id(object)
  }
  out
}

#' Approximate leave-one-out cross-validation
#'
#' `loo()` runs Pareto-smoothed importance-sampling LOO and `waic()` the
#' widely applicable information criterion, both on the [log_lik()]
#' matrix, by handing it to `loo::loo.matrix()` and `loo::waic.matrix()`
#' unchanged. The returned objects are the loo package's own, so
#' `print()` and `loo::pareto_k_table()` work on them directly.
#' `loo_compare()` computes the criterion for each draws object it is
#' given and ranks them; handed criteria instead of draws, it is
#' `loo::loo_compare()` itself. `psis()` returns the smoothed importance
#' weights alone.
#'
#' `LOO()` and `WAIC()` are brms's deprecated capitalized spellings and
#' are defined only to name their replacements.
#'
#' @section Priors, and what these numbers mean:
#' These are posterior quantities, and they inherit the standing of the
#' draws they are computed from. Both of [frm_sample()]'s routes apply
#' brms's default priors, so these numbers are regularized the way brms
#' regularizes them unless the call opted out with `prior = "flat"`.
#' Under that opt-out the elpd is likelihood-shaped and unregularized:
#' expect Pareto k warnings for models with many group-level
#' parameters, because a flat prior leaves those to be identified by
#' the data alone,
#' and an influential observation then moves them a long way. The
#' maximum-likelihood answer to the same question is `AIC()` on the
#' fits, or [frmtmb::frm_bootstrap()].
#'
#' @section Relative efficiency:
#' `r_eff` defaults to `loo::relative_eff()` on the chain structure of
#' the draws, which is what brms does. Thinning with `ndraws` breaks
#' that structure, so `r_eff` is then dropped and the estimate is the
#' one loo computes without an autocorrelation correction.
#'
#' @param x A `frmtmb_draws` from [frm_sample()], or (for
#'   `loo_compare()`) already-computed criteria.
#' @param ndraws,resp Passed to [log_lik()].
#' @param ... Further models for `loo_compare()`; otherwise passed to
#'   the loo package function.
#' @return A `loo`, `waic`, `compare.loo` or `psis` object from the loo
#'   package.
#' @seealso [log_lik()], [frmtmb::bayes_R2()]
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE) &&
#'     requireNamespace("loo", quietly = TRUE)) {
#'   set.seed(9)
#'   dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
#'   dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
#'
#'   # sample with priors: an elpd is a posterior quantity
#'   d1 <- frm_sample(bf(y ~ x + (1 | g)), family = gaussian(),
#'                    data = dd, chains = 1, iter = 500, refresh = 0)
#'   d2 <- frm_sample(bf(y ~ 1 + (1 | g)), family = gaussian(),
#'                    data = dd, chains = 1, iter = 500, refresh = 0)
#'   loo_compare(d1, d2)
#'   # the same thing, one step at a time
#'   loo_compare(loo(d1), loo(d2))
#' }
#' }
#' @name sample-loo
NULL

#' @rdname sample-loo
#' @exportS3Method loo::loo
#' @export
loo.frmtmb_draws <- function(x, ndraws = NULL, resp = NULL, ...) {
  # a second model arrives positionally as ndraws, so check the named
  # formals too
  loo_one_model(c(list(ndraws, resp), list(...)), "loo")
  ll <- loo_matrix(x, ndraws, resp, "loo()")
  loo::loo.matrix(ll, r_eff = loo_r_eff(ll), ...)
}


#' @rdname sample-loo
#' @exportS3Method loo::waic
#' @export
waic.frmtmb_draws <- function(x, ndraws = NULL, resp = NULL, ...) {
  loo_one_model(c(list(ndraws, resp), list(...)), "waic")
  # waic needs no importance weights, so it needs no r_eff either
  loo::waic.matrix(loo_matrix(x, ndraws, resp, "waic()"), ...)
}


#' brms's loo(a, b) compares in one call; here comparison is its own
#' verb, and loo::loo.matrix would otherwise die coercing the second
#' model to an integer.
#'
#' @noRd
loo_one_model <- function(dots, what) {
  extra <- vapply(dots, inherits, logical(1L),
                  what = c("frmtmb_draws", "frmtmb_fit"))
  if (any(extra)) {
    stop(what, "() takes one model here, not the several brms ",
         "compares in a single call. Pass them all to loo_compare(), ",
         "which computes one elpd per model and tables the differences",
         call. = FALSE)
  }
}

#' The log-likelihood matrix with the loo package present.
#'
#' @noRd
loo_matrix <- function(x, ndraws, resp, what) {
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop(what, " needs the 'loo' package; install it, or call ",
         "log_lik() and pass the matrix to your own estimator",
         call. = FALSE)
  }
  log_lik(x, ndraws = ndraws, resp = resp)
}

#' Relative efficiency of the likelihood ratios, or NULL when the draws
#' were thinned and the chain structure no longer describes them.
#'
#' @noRd
loo_r_eff <- function(ll) {
  cid <- attr(ll, "chain_id")
  if (is.null(cid)) return(NULL)
  loo::relative_eff(exp(ll), chain_id = cid)
}


#' @rdname sample-loo
#' @param criterion Which criterion `loo_compare()` computes for each
#'   draws object.
#' @param model_names Row names for the comparison; the default deparses
#'   the arguments, as loo does.
#' @exportS3Method loo::loo_compare
#' @export
loo_compare.frmtmb_draws <- function(x, ..., criterion = c("loo", "waic"),
                                     model_names = NULL) {
  criterion <- match.arg(criterion)
  models <- c(list(x), list(...))
  bad <- which(!vapply(models, inherits, TRUE, "frmtmb_draws"))
  if (length(bad)) {
    stop("loo_compare() on draws computes the criterion for every ",
         "model, so every argument has to be a frmtmb_draws object; ",
         "argument ", bad[1L], " is a ",
         paste(class(models[[bad[1L]]]), collapse = "/"),
         ". Compute loo() on each model first and compare those",
         call. = FALSE)
  }
  crit <- lapply(models, function(m) {
    if (identical(criterion, "loo")) loo(m) else waic(m)
  })
  names(crit) <- model_names %||%
    loo_call_names(match.call(), length(models))
  loo::loo_compare(crit)
}

#' Row names for a comparison, deparsed from the model arguments of the
#' call (loo's own convention for the models it is handed). `UseMethod`
#' names the first argument `x` in the matched call, so the options are
#' excluded by name rather than the models being selected by position.
#'
#' @noRd
loo_call_names <- function(cl, n) {
  a <- as.list(cl)[-1L]
  nm <- names(a) %||% rep("", length(a))
  a <- a[!nm %in% c("criterion", "model_names")]
  out <- vapply(a, function(e) paste(deparse(e), collapse = ""), "")
  if (length(out) != n) paste0("model", seq_len(n)) else out
}

#' @rdname sample-loo
#' @param log_ratios For `psis()`, the draws object whose negative
#'   pointwise log-likelihood supplies the importance ratios.
#' @export
psis <- function(log_ratios, ...) UseMethod("psis")

#' @rdname sample-loo
#' @exportS3Method loo::psis
#' @export
psis.frmtmb_draws <- function(log_ratios, ndraws = NULL, resp = NULL,
                              ...) {
  # the leave-one-out importance ratios are the NEGATIVE log-likelihood
  # (down-weighting the draws the held-out point likes), which is what
  # loo::loo.matrix smooths internally and what brms's psis.brmsfit
  # passes on
  ll <- loo_matrix(log_ratios, ndraws, resp, "psis()")
  loo::psis(-ll, r_eff = loo_r_eff(ll), ...)
}

#' @rdname sample-loo
#' @exportS3Method brms::LOO
#' @export
LOO.frmtmb_draws <- function(x, ...) {
  stop("LOO() is the deprecated brms spelling and frmtmb never had it. ",
       "Use loo(x), whose result is a loo-package object that ",
       "loo::loo_compare() and loo::pareto_k_table() read directly",
       call. = FALSE)
}


#' @rdname sample-loo
#' @exportS3Method brms::WAIC
#' @export
WAIC.frmtmb_draws <- function(x, ...) {
  stop("WAIC() is the deprecated brms spelling and frmtmb never had ",
       "it. Use waic(x); note that loo(x) is the better-behaved ",
       "estimator of the same predictive quantity", call. = FALSE)
}


#' Bayesian R-squared
#'
#' The proportion of the outcome's variance the model explains, computed
#' per draw and returned as a posterior distribution. This is the
#' residual-based estimator of Gelman, Goodrich, Gabry and Vehtari
#' (2019), *R-squared for Bayesian regression models*, The American
#' Statistician 73(3), which brms implements as
#' `var(ypred_s) / (var(ypred_s) + var(y - ypred_s))` for each draw `s`,
#' with both variances taken over the observations. frmtmb computes the
#' same expression on [posterior_epred()] draws, so the two agree up to
#' Monte Carlo error on the same model.
#'
#' Because the denominator is a per-draw variance rather than a fixed
#' total sum of squares, the value is bounded in `(0, 1)` by
#' construction and does not have the classical R-squared's habit of
#' exceeding 1 on posterior draws.
#'
#' @param object A `frmtmb_draws` from [frm_sample()].
#' @param resp For a multivariate model, which response.
#' @param summary If `TRUE` (the default, as in brms), summarize the
#'   draws into estimate, error and quantiles; if `FALSE`, return the
#'   `ndraws x 1` matrix of R-squared draws.
#' @param probs Quantiles for the summary.
#' @param ndraws Number of draws to use (default: all).
#' @param ... Unused.
#' @return A one-row summary matrix, or the matrix of draws when
#'   `summary = FALSE`.
#' @seealso [log_lik()], [frmtmb::loo()]
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE)) {
#'   set.seed(9)
#'   dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
#'   dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
#'   ds <- frm_sample(bf(y ~ x + (1 | g)), family = gaussian(),
#'                    data = dd, chains = 1, iter = 500, refresh = 0)
#'   bayes_R2(ds)
#'   quantile(bayes_R2(ds, summary = FALSE), c(0.1, 0.9))
#' }
#' }
#' @name sample-bayes_R2
NULL


#' @rdname sample-bayes_R2
#' @exportS3Method rstantools::bayes_R2
#' @export
bayes_R2.frmtmb_draws <- function(object, resp = NULL, summary = TRUE,
                                  probs = c(0.025, 0.975),
                                  ndraws = NULL, ...) {
  fit <- draws_base_fit(object)
  resp <- resp %||% names(fit$spec$responses)[1L]
  y <- fit$frame[["y"]][[resp]]
  if (is.null(y) || is.matrix(y)) {
    stop("bayes_R2() needs a single numeric response column, and '",
         resp %||% "?", "' is not one. A categorical, multinomial or ",
         "multivariate-item outcome has no residual variance to ",
         "decompose; brms refuses it for the same reason", call. = FALSE)
  }
  ep <- posterior_epred(object, resp = resp, ndraws = ndraws)
  if (length(dim(ep)) > 2L) {
    stop("bayes_R2() is not defined for an ordinal or categorical ",
         "family: posterior_epred() gives a category DISTRIBUTION per ",
         "observation, and treating those probabilities as a ",
         "continuous prediction (which is what brms does, with a ",
         "warning) makes the ratio uninterpretable. Use log_lik() and ",
         "loo() to compare such models", call. = FALSE)
  }
  # in-sample predictions are padded back to the original rows for an
  # na.exclude fit; the response is not
  ep <- ep[, !is.na(ep[1L, ]), drop = FALSE]
  y <- as.numeric(y)
  vp <- apply(ep, 1L, stats::var)
  ve <- apply(sweep(ep, 2L, y), 1L, stats::var)
  # brms names the column "R2" on a univariate model and "R2<resp>" on a
  # multivariate one, because its resp argument is NULL in the first case
  lab <- if (length(fit$spec$responses) > 1L) paste0("R2", resp) else "R2"
  R2 <- matrix(vp / (vp + ve), ncol = 1L, dimnames = list(NULL, lab))
  if (summary) posterior_summary(R2, probs = probs) else R2
}

# ---- refusals -------------------------------------------------------
#
# Every one of these is a brmsfit method a ported script may call. They
# fail with a reason and a pointer rather than "could not find function"
# (see the post-processing column of dev/brms-vignette-port.md), and
# each names itself so the message identifies one line of source.

#' Refusals for the refit-based and marginal-likelihood brmsfit methods
#'
#' These `brmsfit` methods have frmtmb spellings that do not exist yet.
#' They are defined so that a ported script fails with the reason and
#' the alternative rather than with "could not find function", and they
#' are documented so the reason is findable.
#'
#' * `loo_moment_match()`, `loo_subsample()`, `reloo()` and `kfold()`
#'   all need to refit the model on modified data. `frm_bootstrap()` is
#'   the resampling machinery frmtmb does have, and `AIC()` on the
#'   maximum-likelihood fits answers the comparison question directly.
#' * `bridge_sampler()`, `bayes_factor()` and `post_prob()` are
#'   marginal-likelihood quantities. A marginal likelihood is an
#'   integral against the PRIOR, so it is undefined under
#'   `prior = "flat"`, and even under the default priors
#'   the bridge-sampling estimator needs a normalized
#'   log-posterior evaluator that the RTMB tape does not expose.
#'
#' @param x,... Ignored; these methods always stop.
#' @return These functions never return; they signal an error.
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE)) {
#'   set.seed(1)
#'   dd <- data.frame(x = rnorm(40))
#'   dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
#'   fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
#'   ds <- frm_sample(fit, chains = 1, iter = 400, refresh = 0)
#'   # each refusal names its reason and the replacement
#'   try(reloo(ds))
#'   try(bayes_factor(ds, ds))
#' }
#' }
#' @name frmtmb-loo-refusals
NULL

#' @rdname frmtmb-loo-refusals
#' @export
loo_moment_match <- function(x, ...) UseMethod("loo_moment_match")

#' @rdname frmtmb-loo-refusals
#' @exportS3Method loo::loo_moment_match
#' @export
loo_moment_match.frmtmb_draws <- function(x, ...) {
  stop("loo_moment_match() is not implemented for frmtmb draws: it ",
       "moment-matches the posterior toward each problem fold and then ",
       "falls back to refitting the ones that stay bad, and ",
       "frm_sample() has no stored program to refit. Read loo()'s ",
       "Pareto k table instead; if the flagged points are few, ",
       "frm_bootstrap() answers the influence question directly",
       call. = FALSE)
}

#' @rdname frmtmb-loo-refusals
#' @export
loo_subsample <- function(x, ...) UseMethod("loo_subsample")

#' @rdname frmtmb-loo-refusals
#' @exportS3Method loo::loo_subsample
#' @export
loo_subsample.frmtmb_draws <- function(x, ...) {
  stop("loo_subsample() is not implemented for frmtmb draws: its ",
       "point is a cheap approximation for models too large to hold a ",
       "full log_lik() matrix, and it needs a per-observation ",
       "likelihood callback plus the refit machinery to correct the ",
       "subsample. log_lik() builds the whole matrix here, so compute ",
       "loo() on it; thin the draws with loo(ndraws =) if memory is ",
       "the problem", call. = FALSE)
}

#' @rdname frmtmb-loo-refusals
#' @export
reloo <- function(x, ...) UseMethod("reloo")

#' @rdname frmtmb-loo-refusals
#' @exportS3Method brms::reloo
#' @export
reloo.frmtmb_draws <- function(x, ...) {
  stop("reloo() is not implemented for frmtmb draws: it re-runs the ",
       "sampler once per observation with a high Pareto k, and ",
       "frm_sample() has no stored program to re-run on modified ",
       "data. Read loo()'s Pareto k table and treat a bad k as the ",
       "diagnostic it is (usually many group-level parameters left to ",
       "the data alone; see the prior section of ?loo), or compare the ",
       "maximum-likelihood fits with AIC()", call. = FALSE)
}

#' @rdname frmtmb-loo-refusals
#' @export
kfold <- function(x, ...) UseMethod("kfold")

#' @rdname frmtmb-loo-refusals
#' @exportS3Method loo::kfold
#' @export
kfold.frmtmb_draws <- function(x, ...) {
  stop("kfold() is not implemented for frmtmb draws: K refits of the ",
       "sampler are exactly the refit machinery that is out of scope ",
       "here, and a partial version that silently used the ML fits ",
       "instead would not be the cross-validation the name promises. ",
       "Use loo() for the importance-sampling approximation, or ",
       "frm_bootstrap() for a resampling answer", call. = FALSE)
}

#' @rdname frmtmb-loo-refusals
#' @export
bridge_sampler <- function(x, ...) UseMethod("bridge_sampler")

#' @rdname frmtmb-loo-refusals
#' @exportS3Method bridgesampling::bridge_sampler
#' @export
bridge_sampler.frmtmb_draws <- function(x, ...) {
  stop("bridge_sampler() is not available for frmtmb draws. A marginal ",
       "likelihood is an integral of the likelihood against the PRIOR, ",
       "so it does not exist at all under prior = \"flat\"; and even ",
       "under the default priors ",
       "the estimator needs to evaluate the normalized log posterior ",
       "at arbitrary points, which the RTMB tape does not expose. Use ",
       "loo() for predictive comparison", call. = FALSE)
}

#' @rdname frmtmb-loo-refusals
#' @export
bayes_factor <- function(x, ...) UseMethod("bayes_factor")

#' @rdname frmtmb-loo-refusals
#' @exportS3Method bridgesampling::bayes_factor
#' @export
bayes_factor.frmtmb_draws <- function(x, ...) {
  stop("bayes_factor() is not available for frmtmb draws: it is a ",
       "ratio of the marginal likelihoods bridge_sampler() would have ",
       "to estimate, and those are undefined under prior = \"flat\" ",
       "and unavailable from the tape. hypothesis() gives the posterior ",
       "probability of a directional claim, and loo() the predictive ",
       "comparison", call. = FALSE)
}

#' @rdname frmtmb-loo-refusals
#' @export
post_prob <- function(x, ...) UseMethod("post_prob")

#' @rdname frmtmb-loo-refusals
#' @exportS3Method bridgesampling::post_prob
#' @export
post_prob.frmtmb_draws <- function(x, ...) {
  stop("post_prob() is not available for frmtmb draws: a posterior ",
       "model probability is normalized marginal likelihoods, which ",
       "bridge_sampler() would have to estimate and cannot here. ",
       "Compare models with loo() and loo::loo_compare()",
       call. = FALSE)
}
