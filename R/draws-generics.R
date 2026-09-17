# The draws-shaped generics core defines and keeps.
#
# The sampling surface itself lives in the frmtmb.sample package
# (dev/draws-extraction.md). What stays here is the generics core still
# has a method for: `frm_multiple()` answers the whole as_draws /
# draw-count family with a refusal, and `posterior_summary()` has a
# default method that reduces any matrix of draws. frmtmb.sample
# registers its `frmtmb_draws` methods on these, so a session with both
# packages sees one generic per name rather than a masked pair.
#
# A generic with no core method at all left with its methods; see the
# "generics" section of dev/draws-extraction.md for the rule and the
# list.

#' Summaries and intervals of draws
#'
#' `posterior_summary()` reduces draws to estimate, error and quantiles
#' in brms's column layout (`Estimate`, `Est.Error`, `Q2.5`, `Q97.5`).
#' A matrix has variables in columns and draws in rows and gives one row
#' per variable. A three-dimensional array has draws in its first
#' margin, as `ranef(summary = FALSE)` returns them, and gives a
#' `variables x statistics x third-margin` array, as in brms.
#'
#' The computation is brms's `posterior_summary.default()`, statistic
#' for statistic: `mean` and `sd`, or `median` and `mad` when `robust`,
#' then `quantile`, each with `na.rm = TRUE`.
#'
#' The method for posterior draws of a fitted model is in the
#' `frmtmb.sample` package, along with the sampler that produces them.
#' A maximum-likelihood fit has no draws, and its method says so.
#'
#' @param x A matrix or three-dimensional array of draws.
#' @param probs Quantiles to report.
#' @param robust If `TRUE`, median and MAD instead of mean and SD.
#' @param ... Refused: an argument the method does not have is an
#'   error naming it, rather than silently changing nothing.
#' @return A matrix with one row per variable, or an array, as above.
#' @examples
#' # any matrix of draws: rows are draws, columns are variables
#' m <- cbind(a = rnorm(500), b = rnorm(500, 2))
#' posterior_summary(m)
#' posterior_summary(m, robust = TRUE)
#'
#' # a maximum-likelihood fit has no draws to summarize
#' dd <- data.frame(y = rnorm(40), x = rnorm(40))
#' try(posterior_summary(frm(bf(y ~ x) + gaussian(), data = dd)))
#' @export
posterior_summary <- function(x, ...) UseMethod("posterior_summary")

#' @rdname posterior_summary
#' @export
posterior_summary.default <- function(x, probs = c(0.025, 0.975),
                                      robust = FALSE, ...) {
  frm_check_dots(...)
  if (!length(x)) {
    frm_stop("No posterior draws supplied.", call. = FALSE)
  }
  coefs <- if (robust) c("median", "mad", "quantile") else
    c("mean", "sd", "quantile")
  one <- function(d) {
    do.call(cbind, lapply(coefs, draws_estimate, draws = d,
                          probs = probs, na.rm = TRUE))
  }
  x <- if (length(dim(x)) <= 2L) as.matrix(x) else as.array(x)
  if (length(dim(x)) == 2L) {
    out <- one(x)
    rownames(out) <- colnames(x)
  } else if (length(dim(x)) == 3L) {
    dnx <- dimnames(x)
    per <- lapply(seq_len(dim(x)[3L]), function(k) {
      one(array(x[, , k], dim = dim(x)[-3L]))
    })
    out <- array(unlist(per),
                 dim = c(nrow(per[[1L]]), ncol(per[[1L]]), length(per)))
    dimnames(out) <- list(dnx[[2L]], NULL, dnx[[3L]])
  } else {
    frm_stop("'x' must be of dimension 2 or 3.", call. = FALSE)
  }
  colnames(out) <- c("Estimate", "Est.Error", paste0("Q", probs * 100))
  out
}

#' One statistic over the columns of a draws matrix, the way
#' `brms:::get_estimate()` applies it: extra arguments reach the
#' statistic only if it can take them, and a statistic returning more
#' than one number per column (a quantile) is laid out one row per
#' column.
#'
#' @noRd
draws_estimate <- function(coef, draws, ...) {
  dots <- list(...)
  fun <- match.fun(coef)
  if (!"..." %in% names(formals(fun))) {
    dots <- dots[names(dots) %in% names(formals(fun))]
  }
  x <- do.call(apply, c(list(X = draws, MARGIN = 2L, FUN = fun), dots))
  if (is.null(dim(x))) {
    x <- matrix(x, dimnames = list(NULL, coef))
  } else if (coef == "quantile") {
    x <- aperm(x, length(dim(x)):1)
  }
  x
}

# Without these two methods a fit reached posterior_summary.default,
# which died inside as.matrix() with "is.atomic(x) is not TRUE". The
# brms registration matters as much as the local one: with brms loaded,
# the exported generic is brms's (R/generic-owners.R).
#' @rdname posterior_summary
#' @exportS3Method brms::posterior_summary
#' @export
posterior_summary.frmtmb_fit <- function(x, ...) {
  fit_no_draws("posterior_summary")
}

#' @rdname posterior_summary
#' @exportS3Method brms::posterior_summary
#' @export
posterior_summary.frmtmb_multiple <- function(x, ...) {
  multiple_no_draws("posterior_summary")
}

#' Convert to a posterior draws object
#'
#' The `as_draws` family converts an object holding posterior draws into
#' one of the posterior package's draws formats. frmtmb defines the
#' generics so that they work whether or not posterior is attached, and
#' registers methods with posterior so that its own spellings dispatch
#' too.
#'
#' Core has no object that carries draws: `frm()` is maximum likelihood
#' and `frm_multiple()` pools point estimates, so its methods explain
#' that rather than inventing a draws matrix. Install `frmtmb.sample`
#' and sample with `frmtmb.sample::frm_sample()` to get an object these
#' convert.
#'
#' @param x An object holding draws.
#' @param ... Passed to methods.
#' @return A posterior draws object of the requested format.
#' @examples
#' # a maximum-likelihood fit carries no draws, and says so by name
#' dd <- data.frame(y = rnorm(40), x = rnorm(40))
#' try(as_draws_df(frm(bf(y ~ x) + gaussian(), data = dd)))
#'
#' # frm_multiple() pools estimates rather than carrying draws, so it
#' # answers with the reason rather than a matrix
#' fits <- frm_multiple(bf(y ~ x) + gaussian(), data = list(dd, dd))
#' try(as_draws(fits))
#' @name as_draws
NULL

#' @rdname as_draws
#' @export
as_draws <- function(x, ...) UseMethod("as_draws")

#' @rdname as_draws
#' @export
as_draws_matrix <- function(x, ...) UseMethod("as_draws_matrix")

#' @rdname as_draws
#' @export
as_draws_array <- function(x, ...) UseMethod("as_draws_array")

#' @rdname as_draws
#' @export
as_draws_df <- function(x, ...) UseMethod("as_draws_df")

#' @rdname as_draws
#' @export
as_draws_list <- function(x, ...) UseMethod("as_draws_list")

#' @rdname as_draws
#' @export
as_draws_rvars <- function(x, ...) UseMethod("as_draws_rvars")

#' The refusal a `frmtmb_fit` gets from every draws accessor.
#'
#' Core is maximum likelihood: there is one parameter vector and no
#' chains, so there is nothing for these to convert. The method has to
#' EXIST rather than be left to fall through, because the exported
#' generic is posterior's whenever posterior is loaded, and a
#' `frmtmb_fit` is a bare list that posterior's own default walks into
#' and fails inside. Measured before this method existed:
#' `as_draws(fit)` gave "All list elements must be lists themselves".
#'
#' @noRd
fit_no_draws <- function(fn) {
  frm_stop(fn, "() needs posterior draws and a frmtmb_fit has none: ",
           "frm() is maximum likelihood, so it carries one parameter ",
           "vector rather than chains. Install frmtmb.sample and sample ",
           "first, with ", fn, "(frmtmb.sample::frm_sample(fit)); for the ",
           "point estimates and their covariance use fixef(), vcov() or ",
           "confint() on the fit itself", call. = FALSE)
}

#' @rdname as_draws
#' @exportS3Method posterior::as_draws
#' @export
as_draws.frmtmb_fit <- function(x, ...) {
  fit_no_draws("as_draws")
}

#' @rdname as_draws
#' @exportS3Method posterior::as_draws_matrix
#' @export
as_draws_matrix.frmtmb_fit <- function(x, ...) {
  fit_no_draws("as_draws_matrix")
}

#' @rdname as_draws
#' @exportS3Method posterior::as_draws_array
#' @export
as_draws_array.frmtmb_fit <- function(x, ...) {
  fit_no_draws("as_draws_array")
}

#' @rdname as_draws
#' @exportS3Method posterior::as_draws_df
#' @export
as_draws_df.frmtmb_fit <- function(x, ...) {
  fit_no_draws("as_draws_df")
}

#' @rdname as_draws
#' @exportS3Method posterior::as_draws_list
#' @export
as_draws_list.frmtmb_fit <- function(x, ...) {
  fit_no_draws("as_draws_list")
}

#' @rdname as_draws
#' @exportS3Method posterior::as_draws_rvars
#' @export
as_draws_rvars.frmtmb_fit <- function(x, ...) {
  fit_no_draws("as_draws_rvars")
}

#' Size of a draws object
#'
#' `ndraws()` counts the post-warmup draws (all chains pooled),
#' `niterations()` the draws per chain, `nchains()` the chains and
#' `nvariables()` the sampled parameters. The names and meanings are
#' posterior's; frmtmb registers methods with posterior so that the
#' generics work whether or not that package is attached.
#'
#' As with the `as_draws` family, the objects that answer these with a
#' number come from `frmtmb.sample::frm_sample()`.
#'
#' @param x An object holding draws.
#' @param ... Unused. Carried because posterior's `nvariables()`
#'   generic has it and frmtmb hands these generics back to
#'   posterior when posterior is loaded.
#' @return A single integer.
#' @examples
#' dd <- data.frame(y = rnorm(40), x = rnorm(40))
#' fits <- frm_multiple(bf(y ~ x) + gaussian(), data = list(dd, dd))
#' try(ndraws(fits))
#' @name draws-dimensions
NULL

# posterior's nchains(), ndraws() and niterations() generics take x
# alone and nvariables() takes (x, ...); dev/generics-audit2.R measured
# it.  The signatures have to agree because frmtmb hands these generics
# back to posterior when posterior is loaded.
#' @rdname draws-dimensions
#' @export
ndraws <- function(x) UseMethod("ndraws")

#' @rdname draws-dimensions
#' @export
nchains <- function(x) UseMethod("nchains")

#' @rdname draws-dimensions
#' @export
niterations <- function(x) UseMethod("niterations")

#' @rdname draws-dimensions
#' @export
nvariables <- function(x, ...) UseMethod("nvariables")

#' @rdname draws-dimensions
#' @exportS3Method posterior::ndraws
#' @export
ndraws.frmtmb_fit <- function(x) fit_no_draws("ndraws")

#' @rdname draws-dimensions
#' @exportS3Method posterior::nchains
#' @export
nchains.frmtmb_fit <- function(x) fit_no_draws("nchains")

#' @rdname draws-dimensions
#' @exportS3Method posterior::niterations
#' @export
niterations.frmtmb_fit <- function(x) fit_no_draws("niterations")

#' @rdname draws-dimensions
#' @exportS3Method posterior::nvariables
#' @export
nvariables.frmtmb_fit <- function(x, ...) {
  fit_no_draws("nvariables")
}
