# Small accessor methods: the conventional S3 surface that downstream
# packages (insight, performance, step, and friends) dispatch on.

#' Resolve one linear predictor by response / dpar name.
#'
#' @noRd
find_linpred <- function(object, resp = NULL, dpar = "mu") {
  hits <- Filter(function(lp) {
    lp$dpar == dpar && (is.null(resp) || lp$resp == resp)
  }, object$frame$linpreds)
  if (!length(hits)) {
    stop("No linear predictor for dpar '", dpar, "'",
         if (!is.null(resp)) paste0(" of response '", resp, "'") else "",
         call. = FALSE)
  }
  if (length(hits) > 1L) {
    stop("Multiple responses have dpar '", dpar,
         "'; disambiguate with resp = ", call. = FALSE)
  }
  hits[[1L]]
}

#' Residual standard deviation
#'
#' Returns the estimated `sigma` distributional parameter on the
#' response scale when it is constant across observations
#' (intercept-only or fixed). When `sigma` is modeled with covariates
#' the scalar summary does not exist; the method warns and returns `NA`
#' (use `predict(dpar = "sigma")` for the per-observation values).
#' Families without a `sigma` parameter return 1, following glmmTMB.
#'
#' @param object A `frmtmb_fit`.
#' @param ... Unused.
#' @return A scalar, or a named vector for multivariate fits.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#' # the residual SD, on the response scale
#' sigma(fit)
#' # which is what the standardized residuals divide by
#' max(abs(residuals(fit, type = "pearson") -
#'           residuals(fit) / sigma(fit)))
#'
#' # a poisson fit has no dispersion parameter, so sigma() is 1
#' dd$cnt <- rpois(100, exp(0.5 + 0.3 * dd$x))
#' sigma(frm(bf(cnt ~ x) + poisson(), data = dd))
#' @export
sigma.frmtmb_fit <- function(object, ...) {
  out <- vapply(object$spec$responses, function(rsp) {
    if (!"sigma" %in% rsp$family$dpars) return(1)
    av <- object$frame$aterm_values[[rsp$resp_name]]
    if (!is.null(av[["se"]]) && !isTRUE(av[["se_sigma"]])) {
      return(0)   # se() fixes the residual SD per observation
    }
    lp <- NULL
    for (l in object$frame$linpreds) {
      if (l$resp == rsp$resp_name && l$dpar == "sigma") lp <- l
    }
    if (is.null(lp)) return(1)
    if (!is.null(lp$constant)) return(lp$constant)
    if (ncol(lp$X) == 1L && identical(colnames(lp$X), "(Intercept)") &&
        is.null(lp$Z)) {
      return(lp$link$linkinv(object$estimates[[lp$par]][lp$idx]))
    }
    warning("sigma varies by observation; returning NA ",
            "(use predict(dpar = \"sigma\"))", call. = FALSE)
    NA_real_
  }, numeric(1))
  if (length(out) == 1L) unname(out) else out
}

#' @export
terms.frmtmb_fit <- function(x, resp = NULL, dpar = "mu", ...) {
  find_linpred(x, resp, dpar)$terms
}

#' @export
model.matrix.frmtmb_fit <- function(object, resp = NULL, dpar = "mu",
                                    ...) {
  find_linpred(object, resp, dpar)$X
}

#' @export
weights.frmtmb_fit <- function(object, resp = NULL, ...) {
  rn <- if (is.null(resp)) names(object$frame$y)[1L] else resp
  w <- object$frame$aterm_values[[rn]][["weights"]]
  if (is.null(w)) rep(1, object$frame$n_obs) else w
}

#' @export
na.action.frmtmb_fit <- function(object, ...) {
  object$frame$na_action
}

#' @export
deviance.frmtmb_fit <- function(object, ...) {
  -2 * as.numeric(logLik(object))
}

#' @export
extractAIC.frmtmb_fit <- function(fit, scale = 0, k = 2, ...) {
  ll <- logLik(fit)
  edf <- attr(ll, "df")
  c(edf, -2 * as.numeric(ll) + k * edf)
}

#' Number of levels per random-effect grouping factor
#'
#' @param object A `frmtmb_fit`.
#' @param ... Unused.
#' @return A named integer vector (smooth terms are excluded).
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100),
#'                  g = factor(rep(1:10, 10)),
#'                  h = factor(rep(1:4, each = 25)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g) + (1 | h)) + gaussian(), data = dd)
#'
#' # one count per distinct grouping factor
#' ngrps(fit)
#' # the count that decides whether a variance component is trustworthy,
#' # and the unit influence() deletes when given `groups`
#' ngrps(fit)[["h"]]
#' @export
ngrps <- function(object, ...) UseMethod("ngrps")

#' @rdname ngrps
#' @exportS3Method brms::ngrps
#' @rawNamespace S3method(lme4::ngrps,frmtmb_fit)
#' @export
ngrps.frmtmb_fit <- function(object, ...) {
  bks <- Filter(function(bk) bk$covstruct != "smooth",
                object$frame$re_blocks)
  ng <- vapply(bks, `[[`, 0L, "n_levels")
  names(ng) <- vapply(bks, `[[`, "", "group_name")
  ng[!duplicated(names(ng))]
}

#' Priors used in a fit
#'
#' On draws from `frmtmb.sample::frm_sample()` this reports the priors
#' the sampler applied, which on the formula interface includes the brms
#' default
#' priors it chose (see the Default priors section of that
#' function).
#'
#' @param object A `frmtmb_fit`, or a `frmtmb_draws` from
#'   `frmtmb.sample::frm_sample()`.
#' @param ... Unused.
#' @return The `frmtmb_priorlist` the fit was penalized with, or that
#'   the sampler used, or (invisibly) `NULL` when there were none.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#'
#' # what was actually applied, after set_prior() was matched to the
#' # coefficients of this design
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
#'            prior = set_prior("normal(0, 1)", class = "b") +
#'                     set_prior("exponential(1)", class = "sd"))
#' prior_summary(fit)
#'
#' # a plain maximum-likelihood fit reports that it had none
#' prior_summary(frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd))
#' @export
prior_summary <- function(object, ...) UseMethod("prior_summary")

#' @rdname prior_summary
#' @exportS3Method rstantools::prior_summary
#' @export
prior_summary.frmtmb_fit <- function(object, ...) {
  if (is.null(object$prior)) {
    cat("No priors were set (plain maximum likelihood).\n")
    return(invisible(NULL))
  }
  object$prior
}

#' Refit a model to a new response
#'
#' Reuses the assembled design (no formula parsing, no frame assembly)
#' and warm-starts the optimizer at the previous estimates, so a refit
#' costs one re-tape plus the optimization. This is the engine for
#' parametric bootstrap: simulate responses with [simulate()], refit to
#' each.
#'
#' @param object A `frmtmb_fit` for a univariate model.
#' @param newresp Replacement response: a vector of the original length,
#'   or a matrix of the original dimensions for matrix responses.
#' @param start Optional named start list (as in [frm()]); when given it
#'   replaces the warm start.
#' @param ... Unused.
#' @return A new `frmtmb_fit`.
#' @export
refit <- function(object, newresp, ...) UseMethod("refit")

#' @examples
#' set.seed(2)
#' dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
#' dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#' # refit to a simulated response (the parametric-bootstrap step)
#' ysim <- simulate(fit, nsim = 1, re.form = NA)[[1]]
#' rf <- refit(fit, ysim)
#' fixef(rf)
#' @rdname refit
#' @export
refit.frmtmb_fit <- function(object, newresp, start = NULL, ...) {
  frame <- object$frame
  if (length(frame$y) != 1L) {
    stop("refit() supports univariate models", call. = FALSE)
  }
  y0 <- frame$y[[1L]]
  if (is.matrix(y0)) {
    newresp <- as.matrix(newresp)
    if (!identical(dim(newresp), dim(y0))) {
      stop("newresp must be a ", nrow(y0), " x ", ncol(y0), " matrix",
           call. = FALSE)
    }
  } else {
    if (is.data.frame(newresp)) newresp <- newresp[[1L]]
    newresp <- as.vector(newresp)
    if (length(newresp) != length(y0)) {
      stop("newresp must have length ", length(y0), call. = FALSE)
    }
  }
  frame$y[[1L]] <- newresp
  fit_assembled(object$spec, frame, object$bform, object$call,
                REML = object$REML, start = start,
                control = object$control %||% frmtmb_control(),
                se = FALSE, lower = object$lower, upper = object$upper,
                prior = object$prior,
                quadrature = isTRUE(object$quadrature),
                template = if (is.null(start)) object$estimates,
                data2 = object$data2 %||% list())
}
