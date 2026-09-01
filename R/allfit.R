# Refit under every available optimizer (the lme4::allFit analog).

#' Derivative-free optimizers from Suggests, wrapped to the
#' `frmtmb_control()` custom-optimizer contract.
#'
#' @noRd
allfit_optimizers <- function() {
  opts <- list(nlminb = "nlminb", optim = "optim")
  if (requireNamespace("minqa", quietly = TRUE)) {
    opts$bobyqa <- function(par, fn, gr, lower, upper, control) {
      r <- minqa::bobyqa(par, fn, lower = lower, upper = upper)
      list(par = r$par, objective = r$fval,
           convergence = as.integer(r$ierr != 0), message = r$msg %||% "")
    }
  }
  if (requireNamespace("nloptr", quietly = TRUE)) {
    opts$nloptr_lbfgs <- function(par, fn, gr, lower, upper, control) {
      # nloptr inspects its callbacks' formals and refuses any that
      # declare `...` unless the extra arguments are handed to nloptr
      # itself; RTMB's obj$fn/obj$gr are function(x, ...). The
      # single-argument wrappers also drop obj$fn's "logarithm"
      # attribute and flatten obj$gr's 1-row matrix, both of which
      # nloptr wants as plain numerics.
      f <- function(x) as.numeric(fn(x))
      g <- function(x) as.vector(gr(x))
      # nloptr insists on bounds of length(par); the optimizer contract
      # allows a recycled scalar.
      lb <- rep_len(if (length(lower)) lower else -Inf, length(par))
      ub <- rep_len(if (length(upper)) upper else Inf, length(par))
      # ftol_rel is what actually stops L-BFGS here: on an x-tolerance
      # alone the line search collapses at the optimum before the step
      # ever gets small enough, and NLopt then reports NLOPT_FAILURE at
      # the right answer. Relative tolerances keep this scale-free.
      r <- nloptr::nloptr(par, f, g, lb = lb, ub = ub,
                          opts = list(algorithm = "NLOPT_LD_LBFGS",
                                      ftol_rel = 1e-10, xtol_rel = 1e-8,
                                      maxeval = 2000))
      list(par = stats::setNames(r$solution, names(par)),
           objective = r$objective,
           # NLopt reports success as a positive status, but 5 and 6
           # are the evaluation and time budgets running out, which is
           # not convergence
           convergence = as.integer(!r$status %in% 1:4),
           message = r$message %||% "")
    }
  }
  opts
}

#' Refit a model with every available optimizer
#'
#' The lme4 `allFit()` idea: reuse the assembled design and refit under
#' each optimizer, to separate optimizer trouble from model
#' misspecification. Uses nlminb and optim (L-BFGS-B) always, plus
#' bobyqa (minqa) and NLopt L-BFGS (nloptr) when those packages are
#' installed.
#'
#' @param fit A `frmtmb_fit`.
#' @param optimizers Named list of optimizers (names or functions, as
#'   in [frmtmb_control()]). Default: everything available.
#' @param ... Unused.
#' @return A `frmtmb_allfit` object: `$fits` (the refits, `NULL` where
#'   one errored) and a printed comparison of logLik, convergence, and
#'   fixed-effect spread.
#' @examples
#' set.seed(6)
#' dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
#' dd$y <- rpois(60, exp(0.3 + 0.4 * dd$x + rnorm(6, 0, 0.4)[dd$g]))
#' fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
#' frm_allfit(fit)
#' @export
frm_allfit <- function(fit, optimizers = NULL, ...) {
  optimizers <- optimizers %||% allfit_optimizers()
  ctl0 <- fit$control %||% frmtmb_control()
  # print.frmtmb_allfit() already reports the per-optimizer timings, so
  # a verbose original fit must not make every refit narrate itself
  ctl0$verbose <- FALSE
  fits <- list()
  times <- numeric(length(optimizers))
  for (i in seq_along(optimizers)) {
    ctl <- ctl0
    ctl$optimizer <- optimizers[[i]]
    if (!identical(ctl$optimizer, "nlminb") &&
        !identical(ctl$optimizer, "optim")) {
      ctl$optCtrl <- list()
    }
    t0 <- proc.time()[3]
    fits[i] <- list(tryCatch(
      suppressWarnings(fit_assembled(
        fit$spec, fit$frame, fit$bform, fit$call,
        REML = fit$REML, start = NULL, control = ctl, se = FALSE,
        lower = fit$lower, upper = fit$upper, priors = fit$priors,
        quadrature = isTRUE(fit$quadrature),
        data2 = fit$data2 %||% list()
      )),
      error = function(e) NULL
    ))
    times[i] <- proc.time()[3] - t0
  }
  names(fits) <- names(optimizers)
  structure(list(fits = fits, times = times), class = "frmtmb_allfit")
}

#' @export
print.frmtmb_allfit <- function(x, ...) {
  ok <- !vapply(x$fits, is.null, TRUE)
  tab <- data.frame(
    optimizer = names(x$fits),
    logLik = vapply(x$fits, function(f) {
      if (is.null(f)) NA_real_ else as.numeric(logLik(f))
    }, numeric(1)),
    convergence = vapply(x$fits, function(f) {
      if (is.null(f)) NA_integer_ else f$opt$convergence
    }, integer(1)),
    seconds = round(x$times, 2)
  )
  print(tab, row.names = FALSE)
  if (sum(ok) > 1L) {
    ll <- tab$logLik[ok]
    cat("\nlogLik spread:", format(diff(range(ll)), digits = 3), "\n")
    fe <- vapply(x$fits[ok], function(f) unlist(fixef(f)),
                 unlist(fixef(x$fits[ok][[1L]])))
    fe <- matrix(fe, ncol = sum(ok))
    cat("max fixed-effect spread:",
        format(max(apply(fe, 1, function(r) diff(range(r)))), digits = 3),
        "\n")
  }
  invisible(x)
}
