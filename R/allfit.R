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
#' Every optimizer starts from the same point, as in lme4: by default
#' the original fit's estimates (lme4's `start_from_mle = TRUE`), or,
#' with `start_from_mle = FALSE`, the point the original fit started
#' from, which includes its `start =` and the starting values a located
#' prior places on a nonlinear parameter. A nonlinear model is rarely
#' defined at the default start of zero, so refits from there fail or
#' stop somewhere else, and their disagreement then says nothing about
#' the original fit.
#'
#' The printed table compares each refit with the best log-likelihood
#' any fit reached, the original included. An optimizer that reports
#' convergence more than `grad_tol` (see [frmtmb_control()]) below that
#' best is marked as having stopped elsewhere, so a success code at a
#' worse optimum does not read as agreement. One that reaches the best
#' while reporting a failure code agrees, with its code shown: started
#' at the optimum, `optim`'s line search can fail there (code 52).
#'
#' @param fit A `frmtmb_fit`.
#' @param optimizers Named list of optimizers (names or functions, as
#'   in [frmtmb_control()]). Default: everything available.
#' @param start_from_mle If `TRUE` (default), start every optimizer at
#'   the original fit's estimates; if `FALSE`, where the original fit
#'   started.
#' @param ... Refused: an argument the method does not have is an
#'   error naming it, rather than silently changing nothing.
#' @return A `frmtmb_allfit` object: `$fits` (the refits, `NULL` where
#'   one errored), `$errors` (the error message of each refit that
#'   failed, `NA` elsewhere), `$times`, and `$original`, the
#'   log-likelihood and convergence code of the fit that was checked.
#'   Printing it compares log-likelihoods, convergence and fixed-effect
#'   spread.
#' @examples
#' set.seed(6)
#' dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
#' dd$y <- rpois(60, exp(0.3 + 0.4 * dd$x + rnorm(6, 0, 0.4)[dd$g]))
#' fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
#' frm_allfit(fit)
#' @export
frm_allfit <- function(fit, optimizers = NULL, start_from_mle = TRUE, ...) {
  frm_check_dots(...)
  check_flag(start_from_mle, "start_from_mle")
  optimizers <- optimizers %||% allfit_optimizers()
  ctl0 <- fit$control %||% frmtmb_control()
  # print.frmtmb_allfit() already reports the per-optimizer timings, so
  # a verbose original fit must not make every refit narrate itself
  ctl0$verbose <- FALSE
  fits <- list()
  errors <- rep(NA_character_, length(optimizers))
  times <- numeric(length(optimizers))
  for (i in seq_along(optimizers)) {
    ctl <- ctl0
    ctl$optimizer <- optimizers[[i]]
    if (!identical(ctl$optimizer, "nlminb") &&
        !identical(ctl$optimizer, "optim")) {
      ctl$optCtrl <- list()
    }
    t0 <- proc.time()[3]
    r <- tryCatch(
      suppressWarnings(fit_assembled(
        fit$spec, fit$frame, fit$bform, fit$call,
        REML = fit$REML, start = if (!start_from_mle) fit$start,
        control = ctl, se = FALSE,
        lower = fit$lower, upper = fit$upper, prior = fit$prior,
        quadrature = isTRUE(fit$quadrature),
        importance = fit$importance$draws %||% 0L,
        template = if (start_from_mle) fit$estimates,
        data2 = fit$data2 %||% list()
      )),
      error = function(e) e
    )
    if (inherits(r, "error")) {
      errors[i] <- conditionMessage(r)
      r <- NULL
    }
    fits[i] <- list(r)
    times[i] <- proc.time()[3] - t0
  }
  names(fits) <- names(errors) <- names(optimizers)
  structure(list(fits = fits, errors = errors, times = times,
                 original = list(logLik = as.numeric(logLik(fit)),
                                 convergence = fit$opt$convergence),
                 tol = ctl0$grad_tol %||% 1e-3,
                 start_from_mle = start_from_mle),
            class = "frmtmb_allfit")
}

#' @export
print.frmtmb_allfit <- function(x, ...) {
  frm_check_dots(...)
  ok <- !vapply(x$fits, is.null, TRUE)
  ll <- vapply(x$fits, function(f) {
    if (is.null(f)) NA_real_ else as.numeric(logLik(f))
  }, numeric(1))
  code <- vapply(x$fits, function(f) {
    if (is.null(f)) NA_integer_ else as.integer(f$opt$convergence)
  }, integer(1))
  # an object from before `original` was recorded has no reference row
  ll0 <- x$original$logLik %||% NA_real_
  best <- max(c(ll, ll0), na.rm = TRUE)
  tol <- x$tol %||% 1e-3
  err <- x$errors %||% rep(NA_character_, length(x$fits))
  # an optimizer started AT the optimum can fail its own line search
  # there (optim's code 52) while standing on the best point, so the
  # log-likelihood decides and the code is shown beside it
  below <- ll < best - tol
  verdict <- ifelse(!ok, paste("failed:", substr(err, 1L, 60L)),
             ifelse(below & code != 0L, "did not converge",
             ifelse(below, "converged elsewhere",
             ifelse(code != 0L, paste0("agrees; code ", code), "agrees"))))
  if (!is.na(ll0)) {
    cat("original fit: logLik ", format(ll0, digits = 10), ", code ",
        x$original$convergence, if (ll0 < best - tol) {
          paste0(", ", format(best - ll0, digits = 3), " below the best")
        }, "\n", sep = "")
  }
  cat("refits started at ", if (isFALSE(x$start_from_mle)) {
    "the original fit's start"
  } else "the original fit's estimates", "\n\n", sep = "")
  tab <- data.frame(
    optimizer = names(x$fits),
    logLik = ll,
    `vs best` = signif(ll - best, 3),
    convergence = code,
    verdict = verdict,
    seconds = round(x$times, 2),
    check.names = FALSE
  )
  print(tab, row.names = FALSE)
  if (any(verdict == "converged elsewhere")) {
    cat("\nAn optimizer that reports convergence below the best ",
        "log-likelihood stopped at a different point: the likelihood ",
        "has more than one optimum, or that optimizer stops early ",
        "on this model.\n", sep = "")
  }
  if (sum(ok) > 1L) {
    cat("\nlogLik spread:", format(diff(range(ll[ok])), digits = 3), "\n")
    # the flat vector, not the summary matrix: what is compared across
    # restarts is the estimates, one per coefficient
    fe <- vapply(x$fits[ok], function(f) fixef(f, flatten = TRUE),
                 fixef(x$fits[ok][[1L]], flatten = TRUE))
    fe <- matrix(fe, ncol = sum(ok))
    cat("max fixed-effect spread:",
        format(max(apply(fe, 1, function(r) diff(range(r)))), digits = 3),
        "\n")
  }
  invisible(x)
}
