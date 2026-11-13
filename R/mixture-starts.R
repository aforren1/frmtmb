#' Refit a mixture from jittered starts and keep the best maximum at
#' which no component has degenerated.
#'
#' A family takes part by declaring `post$start_degenerate`, a function
#' of a fit and a response that returns its components at a degenerate
#' boundary in two kinds (`mixture_ord_start_class()` for an ordinal
#' mixture): `step`, a component that is a step function of its
#' predictors, the boundary that buys likelihood; and `boundary`, a flat
#' one (collapsed thresholds, a threshold no row is near). Maximum
#' likelihood for such a mixture often has its supremum at a step
#' function, so the best of several starts ranked by the likelihood alone
#' picks that point more often than it escapes it: on the 40 fits of
#' `dev/optima-ordmix-starts.R`, 12 of the 20 better points of 10 starts
#' were degenerate, every one a step function (dev/ordstart-findings.md).
#'
#' The default start's optimum is a candidate with the others, and only
#' runs the optimizer called converged are candidates. Each is proper, at
#' a flat boundary only, or a step function, and the choice is
#' `mixture_starts_pick()`'s:
#'
#' - A proper default is replaced only by a proper point that is better
#'   by more than `grad_tol` in log-likelihood, so that no fit the
#'   fit-end check passes gains its warning.
#' - A default at a flat boundary stays: the rule flags one of two points
#'   of one likelihood by how far the optimizer ran, so a move could take
#'   its warning away (cum_cum seed 10 under autoscale did).
#' - A step-function default is replaced by the best proper point, or by
#'   the best at a flat boundary when no start reached a proper one,
#'   whatever its likelihood: the step function is higher only because
#'   it lies on the boundary. Under brms 2.23.0's default priors the
#'   posterior density was 53 to 56 units lower at such points than at
#'   the proper maximum (cum_cum seed 4, cum_sratio seed 14,
#'   dev/ordstart-brms.R); brms's posterior mode is near the proper
#'   point, not at it.
#' - When a step function is above the point returned by more than
#'   `grad_tol`, a warning of class `frmtmb_mixture_step_start` says so.
#'
#' A default at a flat boundary is left as it is, and no starts run. A
#' default the optimizer did not call converged is left as it is too,
#' unless it is a step function (or its sharpening cannot be measured):
#' then the best converged proper start replaces it, when one exists. The
#' code-0 rule keeps unconverged runs out of the ranking; it was not meant
#' to protect an unconverged step function, and with it the unscaled and
#' autoscaled fits of cum_sratio 4, 14, 15 and 18 disagreed.
#'
#' The starts are drawn from a private stream with a named generator and
#' a fixed seed, so the fit neither reads nor changes the session's
#' random state, and the same call gives the same fit under any
#' `RNGkind()`. Another seed is another set of 10 starts, and on 11 of
#' the 40 fits another answer: the choice is the best of the points
#' these starts reach, not the best there is. Returns NULL when no starts
#' ran, and
#' otherwise `list(opt, moved)`: the optimum kept, with its record
#' `mixture_starts`, and whether it is a start's. The objective is left
#' at it.
#'
#' @noRd
mixture_starts <- function(fit, control, bounds, par_units = NULL,
                           tally = NULL, est_at, cold = TRUE,
                           verbose = 0L, center = NULL) {
  obj <- fit$obj
  resps <- Filter(function(r) is.function(r$family[["post"]][[
    "start_degenerate"]]), fit$spec$responses)
  n <- mixture_starts_count(control, cold, length(obj$env$random) > 0L)
  # A default the optimizer did not call converged is not a measured
  # maximum (on three cumulative components such points stood for
  # coefficients of 1e3 to 2e4, dev/reviews/2026-11-11-ordstart.md,
  # finding 3). It is replaced only when it is a step function and a
  # converged proper start exists; otherwise it stays the fit 0.81.0
  # returned
  converged <- isTRUE(fit$opt$convergence == 0)
  if (!length(resps) || !n || !length(obj$par) ||
        !all(is.finite(fit$opt$par)) || !is.finite(fit$opt$objective)) {
    return(NULL)
  }
  if (verbose) t0 <- vb_now()
  saved <- obj_state_save(obj)
  n_before <- if (is.null(tally)) 0L else tally$n
  np <- length(obj$par)
  lower <- rep_len(bounds$lower %||% -Inf, np)
  upper <- rep_len(bounds$upper %||% Inf, np)
  # A run from a start describes a point the fit may not report, and the
  # point it keeps is judged afterwards by the convergence and fit-end
  # checks, so a run's own warnings are not the fit's: nlminb's
  # "NA/NaN function evaluation" from a start where the objective is
  # undefined warned 10 times on one cs() mixture of
  # test-linpred-layers.R
  quiet <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
      invokeRestart("muffleWarning")
    })
  }
  # each response's degenerate components at `o`, by kind; NULL when a
  # family's check fails there, so that an unread point is never taken
  # for a proper one
  degen_at <- function(o) {
    f <- fit
    # parList() reads the random effects' modes from the last point
    if (length(obj$env$random)) quiet(obj$fn(o$par))
    f$opt <- o
    f$estimates <- est_at(o)
    out <- lapply(resps, function(r) {
      tryCatch(quiet(r$family[["post"]][["start_degenerate"]](f,
                                                             r$resp_name)),
               error = function(e) NULL)
    })
    if (any(vapply(out, is.null, NA))) NULL else out
  }
  # 2 with a step function anywhere, 1 with only flat boundaries, 0
  tier_of <- function(dg) {
    if (any(lengths(lapply(dg, `[[`, "step")) > 0L)) return(2L)
    if (any(lengths(lapply(dg, `[[`, "boundary")) > 0L)) return(1L)
    0L
  }
  cands <- list(fit$opt)
  degs <- list(degen_at(fit$opt))
  # a default at a flat boundary stays (mixture_starts_pick()), so its
  # starts would buy nothing
  if (is.null(degs[[1L]]) || tier_of(degs[[1L]]) == 1L ||
        (!converged && tier_of(degs[[1L]]) != 2L)) {
    obj_settle(obj, fit$opt, fit$opt, saved)
    return(NULL)
  }
  Z <- mixture_starts_draws(np, n)
  u <- mixture_starts_units(fit$frame, names(obj$par))
  # the point the jitter is about: the fit's own start, or the one a
  # caller names (the plain fit's start under an autoscale warm start)
  p0 <- if (is.function(center)) {
    tryCatch(center(), error = function(e) NULL)
  } else obj$par
  if (length(p0) != np || !all(is.finite(p0))) p0 <- obj$par
  for (i in seq_len(n)) {
    p <- pmin(pmax(p0 + mixture_starts_sd * u * Z[, i], lower), upper)
    v0 <- tryCatch(quiet(obj$fn(p)), error = function(e) NA_real_)
    if (!is.finite(v0)) next
    o <- tryCatch(quiet(optimize_obj(obj, control, bounds, par_units,
                                     start_par = p, tally = tally)),
                  error = function(e) NULL)
    # only a run the optimizer called converged is a maximum to rank; a
    # stop at the evaluation limit won on one three-component fit and
    # turned its code 0 into 1
    if (is.null(o) || !isTRUE(o$convergence == 0) ||
          !all(is.finite(o$par)) || !is.finite(o$objective)) {
      next
    }
    dg <- degen_at(o)
    if (is.null(dg)) next
    cands[[length(cands) + 1L]] <- o
    degs[[length(degs) + 1L]] <- dg
  }
  obj_v <- vapply(cands, `[[`, 0, "objective")
  tier <- vapply(degs, tier_of, 0L)
  tol <- control$grad_tol %||% 1e-3
  keep <- mixture_starts_pick(obj_v, tier, tol, converged)
  best <- cands[[keep]]
  above <- which(tier == 2L & obj_v < obj_v[keep] - tol)
  rec <- c(starts = n, reached = length(cands) - 1L,
           gain = obj_v[1L] - obj_v[keep],
           degenerate_above = length(above), tier = tier[keep])
  if (length(above) && tier[keep] < 2L) {
    a <- above[which.min(obj_v[above])]
    mixture_starts_warn(fit, lapply(degs[[a]], `[[`, "step"), -obj_v[a],
                        -obj_v[keep], n, length(above), keep == 1L)
  }
  if (verbose) {
    vb_stage("mixture starts", t0,
             paste0(n, " start", if (n != 1L) "s", ", ",
                    sum(tier[-1L] < 2L), " without a step function, gain ",
                    format(rec[["gain"]], digits = 3)))
  }
  obj_settle(obj, fit$opt, best, saved)
  if (keep > 1L) {
    # a start's run replaces the optimum, whose records would go with it
    for (nm in c("stationary_escape", "mo_search")) {
      if (is.null(best[[nm]])) best[[nm]] <- fit$opt[[nm]]
    }
  }
  # recorded whether or not the optimum moved: what the starts cost
  best$mixture_starts <- rec
  # the total over every run of the objective, as escape_stationary()
  if (!is.null(tally) && !is.null(best$nonfinite_trials)) {
    best$nonfinite_trials <- (fit$opt$nonfinite_trials %||% 0L) +
      max(tally$n - n_before, 0L)
  }
  if (!is.null(tally)) best$evals <- tally$evals
  list(opt = best, moved = keep > 1L)
}

#' Which candidate `mixture_starts()` keeps: the index into `obj_v`, the
#' candidates' objectives with the default's first, given each one's
#' `tier` (0 proper, 1 a flat boundary only, 2 a step function) and the
#' gain `tol` that counts as better.
#'
#' A proper default is replaced only by a better point that is proper
#' too, so no fit the fit-end check passes gains its warning: on
#' cum_sratio seed 10 a start 0.23 higher had a threshold no row is near.
#' A flat-boundary default stays, so that its warning does: on cum_cum
#' seed 10 under autoscale the better point (+2.78) had a gap of 1e-10
#' where the plain fit's had one of 0, and the warning went. A step
#' function is left for the best point without one, whatever its
#' likelihood. A default the optimizer did not call converged
#' (`converged = FALSE`) gives way only if it is a step function, and
#' only to the best proper point.
#'
#' @noRd
mixture_starts_pick <- function(obj_v, tier, tol, converged = TRUE) {
  if (tier[1L] == 1L) return(1L)
  # an unconverged default: only a step function gives way, and only to
  # a converged proper point (every candidate after it converged)
  if (!converged) {
    if (tier[1L] != 2L) return(1L)
    pool <- which(tier == 0L)
    if (!length(pool)) return(1L)
    return(pool[which.min(obj_v[pool])])
  }
  pool <- if (tier[1L] == 0L) which(tier == 0L) else
    if (any(tier == 0L)) which(tier == 0L) else which(tier == 1L)
  if (!length(pool)) return(1L)
  j <- pool[which.min(obj_v[pool])]
  if (tier[1L] == 2L || obj_v[1L] - obj_v[j] > tol) j else 1L
}

#' The standard deviation of a start's jitter, in each coordinate's unit
#' (`mixture_starts_units()`): the 0.5 of dev/optima-ordmix-starts.R.
#'
#' @noRd
mixture_starts_sd <- 0.5

#' The number of jittered starts: `starts` when it is a number, and under
#' the default (`NULL`) `mixture_starts_default` on a cold fit without
#' random effects. None on a fit given a `start` or a refit's warm start,
#' which already says where to look, and none with random effects, where
#' each start's Laplace fit made the starts cost 15 to 23 times the fit
#' (n = 1000, 50 groups, dev/ordstart-cost.R) against 5 to 10 times
#' without.
#'
#' @noRd
mixture_starts_count <- function(control, cold, random = FALSE) {
  n <- control[["starts"]]
  if (is.null(n)) {
    return(if (isTRUE(cold) && !random) mixture_starts_default else 0L)
  }
  as.integer(n)
}

#' The default number of jittered starts (see dev/ordstart-findings.md
#' for the count measured against 40).
#'
#' @noRd
mixture_starts_default <- 10L

#' Standard normal draws for `n` starts of `np` coordinates, from a
#' private stream: a fit must not depend on or disturb the session's
#' random state (as imp_draws()).
#'
#' The generator is named, not inherited: under the session's kind, a
#' session with `RNGkind("L'Ecuyer-CMRG")`, which parallel work sets, drew
#' other starts and got another answer on 11 of the 40 fits, by up to
#' 10.4 in log-likelihood (dev/reviews/2026-11-11-ordstart.md, finding
#' 1). The caller's kinds (all of them, however many the R version has)
#' and seed are put back, and a seed that did not exist is removed again.
#'
#' @noRd
mixture_starts_draws <- function(np, n) {
  ge <- globalenv()
  kinds <- RNGkind()
  had <- exists(".Random.seed", envir = ge, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = ge)
  on.exit({
    # setting a kind seeds the stream anew, so the seed goes back after;
    # a "Rounding" sample.kind warns that it is non-uniform, which the
    # caller chose and has been told. Every kind RNGkind() reports goes
    # back, by position: R-devel added a fourth (binom.kind) after the
    # three of R 4.6
    suppressWarnings(do.call(RNGkind, as.list(kinds)))
    if (had) {
      assign(".Random.seed", old, envir = ge)
    } else if (exists(".Random.seed", envir = ge, inherits = FALSE)) {
      rm(".Random.seed", envir = ge)
    }
  })
  set.seed(20261111L, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  matrix(stats::rnorm(np * n), np, n)
}

#' The unit of each outer coordinate for the jitter: one over the root
#' mean square of its design column for a coefficient, so that a step
#' moves the linear predictor by about the same amount whatever the
#' column's scale and offset, and 1 for every other coordinate.
#'
#' @noRd
mixture_starts_units <- function(frame, pn) {
  u <- rep(1, length(pn))
  for (lp in frame[["linpreds"]]) {
    X <- lp[["X"]]
    if (is.null(X) || !ncol(X) || !is.null(lp[["constant"]])) next
    idx <- outer_index_of(lp, pn, frame)
    if (is.null(idx)) next
    p <- min(ncol(X), length(idx))
    rms <- sqrt(Matrix::colMeans(X[, seq_len(p), drop = FALSE]^2))
    ok <- is.finite(rms) & rms > 0
    u[idx[seq_len(p)][ok]] <- 1 / rms[ok]
  }
  u
}

#' The warning that a start reached a higher point where a component is
#' degenerate.
#'
#' @noRd
mixture_starts_warn <- function(fit, deg, ll_deg, ll_kept, n, n_above,
                                kept_default) {
  deg <- deg[lengths(deg) > 0L]
  fams <- vapply(names(deg), function(r) {
    fit$spec$responses[[r]]$family[["family"]]
  }, "")
  what <- paste(vapply(seq_along(deg), function(i) {
    many <- length(deg[[i]]) > 1L
    paste0("in ", fams[i], if (length(fit$spec$responses) > 1L) {
      paste0(" (response ", names(deg)[i], ")")
    }, ", component", if (many) "s", " ",
    paste(deg[[i]], collapse = " and "),
    if (many) " are step functions of their predictors" else
      " is a step function of its predictors")
  }, ""), collapse = "; ")
  # under a prior the optimizer ranks the log posterior, not logLik()
  what_ll <- if (length(fit$resolved_prior)) {
    "log-likelihood plus log prior"
  } else {
    "log-likelihood"
  }
  mixture_starts_signal(paste0(
    n_above, " of ", n, " jittered starts reached a higher ", what_ll,
    " than the fit reports, the highest ",
    format(round(ll_deg, 4), nsmall = 4), " against ",
    format(round(ll_kept, 4), nsmall = 4), ", where ", what, ". That ",
    "point is a degenerate boundary, where maximum likelihood for an ",
    "ordinal mixture often has its supremum, so the fit reports the best ",
    "converged point ", if (kept_default) "(the default start's) ",
    "without one. logLik() and the information criteria are those of this ",
    "point. Priors hold the components, on class Intercept and b with ",
    "dpar = \"mu<k>\" (brms's default prior on b is flat); ",
    "frmtmb_control(starts = 0) turns the starts off"))
}

#' How the step-start finding reaches the caller. A warning or a message
#' is the user's call (dev/ordstart-findings.md); the class is the same
#' either way, so a switch changes this function alone.
#'
#' @noRd
mixture_starts_signal <- function(msg) {
  frm_warning(msg, call. = FALSE, class = "frmtmb_mixture_step_start")
}
