#' Refit a mixture from other starts and keep the best maximum at which
#' no component has degenerated: jittered starts when `starts` asks for
#' them, and the swapped starts of a mixture whose components are of
#' different families.
#'
#' A family takes part by declaring `post$start_degenerate`, a function
#' of a fit and a response that returns its components at a degenerate
#' boundary in two kinds (`mixture_ord_start_class()` for an ordinal
#' mixture): `step`, a component that is a step function of its
#' predictors, the boundary that buys likelihood; and `boundary`, a flat
#' one (collapsed thresholds, a threshold no row is near, a threshold
#' without a standard error or with one above the degenerate cut).
#' Maximum likelihood for such a mixture often has its supremum at a step
#' function, so the best of several starts ranked by the likelihood alone
#' picks that point more often than it escapes it: on the 40 fits of
#' `dev/optima-ordmix-starts.R`, 12 of the 20 better points of 10 starts
#' were degenerate, every one a step function (dev/ordstart-findings.md).
#'
#' Jittered starts (`starts = n`): the default start's optimum is a
#' candidate with the others, and only runs the optimizer called
#' converged are candidates. Each is proper, at a flat boundary only, or
#' a step function, and the choice is `mixture_starts_pick()`'s:
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
#'   `grad_tol`, a message of class `frmtmb_mixture_step_start` says so.
#'
#' A default at a flat boundary runs no jittered starts. A default the
#' optimizer did not call converged runs none either, unless it is a step
#' function (or its sharpening cannot be measured): then the best
#' converged proper start replaces it, when one exists.
#'
#' Swapped starts (`mixture_swap_candidates()`): on a fit from scratch,
#' unless `starts = 0`, one run from the default optimum with the
#' components' parameters moved between components, for each assignment
#' whose start is another model (`mixture_swap_distinct()`). Without
#' jittered starts a swapped run is kept only when it converged, is
#' proper, and is above the default by more than `grad_tol` per unit of
#' mean case weight; with them, it is one more candidate of the pick. A
#' kept run returns the best proper end of its check
#' (`mixture_swap_settle()`), and is said in a message of class
#' `frmtmb_mixture_swap`. The default
#' runs no jittered starts (backlog 6136, the user's decision of
#' 2026-11-17): with random effects they cost a median 18 times the fit,
#' without 5.6 to 7.1 times. With the swap, the silent shortfalls of 320
#' fits with random effects fell from 19 to 3, no fit moved down, and the
#' fit cost a median 2.0 to 3.8 times the fit without it on mixtures of
#' two families (6.9 times where the swap was kept), 1.0 on mixtures of
#' one (dev/ordstarts-sum6.R, dev/ordstarts-findings.md).
#'
#' The jitter is drawn from a private stream with a named generator and a
#' fixed seed, so the fit neither reads nor changes the session's random
#' state, and the same call gives the same fit under any `RNGkind()`.
#' Returns NULL when nothing ran, and otherwise `list(opt, moved)`: the
#' optimum kept, with its record `mixture_starts`, and whether it is
#' another run's. The objective is left at it.
#'
#' @noRd
mixture_starts <- function(fit, control, bounds, par_units = NULL,
                           tally = NULL, est_at, cold = TRUE,
                           verbose = 0L, center = NULL) {
  obj <- fit$obj
  resps <- Filter(function(r) is.function(r$family[["post"]][[
    "start_degenerate"]]), fit$spec$responses)
  n <- mixture_starts_count(control, cold)
  # A default the optimizer did not call converged is not a measured
  # maximum (on three cumulative components such points stood for
  # coefficients of 1e3 to 2e4, dev/reviews/2026-11-11-ordstart.md,
  # finding 3). The jittered starts replace it only when it is a step
  # function and a converged proper start exists
  converged <- isTRUE(fit$opt$convergence == 0)
  if (!length(resps) || !length(obj$par) ||
        !all(is.finite(fit$opt$par)) || !is.finite(fit$opt$objective)) {
    return(NULL)
  }
  swaps <- if (isTRUE(cold) && !identical(control[["starts"]], 0L)) {
    mixture_swap_distinct(fit, mixture_swap_candidates(fit, resps))
  } else list()
  if (!n && !length(swaps)) return(NULL)
  if (verbose) t0 <- vb_now()
  saved <- obj_state_save(obj)
  n_before <- if (is.null(tally)) 0L else tally$n
  np <- length(obj$par)
  lower <- rep_len(bounds$lower %||% -Inf, np)
  upper <- rep_len(bounds$upper %||% Inf, np)
  # a gain counts per unit of mean case weight, as the sharpening does:
  # weights 1/1000 divide every gain by 1000, and the swap went on 5 of 6
  # fits whose unweighted gain was 0.14 to 0.93 (review of 2026-11-18,
  # item 2)
  tol <- (control$grad_tol %||% 1e-3) * mixture_starts_weight(fit, resps)
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
  # for a proper one. A cache of its own: the class builds the
  # covariance at the point it reads, and the fit's cache holds the
  # default's
  degen_at <- function(o) {
    f <- fit
    # parList() reads the random effects' modes from the last point
    if (length(obj$env$random)) quiet(obj$fn(o$par))
    f$opt <- o
    f$estimates <- est_at(o)
    f$cache <- new.env(parent = emptyenv())
    out <- lapply(resps, function(r) {
      tryCatch(quiet(r$family[["post"]][["start_degenerate"]](f,
                                                             r$resp_name)),
               error = function(e) NULL)
    })
    if (any(vapply(out, is.null, NA))) NULL else out
  }
  # 2 with a step function anywhere, 1 with only flat boundaries (by the
  # thresholds' travel or by their standard errors), 0
  tier_of <- function(dg) {
    if (any(lengths(lapply(dg, `[[`, "step")) > 0L)) return(2L)
    if (any(lengths(lapply(dg, `[[`, "boundary")) > 0L) ||
          any(lengths(lapply(dg, `[[`, "unplaced")) > 0L)) {
      return(1L)
    }
    0L
  }
  run_from <- function(p) {
    p <- pmin(pmax(p, lower), upper)
    v0 <- tryCatch(quiet(obj$fn(p)), error = function(e) NA_real_)
    if (!is.finite(v0)) return(NULL)
    o <- tryCatch(quiet(optimize_obj(obj, control, bounds, par_units,
                                     start_par = p, tally = tally)),
                  error = function(e) NULL)
    # only a run the optimizer called converged is a maximum to rank; a
    # stop at the evaluation limit won on one three-component fit and
    # turned its code 0 into 1
    if (is.null(o) || !isTRUE(o$convergence == 0) ||
          !all(is.finite(o$par)) || !is.finite(o$objective)) {
      return(NULL)
    }
    o
  }
  proper_at <- function(o) {
    dg <- degen_at(o)
    !is.null(dg) && tier_of(dg) == 0L
  }
  cands <- list(fit$opt)
  degs <- list(degen_at(fit$opt))
  kind <- "default"
  if (is.null(degs[[1L]])) {
    obj_settle(obj, fit$opt, fit$opt, saved)
    return(NULL)
  }
  # A default flat by its thresholds' travel (run off, collapsed) stays
  # under the pick (mixture_starts_pick()), so jittered starts would buy
  # nothing there. One flat only by a threshold's standard error, or a
  # lost one, does not: there the starts ran on none of 30 of 120 fits
  # whose warning named them (review of 2026-11-18, item 1)
  t1 <- tier_of(degs[[1L]])
  travel <- t1 == 1L &&
    any(lengths(lapply(degs[[1L]], `[[`, "boundary")) > 0L)
  run_jit <- n > 0L && !travel && (converged || t1 == 2L)
  if (!run_jit && !length(swaps)) {
    obj_settle(obj, fit$opt, fit$opt, saved)
    return(NULL)
  }
  if (run_jit) {
    Z <- mixture_starts_draws(np, n)
    u <- mixture_starts_units(fit$frame, names(obj$par))
    # the point the jitter is about: the fit's own start, or the one a
    # caller names (the plain fit's start under an autoscale warm start)
    p0 <- if (is.function(center)) {
      tryCatch(center(), error = function(e) NULL)
    } else obj$par
    if (length(p0) != np || !all(is.finite(p0))) p0 <- obj$par
    for (i in seq_len(n)) {
      o <- run_from(p0 + mixture_starts_sd * u * Z[, i])
      if (is.null(o)) next
      dg <- degen_at(o)
      if (is.null(dg)) next
      cands[[length(cands) + 1L]] <- o
      degs[[length(degs) + 1L]] <- dg
      kind <- c(kind, "start")
    }
  }
  swap_label <- character(length(cands))
  for (sw in swaps) {
    mixture_swap_inner(obj, sw$b)
    o <- run_from(sw$par)
    if (is.null(o)) next
    dg <- degen_at(o)
    if (is.null(dg)) next
    cands[[length(cands) + 1L]] <- o
    degs[[length(degs) + 1L]] <- dg
    kind <- c(kind, "swap")
    swap_label <- c(swap_label, sw$label)
  }
  obj_v <- vapply(cands, `[[`, 0, "objective")
  tier <- vapply(degs, tier_of, 0L)
  keep <- if (run_jit) {
    mixture_starts_pick(obj_v, tier, tol, converged, sticky = travel)
  } else {
    # the swap alone: a higher converged proper run, else the default
    up <- which(kind == "swap" & tier == 0L & obj_v < obj_v[1L] - tol)
    if (length(up)) up[which.min(obj_v[up])] else 1L
  }
  if (kind[keep] == "swap") {
    # the best proper end of its jittered runs (mixture_swap_settle())
    cands[[keep]] <- mixture_swap_settle(cands[[keep]], run_from, np, tol,
                                         proper_at)
    obj_v[keep] <- cands[[keep]]$objective
  }
  best <- cands[[keep]]
  n_swaps <- sum(kind == "swap")
  above <- which(tier == 2L & obj_v < obj_v[keep] - tol)
  rec <- c(starts = if (run_jit) n else 0L, swaps = length(swaps),
           reached = length(cands) - 1L, gain = obj_v[1L] - obj_v[keep],
           degenerate_above = length(above), tier = tier[keep],
           swapped = as.numeric(kind[keep] == "swap"))
  if (length(above) && tier[keep] < 2L) {
    a <- above[which.min(obj_v[above])]
    mixture_starts_warn(fit, lapply(degs[[a]], `[[`, "step"), -obj_v[a],
                        -obj_v[keep], if (run_jit) n else 0L,
                        length(swaps), length(above), keep == 1L)
  }
  if (kind[keep] == "swap") {
    mixture_swap_signal(paste0(
      "The components of ", paste(vapply(resps, function(r) {
        r$family[["family"]]
      }, ""), collapse = ", "), " were refitted from the default fit ",
      "with their parameters swapped (", swap_label[keep], "), and that ",
      "fit is higher, log-likelihood ",
      format(round(-obj_v[keep], 4), nsmall = 4), " against ",
      format(round(-obj_v[1L], 4), nsmall = 4), ", and passes the ",
      "degenerate test, so frm() returns it. ",
      "frmtmb_control(starts = 0) turns the swap off"))
  }
  if (verbose) {
    vb_stage("mixture starts", t0,
             paste0(if (run_jit) n else 0L, " jittered, ", n_swaps,
                    " swapped, ", sum(tier[-1L] < 2L),
                    " without a step function, gain ",
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

#' Which candidate `mixture_starts()` keeps among the default and the
#' jittered (and swapped) runs: the index into `obj_v`, the candidates'
#' objectives with the default's first, given each one's `tier` (0
#' proper, 1 a flat boundary only, 2 a step function) and the gain `tol`
#' that counts as better.
#'
#' A proper default is replaced only by a better point that is proper
#' too, so no fit the fit-end check passes gains its warning: on
#' cum_sratio seed 10 a start 0.23 higher had a threshold no row is near.
#' A default flat by its thresholds' travel (`sticky`: run off,
#' collapsed) stays, so that its warning does: on cum_cum seed 10 under
#' autoscale the better point (+2.78) had a gap of 1e-10 where the plain
#' fit's had one of 0, and the warning went. A default flat only by a
#' threshold's standard error, or a lost one, is no such flag of the
#' optimizer's travel, and a better proper point replaces it as it does
#' a proper default (test-compat.R's sum_to_zero cell: 1.77 higher, no
#' warning). A step function is left for the best point without one,
#' whatever its likelihood. A default the optimizer did not call
#' converged (`converged = FALSE`) gives way only if it is a step
#' function, and only to the best proper point.
#'
#' @noRd
mixture_starts_pick <- function(obj_v, tier, tol, converged = TRUE,
                                sticky = TRUE) {
  if (tier[1L] == 1L && sticky) return(1L)
  # an unconverged default: only a step function gives way, and only to
  # a converged proper point (every candidate after it converged)
  if (!converged) {
    if (tier[1L] != 2L) return(1L)
    pool <- which(tier == 0L)
    if (!length(pool)) return(1L)
    return(pool[which.min(obj_v[pool])])
  }
  pool <- if (tier[1L] <= 1L) setdiff(which(tier == 0L), 1L) else
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

#' The number of jittered starts: `starts` when it is a number, and none
#' under the default (`NULL`), on any fit. The default ran 10 on a fit
#' from scratch without random effects until the user's decision of
#' 2026-11-17 (backlog 6136): they cost 5.6 to 7.1 times the fit there
#' and a median 18 times with random effects, and a uniform default
#' keeps the cost the same for every model form. A number runs that many
#' on every fit, a refit's warm-started one too.
#'
#' @noRd
mixture_starts_count <- function(control, cold) {
  n <- control[["starts"]]
  if (is.null(n)) return(0L)
  as.integer(n)
}

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
mixture_starts_draws <- function(np, n, seed = 20261111L) {
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
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
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

#' The message that a start reached a higher point where a component is
#' degenerate.
#'
#' @noRd
mixture_starts_warn <- function(fit, deg, ll_deg, ll_kept, n, n_swaps,
                                n_above, kept_default) {
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
  runs <- if (n && n_swaps) {
    paste0(n, " jittered and ", n_swaps, " swapped starts")
  } else if (n) {
    paste0(n, " jittered starts")
  } else {
    paste0(n_swaps, " swapped start", if (n_swaps != 1L) "s")
  }
  mixture_starts_signal(paste0(
    n_above, " of ", runs, " reached a higher ", what_ll,
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

#' How the step-start finding reaches the caller: a message, by the
#' user's decision of 2026-11-13 (dev/ordstart-findings.md had left it
#' to the user). The point returned is the better answer, so nothing
#' about it calls for a warning; the class stays for a handler.
#'
#' @noRd
mixture_starts_signal <- function(msg) {
  frm_message(msg, class = "frmtmb_mixture_step_start")
}

#' How a kept swap reaches the caller: a message of class
#' `frmtmb_mixture_swap`, as the step-start finding is one. The point
#' returned is the better answer, a higher maximum that passes the
#' degenerate test, so nothing about it calls for a warning.
#'
#' @noRd
mixture_swap_signal <- function(msg) {
  frm_message(msg, class = "frmtmb_mixture_swap")
}

#' The point a kept swap returns: the best of its optimum `o` and four
#' runs from it, every outer coordinate moved by a normal draw with
#' standard deviation 0.1 (a private stream of its own seed), counting a
#' run's end only when it converged, is above the best so far by more
#' than `tol` and is proper (`proper(end)`). This is the check of
#' dev/ordstarts-refits.R ("real"), which held on 45 of the 46 points the
#' jittered starts had moved to with random effects and found the one
#' that is not a maximum (cum_sratio n = 300 seed 40, a group sd running
#' to zero, 0.0099 more) only with its fourth run, so fewer runs would
#' have missed it. A climbing end is the better answer, not a reason to
#' return the lower default; an end at a step function is the supremum
#' the pick exists to avoid, and is not counted (review of 2026-11-18,
#' item 5). Its cost falls only on a fit whose swap is kept.
#'
#' @noRd
mixture_swap_settle <- function(o, run_from, np, tol, proper) {
  Z <- mixture_starts_draws(np, mixture_swap_jitter_runs, seed = 20261117L)
  best <- o
  for (j in seq_len(ncol(Z))) {
    r <- run_from(o$par + mixture_swap_jitter_sd * Z[, j])
    if (!is.null(r) && r$objective < best$objective - tol && proper(r)) {
      best <- r
    }
  }
  best
}

#' The mean case weight of the responses `resps` of `fit`, the unit in
#' which `mixture_starts()` counts a gain: constant weights scale the
#' objective, not the fit.
#'
#' @noRd
mixture_starts_weight <- function(fit, resps) {
  w <- vapply(resps, function(r) {
    v <- fit$frame[["aterm_values"]][[r$resp_name]][["weights"]]
    if (is.null(v)) 1 else mean(as.numeric(v))
  }, 0)
  m <- mean(w)
  if (is.finite(m) && m > 0) m else 1
}

#' The swapped starts that are another model: each start's objective is
#' evaluated once, and a start is dropped when it equals the fit's, or an
#' earlier start's, to rounding. Permuting exchangeable components is
#' label switching: the objective at the permuted point differed from the
#' fit's by at most 2.8e-9 on 78 cum_cum fits (dev/ordstarts-sum4.R),
#' and by none on sratio against cratio under the logit, two names of one
#' model, while a prior, an offset or a fixed parameter on one component
#' moved it by 7.7 to 119 (review of 2026-11-18, item 3). The objective is
#' the definition, so it decides, not the families' names. The objective
#' is left as it was.
#'
#' @noRd
mixture_swap_distinct <- function(fit, swaps) {
  if (!length(swaps)) return(swaps)
  obj <- fit$obj
  saved <- obj_state_save(obj)
  on.exit(obj_state_restore(obj, saved), add = TRUE)
  f0 <- fit$opt$objective
  seen <- f0
  out <- list()
  for (s in swaps) {
    mixture_swap_inner(obj, s$b)
    v <- tryCatch(suppressWarnings(as.numeric(obj$fn(s$par))),
                  error = function(e) NA_real_)
    if (!is.finite(v)) next
    eps <- sqrt(.Machine$double.eps) * max(1, abs(f0))
    if (any(abs(v - seen) <= eps)) next
    seen <- c(seen, v)
    out[[length(out) + 1L]] <- s
  }
  out
}

#' Start the inner (Laplace) problem of the next evaluation of `obj` at
#' the swapped random-effect modes `b`, so that the objective at a
#' swapped start does not climb from the other component's modes: TMB
#' starts the inner optimization at `last.par.best[random]`. Nothing
#' without random effects or without `b`.
#'
#' @noRd
mixture_swap_inner <- function(obj, b) {
  rnd <- obj$env$random
  lp <- obj$env$last.par.best
  if (!length(rnd) || is.null(b) || length(b) != length(rnd) ||
        is.null(lp)) {
    return(invisible(NULL))
  }
  lp[rnd] <- b
  assign("last.par.best", lp, envir = obj$env)
  invisible(NULL)
}

#' The runs and the jitter of `mixture_swap_settle()`.
#'
#' @noRd
mixture_swap_jitter_runs <- 4L

#' @rdname mixture_swap_jitter_runs
#' @noRd
mixture_swap_jitter_sd <- 0.1

#' The swapped starts of every response whose family declares
#' `post$start_swaps`: a list of `list(par, label)`, one outer vector per
#' assignment, the response named in the label of a fit with more than
#' one.
#'
#' @noRd
mixture_swap_candidates <- function(fit, resps) {
  out <- list()
  for (r in resps) {
    fn <- r$family[["post"]][["start_swaps"]]
    if (!is.function(fn)) next
    sw <- tryCatch(fn(fit, r$resp_name), error = function(e) list())
    for (s in sw) {
      if (length(fit$spec$responses) > 1L) {
        s$label <- paste0("response ", r$resp_name, ": ", s$label)
      }
      out[[length(out) + 1L]] <- s
    }
  }
  out
}

#' The swapped starts of an ordinal mixture: from the fit's optimum, one
#' outer vector per assignment of components to components, as
#' `list(par, label, sigma)`.
#'
#' Assignment `sigma` gives component k the parameters component
#' `sigma[k]` has at the optimum: its coefficients of every dpar
#' (`mu<k>_*`, `disc<k>_*`), its random-effect blocks (standard
#' deviations and correlations, block for block), and its raw threshold
#' vector as stored. The mixing weights are a multinomial logit against
#' the family's reference component `ref` (the last, or the one a
#' formula on every other `theta<k>` leaves out,
#' `mixture_theta_reference()`), and they move as
#' `eta'_k = eta_sigma(k) - eta_sigma(ref)` for each of their
#' coefficients, so that each component keeps its weight. The raw vector
#' moves as it is stored, the review's "raw" swap: between cumulative
#' (first threshold and log gaps) and sratio (the values) that is not the
#' same thresholds, and it still reached the jittered starts' better
#' point on 16 of 18 upward moves of cum_sratio with `(1 | g)` at
#' n = 1000, where moving the values reached 15 (dev/ordstarts-sum4.R).
#' The start is only a start: `mixture_swap_distinct()` drops a start
#' that is the same model, and `mixture_starts()` judges what is kept.
#'
#' An assignment between components whose coefficients, blocks or
#' threshold counts differ cannot carry the parameters over and is
#' skipped. A fit with more than one response runs none.
#'
#' @noRd
mixture_ord_swap_starts <- function(fit, resp) {
  if (length(fit$spec$responses) != 1L) return(list())
  fam <- fit$spec$responses[[resp]]$family
  mx <- fam[["mix"]][["ord"]]
  K <- mx$K
  p <- fit$opt$par
  nm <- outer_par_map(fit)$names
  if (length(nm) != length(p) || K < 2L) return(list())
  ref <- fam[["mix"]][["ref"]] %||% K
  dps <- unique(sub("[0-9]+$", "",
                    setdiff(fam[["dpars"]], paste0("theta", seq_len(K)))))
  # per component: its parameters by dpar and suffix, its blocks and its
  # thresholds
  comp <- lapply(seq_len(K), function(k) {
    idx <- integer(0)
    key <- character(0)
    for (dp in dps) {
      pre <- paste0("^", dp, k, "_")
      i <- grep(pre, nm)
      if (!length(i)) next
      idx <- c(idx, i)
      key <- c(key, paste0(dp, "_", sub(pre, "", nm[i])))
    }
    tn <- mx$tau_names[k]
    tpos <- if (!isTRUE(mx$shared)) which(names(p) == tn) else integer(0)
    blk <- Filter(function(b) {
      d <- b$dpar %||% ""
      grepl("[0-9]+$", d) && sub("[0-9]+$", "", d) %in% dps &&
        as.integer(sub("^[a-z]+", "", d)) == k
    }, fit$frame[["re_blocks"]])
    bsig <- vapply(blk, function(b) {
      paste(sub("[0-9]+$", "", b$dpar), b$group_name,
            paste(b$cnms, collapse = ","), b$covstruct %||% "",
            length(b$theta_idx), length(b$b_idx))
    }, "")
    list(idx = idx[order(key)], key = sort(key), tpos = tpos, blk = blk,
         bsig = bsig)
  })
  # can component k take component s's parameters?
  fits_in <- function(k, s) {
    identical(comp[[k]]$key, comp[[s]]$key) &&
      identical(comp[[k]]$bsig, comp[[s]]$bsig) &&
      length(comp[[k]]$tpos) == length(comp[[s]]$tpos)
  }
  thpos <- which(names(p) == "theta")
  # the weights' coefficients per component, none for the reference
  wt <- lapply(seq_len(K), function(k) {
    if (k == ref) integer(0) else grep(paste0("^theta", k, "_"), nm)
  })
  wkey <- lapply(setdiff(seq_len(K), ref), function(k) {
    sub("^theta[0-9]+_", "", nm[wt[[k]]])
  })
  same_w <- length(wkey) && all(vapply(wkey, identical, NA, wkey[[1L]])) &&
    length(wkey[[1L]]) > 0L
  eta <- function(k) {
    if (k == ref) rep(0, length(wkey[[1L]])) else p[wt[[k]]]
  }
  # the random-effect modes, moved block for block with their sds, as the
  # inner problem's start (mixture_swap_inner()); only where the random
  # vector is the modes `b` alone
  rnd <- fit$obj$env$random
  b0 <- as.numeric(fit$estimates[["b"]])
  move_b <- length(rnd) > 0L && length(b0) == length(rnd)
  out <- list()
  for (sg in mixture_swap_assignments(K)) {
    moved <- which(sg != seq_len(K))
    if (!all(vapply(moved, function(k) fits_in(k, sg[k]), NA))) next
    s <- p
    b <- b0
    for (k in moved) {
      s[comp[[k]]$idx] <- p[comp[[sg[k]]]$idx]
      s[comp[[k]]$tpos] <- p[comp[[sg[k]]]$tpos]
      for (j in seq_along(comp[[k]]$blk)) {
        to <- comp[[k]]$blk[[j]]
        from <- comp[[sg[k]]]$blk[[j]]
        s[thpos[to$theta_idx]] <- p[thpos[from$theta_idx]]
        if (move_b) b[to$b_idx] <- b0[from$b_idx]
      }
    }
    # the weights: each coefficient of eta_k, against the reference
    if (same_w) {
      for (k in setdiff(seq_len(K), ref)) {
        s[wt[[k]]] <- eta(sg[k]) - eta(sg[ref])
      }
    }
    out[[length(out) + 1L]] <- list(
      par = s, sigma = sg, b = if (move_b) b,
      label = paste0(paste0("component ", seq_len(K)[moved], " from ",
                            sg[moved]), collapse = ", "))
  }
  out
}

#' The assignments of K components to try, as permutation vectors, the
#' identity left out. Up to K = 3, all K! - 1 of them; from K = 4, the
#' K (K - 1) / 2 transpositions (two components swapped, the rest in
#' place).
#'
#' On the K = 3 design of dev/ordstarts-k3probe.R (cumulative, sratio
#' and acat, three classes, n = 3000, 20 seeds, dev/ordstarts-sum5.R) a
#' swapped run ended at a better proper point on 11 seeds with all 5
#' assignments and on 9 with the 3 transpositions alone; on 7 seeds the
#' best end came from a 3-cycle, above the transpositions'. So K = 3 runs
#' all 5 (a median 4.75 times the fit, against 7.78 for 10 jittered
#' starts). From K = 4 all assignments are 23 refits and more, more than
#' the jittered starts cost, and no K = 4 design was measured, so the
#' count is held to the transpositions (6 at K = 4).
#'
#' @noRd
mixture_swap_assignments <- function(K) {
  if (K > 3L) {
    pairs <- utils::combn(K, 2L)
    return(lapply(seq_len(ncol(pairs)), function(j) {
      s <- seq_len(K)
      s[pairs[, j]] <- pairs[2:1, j]
      s
    }))
  }
  perm <- function(v) {
    if (length(v) <= 1L) return(list(v))
    out <- list()
    for (i in seq_along(v)) {
      for (r in perm(v[-i])) out[[length(out) + 1L]] <- c(v[i], r)
    }
    out
  }
  Filter(function(s) !identical(s, seq_len(K)), perm(seq_len(K)))
}
