# Multistart for an hmm() fit.
#
# WHY THIS EXISTS. An HMM likelihood is multimodal, the way a mixture
# likelihood is, and NO convergence diagnostic sees it: on the probe in
# dev/hmm-feasibility.md the default cold start converges 8.1
# log-likelihood units below the global optimum with `convergence == 0`,
# a positive definite Hessian, a gradient of 7e-4, and `diagnose()`
# reporting nothing at all. That is a silent wrong answer produced by
# starting values alone, so the family ships the remedy beside itself.
#
# WHAT IT REPORTS, AND WHY THE SPREAD IS THE SPREAD OF OPTIMA. A
# multistart that reports the spread of its STARTING values says
# nothing; a multistart that quietly drops the refits that did not
# converge reports a spread that is tight because the disagreeing runs
# were removed. So the summary carries three counts (converged, not
# converged, errored) and TWO spreads: over the converged optima, which
# is the number to quote, and over every refit that finished at all,
# which is what says whether dropping any of them mattered.
#
# WHY EACH REFIT RE-EVALUATES THE CALL. Reusing `fit$obj` would be
# cheaper by the tape build, and it would also leave RTMB's
# `last.par.best` at the last refit's parameters. The deferred
# `sdreport()` behind `summary()`, `vcov()` and `confint()` reads
# exactly that, so a `hmm_starts()` call would silently move the
# standard errors of the fit it was handed. A fresh objective per refit
# costs 839 ms at 25 000 rows against 118 s of optimization there
# (dev/scale-findings.md), which is under one percent.
#
# `[[ ]]`, NEVER `$`, on the fit: see the header of hmm.R.

#' The per-parameter perturbation scale, and the label that says where
#' it came from.
#'
#' The fit's own standard errors are the unit a user can reason about:
#' "two standard errors away" is a distance in the answer, not in an
#' arbitrary internal coordinate. `confint()` reports one row per outer
#' parameter, in the order `fit[["opt"]][["par"]]` holds them.
#'
#' WHAT IS AND IS NOT CHECKED HERE, stated exactly, because a first
#' version of this comment claimed protection the code did not give. It
#' said the row order was CHECKED against the `est` column;
#' `confint.frmtmb_fit()` builds that column as `est <- object$opt$par`,
#' so the test compared a vector with itself and its measured difference
#' was exactly 0. It could never have failed.
#'
#' What is checked now:
#'
#' \itemize{
#'   \item the ROW COUNT, which can differ: a `confint()` that returned
#'     a different number of rows than the optimizer has parameters
#'     cannot be read by position at all;
#'   \item the standard errors themselves, against `vcov()`, on the
#'     rows the two have in common. `vcov()` reaches the same covariance
#'     by its own path and names its rows itself, so a mismatch there is
#'     a real disagreement about which standard error belongs to which
#'     coefficient. This is the check that can fire.
#' }
#'
#' What is NOT checked, and what stands in for it: that row `i` of
#' `confint()` is parameter `i` of `opt$par`. Nothing public names the
#' optimizer's positions, so the guarantee for the fixed-effect rows is
#' the `vcov()` cross-check above by NAME, and for the variance rows
#' (`theta_1 ...`), which have no second path, it is the row count plus
#' `confint()`'s documented ordering. A permutation confined to those
#' rows would not be caught, and the cost of it would be jittering one
#' log standard deviation by another's standard error.
#'
#' @noRd
hmm_starts_scale <- function(fit) {
  par <- fit[["opt"]][["par"]]
  fallback <- list(s = rep(1, length(par)), how = "unit")
  ci <- tryCatch(suppressWarnings(stats::confint(fit)),
                 error = function(e) NULL)
  if (is.null(ci) || nrow(ci) != length(par) ||
      !all(c("lwr", "upr") %in% colnames(ci))) {
    return(fallback)
  }
  s <- (as.numeric(ci[, "upr"]) - as.numeric(ci[, "lwr"])) /
    (2 * stats::qnorm(0.975))
  nm <- rownames(ci)
  v <- tryCatch(suppressWarnings(sqrt(diag(stats::vcov(fit)))),
                error = function(e) NULL)
  if (is.null(v) || !length(v) || is.null(names(v))) return(fallback)
  shared <- intersect(names(v), nm %||% character(0))
  if (!length(shared)) return(fallback)
  a <- s[match(shared, nm)]
  b <- as.numeric(v[shared])
  cmp <- is.finite(a) & is.finite(b) & b > 0
  if (!any(cmp)) return(fallback)
  if (max(abs(a[cmp] - b[cmp]) / b[cmp]) > 1e-6) return(fallback)
  ok <- is.finite(s) & s > 0
  if (!any(ok)) return(fallback)
  # a parameter whose standard error did not come back is moved by the
  # typical distance of the ones that did, rather than left alone
  s[!ok] <- stats::median(s[ok])
  list(s = s, how = if (all(ok)) "se" else "se-partial")
}

#' How much better one log-likelihood has to be than another before
#' this function will call it a DIFFERENT optimum.
#'
#' There are three places that ask that question: whether a refit
#' displaces `$best`, whether the printed summary says the original
#' found a local optimum, and how `hmm_starts_modes()` merges. They used
#' to answer it three ways, and two of the three carried no tolerance at
#' all, so on six of six unimodal fits the summary announced a local
#' optimum on gaps of 6e-11 to 1.5e-09 while the modes table beside it
#' correctly reported ONE optimum. The two thresholds differed by a
#' factor of 1.7e05. The fix is not a better constant, it is having one
#' definition: a check that fires on a correct model is worse than no
#' check, and three copies of a threshold is how one of them comes to.
#'
#' The scale is `grad_tol^2`, relative to the log-likelihood's own
#' magnitude. Squaring is deliberate: a refit is declared CONVERGED when
#' its relative gradient is at most `grad_tol`, and near a smooth
#' optimum the objective is second order in the distance the gradient
#' measures, so `grad_tol^2` is the log-likelihood difference that
#' convergence test can no longer resolve.
#'
#' ONE FUNCTION IS NOT ONE THRESHOLD. A first version of this fix had
#' all three sites calling this helper and passing DIFFERENT references:
#' the displacement passed the incumbent's log-likelihood and the mode
#' merge passed each cluster head. For a negative log-likelihood the
#' better optimum always has the smaller magnitude, so the merge
#' threshold was always the smaller of the two and a band of `grad_tol`
#' always existed where the modes table reported two optima while the
#' verdict three lines above said the incumbent was the best found. At
#' `grad_tol = 0.086` on the probe in `dev/hmm-feasibility.md` that band
#' is real and the printout contradicts itself. So the threshold is
#' computed ONCE per call, from the original fit's log-likelihood, and
#' carried on the returned object as `$mode_tol`; the print method reads
#' it rather than recomputing it, which is what makes the two
#' impossible to separate again.
#'
#' @noRd
hmm_starts_tol <- function(grad_tol, ref) {
  grad_tol^2 * max(abs(ref), 1)
}

#' Objective and gradient evaluations an optimizer took, as a length-2
#' numeric.
#'
#' A wall-clock second is a measurement of the machine as much as of the
#' model, and this project has had two timing claims evaporate under
#' replication. An evaluation count does not move with the load, so the
#' cost of a refit is reported both ways. Only nlminb reports the pair;
#' `optim` counts differently and a custom optimizer need not count at
#' all, so an unrecognized shape is NA rather than a guess.
#'
#' @noRd
hmm_starts_evals <- function(fit) {
  ev <- fit[["opt"]][["evaluations"]]
  if (is.numeric(ev) && length(ev) == 2L) return(as.numeric(ev))
  cn <- fit[["opt"]][["counts"]]
  if (is.numeric(cn) && length(cn) == 2L) return(as.numeric(cn))
  c(NA_real_, NA_real_)
}

#' Split a flat outer vector back into the `start =` list shape.
#'
#' `fit[["estimates"]]` already has the right names and includes the
#' conditional modes; only the outer components are replaced, so the
#' random effects start from where the original fit left them.
#'
#' @noRd
hmm_starts_relist <- function(fit, v) {
  st <- fit[["estimates"]]
  pn <- names(fit[["opt"]][["par"]])
  for (cp in unique(pn)) {
    idx <- which(pn == cp)
    tgt <- st[[cp]]
    if (is.null(tgt) || length(tgt) != length(idx)) {
      # a component the optimizer only partly owns (a fixed betad cell,
      # say) cannot be written back by position; leave it at the fit's
      # own values rather than guess
      next
    }
    st[[cp]][] <- v[idx]
  }
  st
}

#' Refits of an `hmm()` model from jittered starting values
#'
#' A hidden Markov likelihood has more than one mode, and no convergence
#' diagnostic reports which one a fit reached. `hmm_starts()` refits the
#' model `n` times from starting values scattered around its estimates,
#' returns the best fit it found, and reports the spread of the OPTIMA
#' the refits reached, so that "the optimizer converged" and "the
#' optimizer converged to the answer" stop being the same sentence.
#'
#' @section What the numbers mean:
#' \describe{
#'   \item{the spread}{The range of the log-likelihood over the refits
#'     that CONVERGED. A second spread is printed over every refit that
#'     finished, converged or not, because a summary that silently drops
#'     the runs that disagreed reports a spread that is tight for the
#'     wrong reason.}
#'   \item{converged}{A refit counts as converged when it returns
#'     without error and its largest absolute gradient, divided by the
#'     absolute log-likelihood it reached, is at most `grad_tol`. The
#'     ratio rather than the gradient itself is what is tested, because
#'     on a long chain the objective's own magnitude sets the scale the
#'     optimizer stops at.}
#'   \item{the modes}{Refits whose log-likelihoods agree to within
#'     `grad_tol^2` times the original fit's log-likelihood are counted
#'     as the same optimum. Two modes with a wide gap between them is
#'     the finding; one mode found `n` times is the reassurance.}
#' }
#'
#' That last threshold is returned as `$mode_tol`, and the same value
#' decides all three questions this function answers: whether a refit
#' displaces `$best`, whether the summary says the original found a
#' local optimum, and how the modes are merged. **It moves as the SQUARE
#' of `grad_tol`**, so loosening the convergence test by a factor of 100
#' loosens the merge threshold by 10 000. That is deliberate and it has
#' a range: on the probe `hmm_starts()` exists for, the 8.099-unit
#' detection survives every `grad_tol` up to about 0.08 and is gone at
#' 0.086, where the threshold has grown past the gap itself. The default
#' of 1e-3 puts the threshold at about 1e-03 log-likelihood units on a
#' fit of -1096, which is seven orders of magnitude below that gap.
#'
#' @section What it costs, and what `n` to use:
#' A refit is a whole fit: the frame, the tape and the optimization
#' again. One call therefore costs `n` fits plus one `sdreport()`, which
#' is paid once for the jitter scale and is the same one `summary()`
#' would have paid.
#'
#' Measured at the 25 000-row design of `dev/scale-findings.md` (50
#' sequences of 500, `K = 3` gaussian, `tr12 ~ (1 | id)`), alone on the
#' machine, three blocks of the same work from the same seed: one refit
#' is **175.7 s and 26.5 gradient evaluations**, against 183.5 s and 16
#' gradient evaluations for the plain fit it started from, and the
#' `sdreport()` is 37.2 s. So `n = 4` is 12.3 minutes there, `n = 8` is
#' 24.0 and `n = 16` is 47.5. The instrument's control (the largest
#' block over the smallest) reported 1.023 and the gradient counts were
#' identical across blocks, which is the check a clock cannot make.
#' `dev/latent-2p3-starts-cost.R` is the script.
#'
#' Twenty-four minutes is not something to do by accident on a model
#' that took three to fit, so on anything that size set `n` from the
#' `seconds` line of a first, small call rather than from the default.
#'
#' `n` cannot have a default that is right at every size, and the one
#' here is a starting point rather than a recommendation. What it is
#' chosen from: on the probe in `dev/hmm-feasibility.md` a single
#' jittered refit at `jitter = 2` reaches the better optimum about 6
#' times in 10, so `n = 8` misses it with probability about 4e-4 and
#' `n = 4` with probability about 2e-2. That rate is a property of that
#' surface, not a constant, so treat `n` as something to raise when the
#' printed summary shows more than one optimum and to price from the
#' `seconds` line when the model is large.
#'
#' @section Where the starts come from:
#' Each start is the fit's own estimates plus a normal draw per outer
#' parameter, with standard deviation `jitter` times that parameter's
#' standard error. Standard errors make `jitter` a distance in the
#' answer rather than in an internal coordinate. Where the fit reports
#' no usable standard error, the median of the ones it does report is
#' used, and where it reports none at all every parameter is moved by
#' `jitter` on the internal scale; the printed summary names which of
#' the three applied. Random effects start from the original fit's
#' conditional modes, since the inner problem is re-solved either way.
#'
#' @param fit A `frmtmb_fit` with an [hmm()] family.
#' @param n Number of refits. Each one costs a whole fit.
#' @param jitter Perturbation size, in standard errors of the fit's own
#'   estimates. Larger values explore further and converge less often.
#' @param seed Optional integer. Given, the starts are reproducible and
#'   the session's random stream is restored afterwards.
#' @param grad_tol Largest `max|gradient| / |logLik|` a refit may have
#'   and still count as converged.
#' @param keep If `TRUE`, every refit that returned a fit at all is
#'   kept in `$fits`, in the order they ran, so it is shorter than `n`
#'   when a refit errored. These are whole fitted models with their
#'   tapes; `n` of them at a realistic size is gigabytes, which is why
#'   the default discards all but the best.
#' @return An object of class `frmtmb_hmm_starts`:
#'   `$best` (the highest-likelihood fit, which is the original when no
#'   refit beat it, and whose own `call` carries the winning start by
#'   value so that re-evaluating it reproduces the fit), `$table` (one
#'   row per start: `logLik`, `grad_rel`, `status`, `seconds`, and the
#'   `fn_evals` / `gr_evals` the optimizer took, which is a cost measure
#'   that does not move with the machine's load), `$spread` and
#'   `$spread_all` (`NA` when only one fit is in the set, since a
#'   comparison that was never made is not a spread of zero), `$modes`,
#'   `$mode_tol` (the one threshold above), and the counts.
#' @seealso [hmm()], [frmtmb::frm_allfit()], which varies the OPTIMIZER
#'   from one fixed start and answers a different question.
#' @examples
#' set.seed(4)
#' n_seq <- 12; len <- 20
#' G <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2, byrow = TRUE)
#' mu <- c(0, 3); sg <- c(0.6, 0.6)
#' dd <- do.call(rbind, lapply(seq_len(n_seq), function(id) {
#'   s <- integer(len); s[1] <- 1L
#'   for (t in 2:len) s[t] <- sample.int(2, 1, prob = G[s[t - 1], ])
#'   data.frame(id = id, t = seq_len(len), y = rnorm(len, mu[s], sg[s]))
#' }))
#' fit <- frm(bf(y ~ 1),
#'            family = hmm(K = 2, gaussian(), time = t, group = id),
#'            data = dd)
#' ms <- hmm_starts(fit, n = 3, seed = 1)
#' ms
#' logLik(ms$best)
#' @export
hmm_starts <- function(fit, n = 8L, jitter = 2, seed = NULL,
                       grad_tol = 1e-3, keep = FALSE) {
  hmm_rspec(fit, "hmm_starts()")
  check_flag(keep, "keep")
  if (length(n) != 1L || !is.numeric(n) || is.na(n) || n != round(n) ||
      n < 1) {
    stop("hmm_starts(): `n` must be a whole number of refits, at ",
         "least 1, not ", arg_desc(n), call. = FALSE)
  }
  n <- as.integer(n)
  if (length(jitter) != 1L || !is.numeric(jitter) || !is.finite(jitter) ||
      jitter <= 0) {
    stop("hmm_starts(): `jitter` must be one positive number, the ",
         "perturbation size in standard errors, not ", arg_desc(jitter),
         call. = FALSE)
  }
  if (length(grad_tol) != 1L || !is.numeric(grad_tol) ||
      !is.finite(grad_tol) || grad_tol <= 0) {
    stop("hmm_starts(): `grad_tol` must be one positive number, the ",
         "largest max|gradient| / |logLik| a refit may have and still ",
         "count as converged, not ", arg_desc(grad_tol), call. = FALSE)
  }
  par <- fit[["opt"]][["par"]]
  if (!length(par)) {
    stop("hmm_starts(): this fit has no free outer parameters, so ",
         "there is nothing to restart from", call. = FALSE)
  }

  cl <- fit[["call"]]
  if (is.null(cl)) {
    stop("hmm_starts(): this fit carries no call to re-evaluate, so it ",
         "cannot be refit", call. = FALSE)
  }
  env <- new.env(parent = parent.frame())
  # a `[[` READ of an absent name on a call is an error rather than
  # NULL, which is why this asks names(cl) instead
  if (length(fit[["data2"]]) && !("data2" %in% names(cl))) {
    assign(".hmm_starts_data2", fit[["data2"]], envir = env)
    cl[["data2"]] <- as.name(".hmm_starts_data2")
  }
  # Assemble the frame once before paying for any refit. The call names
  # its data by SYMBOL, so a fit built somewhere the caller cannot see
  # fails here, once and by name, instead of n times as n identical
  # errors counted as n bad starts.
  probe <- cl
  probe[["dry_run"]] <- "frame"
  # dropping an argument a call does not have is another out-of-bounds
  # `[[`, not a no-op
  if ("start" %in% names(probe)) probe[["start"]] <- NULL
  ok <- tryCatch({
    suppressMessages(suppressWarnings(eval(probe, env)))
    TRUE
  }, error = function(e) conditionMessage(e))
  if (!isTRUE(ok)) {
    stop("hmm_starts(): the model's own call could not be re-evaluated ",
         "from here, so there is nothing to refit. The data and any ",
         "other symbol the call names must be visible from where ",
         "hmm_starts() is called. The call reported: ", ok, call. = FALSE)
  }

  sc <- hmm_starts_scale(fit)
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
              add = TRUE)
    } else {
      on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())),
              add = TRUE)
    }
    set.seed(seed)
  }

  ll0 <- as.numeric(stats::logLik(fit))
  # the ONE threshold, from the ONE reference, for the whole call
  mode_tol <- hmm_starts_tol(grad_tol, ll0)
  g0 <- max(abs(drop(fit[["obj"]][["gr"]](par))))
  ev0 <- hmm_starts_evals(fit)
  rows <- data.frame(
    start = 0L, logLik = ll0, grad_rel = g0 / max(abs(ll0), 1),
    status = "original", seconds = NA_real_,
    fn_evals = ev0[1L], gr_evals = ev0[2L],
    stringsAsFactors = FALSE)
  best <- fit
  best_st <- NULL
  kept <- list()

  for (i in seq_len(n)) {
    v <- as.numeric(par) + stats::rnorm(length(par), 0, jitter * sc$s)
    st <- hmm_starts_relist(fit, v)
    cli <- cl
    assign(".hmm_starts_start", st, envir = env)
    cli[["start"]] <- as.name(".hmm_starts_start")
    cli[["se"]] <- FALSE
    t0 <- Sys.time()
    f <- tryCatch(suppressWarnings(suppressMessages(eval(cli, env))),
                  error = function(e) e)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(f, "error")) {
      rows <- rbind(rows, data.frame(
        start = i, logLik = NA_real_, grad_rel = NA_real_,
        status = "error", seconds = secs,
        fn_evals = NA_real_, gr_evals = NA_real_,
        stringsAsFactors = FALSE))
      next
    }
    lli <- as.numeric(stats::logLik(f))
    gi <- max(abs(drop(f[["obj"]][["gr"]](f[["opt"]][["par"]]))))
    rel <- gi / max(abs(lli), 1)
    conv <- is.finite(lli) && is.finite(rel) && rel <= grad_tol
    evi <- hmm_starts_evals(f)
    rows <- rbind(rows, data.frame(
      start = i, logLik = lli, grad_rel = rel,
      status = if (conv) "converged" else "not converged",
      seconds = secs, fn_evals = evi[1L], gr_evals = evi[2L],
      stringsAsFactors = FALSE))
    if (keep) kept[[length(kept) + 1L]] <- f
    bl <- as.numeric(stats::logLik(best))
    if (conv && is.finite(lli) && lli > bl + mode_tol) {
      best <- f
      best_st <- st
    }
  }

  if (!is.null(best_st)) {
    # the returned fit's own call has to reproduce the returned fit, so
    # the winning start goes in BY VALUE rather than as a symbol that
    # will not exist a moment from now
    bc <- fit[["call"]]
    bc[["start"]] <- best_st
    best[["call"]] <- bc
  }

  # The ORIGINAL is in both spreads whatever its own gradient looks
  # like: it is the point the spread is measured FROM, and dropping it
  # would leave the reader with a spread of refits and no anchor. Its
  # own verdict is reported separately instead, because a summary that
  # holds the refits to a test it does not hold the incumbent to is
  # exactly the guard that fails open.
  original_converged <- is.finite(rows$grad_rel[1L]) &&
    rows$grad_rel[1L] <= grad_tol
  conv_ll <- rows$logLik[rows$status %in% c("original", "converged")]
  fin_ll <- rows$logLik[rows$status != "error"]
  modes <- hmm_starts_modes(conv_ll, mode_tol)
  # A spread over ONE value is not a spread of zero, it is a comparison
  # that was never made, and the two must not print the same. Every
  # refit erroring used to give `spread = 0` and `spread_all = 0`, which
  # a caller reading the fields without the summary would read as
  # perfect agreement.
  spread_of <- function(v) if (length(v) > 1L) diff(range(v)) else NA_real_
  structure(
    list(best = best,
         table = rows,
         spread = spread_of(conv_ll),
         spread_all = spread_of(fin_ll),
         modes = modes,
         n = n, jitter = jitter, seed = seed, grad_tol = grad_tol,
         mode_tol = mode_tol,
         scale_how = sc$how,
         n_converged = sum(rows$status == "converged"),
         n_not_converged = sum(rows$status == "not converged"),
         n_error = sum(rows$status == "error"),
         original_logLik = ll0,
         original_converged = original_converged,
         fits = if (keep) kept else NULL),
    class = "frmtmb_hmm_starts")
}

#' Distinct optima among a set of log-likelihoods, largest first, at
#' the ABSOLUTE threshold the caller computed once. It takes the value
#' and not the ingredients on purpose: a threshold recomputed here from
#' a local reference is how the two came apart before.
#'
#' @noRd
hmm_starts_modes <- function(ll, mode_tol) {
  ll <- ll[is.finite(ll)]
  if (!length(ll)) {
    return(data.frame(logLik = numeric(0), n = integer(0)))
  }
  o <- sort(ll, decreasing = TRUE)
  lab <- integer(length(o))
  lab[1L] <- 1L
  for (i in seq_along(o)[-1L]) {
    ref <- o[which(lab == lab[i - 1L])[1L]]
    same <- abs(o[i] - ref) <= mode_tol
    lab[i] <- lab[i - 1L] + as.integer(!same)
  }
  data.frame(logLik = tapply(o, lab, function(z) z[1L]),
             n = as.integer(table(lab)),
             row.names = NULL)
}

#' @export
print.frmtmb_hmm_starts <- function(x, digits = 10, ...) {
  cat("<hmm_starts> ", x[["n"]], " refits, jitter ", x[["jitter"]],
      " (", switch(x[["scale_how"]],
                   se = "standard errors",
                   `se-partial` = "standard errors, some imputed",
                   "internal scale, no standard errors available"),
      ")\n", sep = "")
  cat("  converged      ", x[["n_converged"]], " of ", x[["n"]], "\n",
      sep = "")
  if (x[["n_not_converged"]]) {
    cat("  not converged  ", x[["n_not_converged"]],
        " (max|grad| / |logLik| above ", format(x[["grad_tol"]]),
        "), and OUT of the converged spread below\n", sep = "")
  }
  if (x[["n_error"]]) {
    cat("  errored        ", x[["n_error"]],
        ", reaching no optimum at all, and out of BOTH spreads\n",
        sep = "")
  }
  bll <- as.numeric(stats::logLik(x[["best"]]))
  cat("\n  best logLik    ", format(bll, digits = digits), "\n",
      sep = "")
  cat("  original       ", format(x[["original_logLik"]],
                                  digits = digits), sep = "")
  gap <- bll - x[["original_logLik"]]
  # the SAME value the modes table and the `$best` test used, read off
  # the object rather than recomputed, so this sentence cannot
  # contradict the count three lines below it
  if (gap > x[["mode_tol"]]) {
    cat("   (", format(gap, digits = 4),
        " BELOW the best: the original fit found a local optimum)\n",
        sep = "")
  } else {
    cat("   (the original fit was the best found)\n")
  }
  if (!isTRUE(x[["original_converged"]])) {
    cat("\n  the ORIGINAL fit does not meet that test either",
        " (max|grad| / |logLik| = ",
        format(x[["table"]]$grad_rel[1L], digits = 3),
        ").\n  It stays in both spreads below, as the point they are",
        " measured from.\n", sep = "")
  }
  nc <- sum(x[["table"]]$status %in% c("original", "converged"))
  nf <- sum(x[["table"]]$status != "error")
  fmt_spread <- function(v, n) {
    if (n < 2L) "not measurable, only one fit is in this set" else
      format(v, digits = digits)
  }
  cat("\n  logLik spread, the original and the ", nc - 1L,
      " converged refits: ", fmt_spread(x[["spread"]], nc),
      "\n", sep = "")
  cat("  logLik spread, every one of the ", nf,
      " that finished    : ", fmt_spread(x[["spread_all"]], nf),
      "\n", sep = "")
  m <- x[["modes"]]
  # "among the converged" would be wrong when the incumbent itself
  # failed the test: it is in this set either way, and the line above
  # has already said so
  cat("\n  distinct optima among the original and the converged: ",
      nrow(m), "\n", sep = "")
  if (nrow(m)) {
    print(data.frame(logLik = format(m$logLik, digits = digits),
                     found = m$n), row.names = FALSE)
  }
  tt <- x[["table"]]$seconds
  gg <- x[["table"]]$gr_evals[-1L]
  cat("\n  seconds: ", format(sum(tt, na.rm = TRUE), digits = 4),
      " total, ", format(stats::median(tt, na.rm = TRUE), digits = 4),
      " median per refit\n", sep = "")
  if (any(is.finite(gg))) {
    cat("  gradient evaluations: ",
        format(sum(gg, na.rm = TRUE), digits = 6), " total, ",
        format(stats::median(gg, na.rm = TRUE), digits = 6),
        " median per refit\n", sep = "")
  }
  invisible(x)
}
