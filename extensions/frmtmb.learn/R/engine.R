# The recursion engine.
#
# WHY ONE ENGINE. Every family in this package is the same walk. Carry a
# value store with one entry per subject, visit the trials in order, at
# each trial turn the store into a choice probability and turn the
# outcome into a new store. What differs between a two-armed delta rule
# and a Kalman filter over four restless arms is the body of those two
# steps, not the walk. So the walk is written once, here, and a family
# supplies a CHOICE RULE and a LEARNING RULE.
#
# WHY IT IS SHAPED THIS WAY. Two constraints, both measured rather than
# assumed (dev/learn-findings.md).
#
# The loop runs over TRIALS and vectorizes over SUBJECTS. Written the
# other way round, one iteration per row with the store updated by
# sub-assignment, RTMB pays for it at tape construction: the `[<-`
# replacement copies the vector it writes into, so the penalty grows
# with the row count.
#
# The same code runs on the tape and off it. `RTMB::logspace_add()` and
# `sign()` both accept plain doubles as well as advectors, so the
# likelihood, the per-trial fitted values and the forward simulator are
# ONE recursion at three settings of `mode`, not three copies that have
# to be kept in step. The worked example this package generalizes,
# `inst/rl/rw-delta.R` in frmtmb, needed two copies.

#' Gather every per-row quantity at one trial's rows.
#'
#' A length-one entry is left alone rather than indexed: a
#' distributional parameter with no formula arrives as one value, and
#' broadcasting it to `n` first would add `n` tape nodes and then read
#' `n_subj` of them back at every trial, for a quantity that is one
#' number. R's own recycling does the rest.
#'
#' @noRd
ln_at <- function(cd, i) {
  for (k in seq_along(cd)) {
    v <- cd[[k]]
    if (length(v) > 1L) cd[[k]] <- v[i]
  }
  cd
}

#' Advance a value store past one trial, with padded cells frozen.
#'
#' The blend is what makes unequal trial counts free of branches: on a
#' padded cell `m` is 0 and the store comes back untouched, whatever the
#' learning rule did with it. Folding the mask into the learning rate
#' instead, as the worked example does, is one operation cheaper and has
#' to be got right again in every rule; the engine's job is that a new
#' family cannot get it wrong.
#'
#' MATCHED BY NAME, and the name matters. The first version matched by
#' POSITION, which is correct only while a rule happens to rebuild its
#' store in the order `init()` declared it. The Kalman filter does not:
#' it fills one loop over the four arms, so it returns
#' `mu1, s1, mu2, s2, ...` where `init()` returned `mu1..mu4, s1..s4`,
#' and every posterior mean was blended with a posterior variance. The
#' fit converged, recovered plausible parameters, and was wrong by 206
#' log-likelihood units. That is why test-reference.R exists.
#'
#' @noRd
ln_blend <- function(old, new, m) {
  if (!setequal(names(old), names(new))) {
    stop("a learning rule returned the value stores ",
         paste(sort(names(new)), collapse = ", "), " where init() ",
         "declared ", paste(sort(names(old)), collapse = ", "),
         ". Every slot must come back under the name it went in with",
         call. = FALSE)
  }
  for (k in names(old)) old[[k]] <- old[[k]] + m * (new[[k]] - old[[k]])
  old
}

#' One draw per row from a matrix of category probabilities.
#'
#' Written as a loop over the (few) options rather than through
#' `apply(cumsum)`, because `apply()` drops a one-row matrix to a vector
#' and the failure would be silent.
#'
#' @noRd
ln_rcat <- function(pm) {
  ns <- nrow(pm)
  u <- stats::runif(ns)
  acc <- numeric(ns)
  out <- rep(ncol(pm), ns)
  done <- rep(FALSE, ns)
  for (k in seq_len(ncol(pm))) {
    acc <- acc + pm[, k]
    hit <- !done & (u <= acc)
    out[hit] <- k
    done <- done | hit
  }
  out
}

#' The larger of two AD quantities, without a branch.
#'
#' RTMB refuses a comparison on an advector outright and has no
#' `CondExp`, so a model-based value that takes the best stage-two
#' option cannot ask which one it is. `abs()` is defined on an advector
#' and its kink is the function's own, not an artifact: the true maximum
#' is not differentiable where the two arguments cross either.
#'
#' @noRd
ln_max2 <- function(a, b) 0.5 * (a + b + abs(a - b))

#' Fold a set of option utilities into their log normalizing constant.
#'
#' @noRd
ln_lse <- function(u) {
  out <- u[[1L]]
  for (k in seq_along(u)[-1L]) out <- RTMB::logspace_add(out, u[[k]])
  out
}

#' Walk a value-learning recursion over every subject at once.
#'
#' @param block The frame block built by `ln_block()`: `idx` and `mask`
#'   are subject-by-trial matrices of row numbers, and of 1 for a real
#'   trial and 0 for padding.
#' @param cd Named list of per-row quantities, each of length `n` or 1:
#'   the distributional parameters on their natural scales and the
#'   addition-term columns.
#' @param y The response, coded 1 to K.
#' @param spec The model, as built by `ln_spec()`.
#' @param mode `"loglik"` returns one AD scalar, `"trace"` returns a
#'   named list of per-row vectors (choice probabilities, value
#'   estimates, prediction errors), `"simulate"` returns a drawn
#'   response vector.
#'
#' @noRd
ln_recurse <- function(block, cd, y, spec,
                       mode = c("loglik", "trace", "simulate")) {
  mode <- match.arg(mode)
  idx <- block[["idx"]]
  msk <- block[["mask"]]
  n <- block[["n"]]
  ns <- nrow(idx)
  nopt <- spec[["n_option"]]
  nd <- length(nopt)
  # The observed-choice indicators are DATA, so they are built once over
  # the whole response rather than compared at every trial. This is also
  # what keeps the chosen option out of a branch: selecting it is a
  # multiplication by a column of zeros and ones.
  ind <- NULL
  if (mode != "simulate") {
    ind <- vector("list", nd)
    for (j in seq_len(nd)) {
      v <- if (j == 1L) y else cd[[spec[["resp"]][[j]]]]
      ind[[j]] <- lapply(seq_len(nopt[j]), function(k) as.numeric(v == k))
    }
  }
  state <- spec[["init"]](ns, ln_at(cd, idx[, 1L]))
  ll <- 0
  tr <- list()
  ysim <- if (mode == "simulate") rep(NA_real_, n)
  # A task whose trial holds more than one choice, or whose outcome is
  # produced by the environment rather than fixed in advance, cannot
  # state a draw as one response vector. Such a family declares the
  # other columns its draw writes, and the caller rebuilds a whole
  # dataset from them.
  sc <- if (mode == "simulate" && length(spec[["sim_cols"]])) {
    stats::setNames(rep(list(rep(NA_real_, n)), length(spec[["sim_cols"]])),
                    spec[["sim_cols"]])
  } else list()
  for (t in seq_len(ncol(idx))) {
    i <- idx[, t]
    m <- msk[, t]
    d <- ln_at(cd, i)
    d[["m"]] <- m
    ch <- vector("list", nd)
    lp <- 0
    pj <- list()
    for (j in seq_len(nd)) {
      u <- spec[["choice"]](state, d, j)
      lse <- ln_lse(u)
      if (mode == "simulate") {
        pm <- matrix(0, ns, nopt[j])
        for (k in seq_len(nopt[j])) {
          pm[, k] <- exp(as.numeric(u[[k]]) - as.numeric(lse))
        }
        drw <- ln_rcat(pm)
        ch[[j]] <- lapply(seq_len(nopt[j]), function(k) as.numeric(drw == k))
        if (j == 1L) {
          ysim[i[m == 1]] <- drw[m == 1]
        } else {
          d[[spec[["resp"]][[j]]]] <- drw
        }
        if (!is.null(spec[["draw_env"]])) d <- spec[["draw_env"]](d, ch, j)
      } else {
        ch[[j]] <- lapply(ind[[j]], function(v) v[i])
      }
      uc <- ch[[j]][[1L]] * u[[1L]]
      for (k in seq_along(u)[-1L]) uc <- uc + ch[[j]][[k]] * u[[k]]
      lpj <- uc - lse
      lp <- lp + lpj
      if (mode == "trace" && nd > 1L) pj[[paste0("p", j)]] <- exp(lpj)
    }
    if (mode == "loglik") ll <- ll + sum(m * lp)
    if (length(sc)) {
      keep <- m == 1
      for (nm in names(sc)) sc[[nm]][i[keep]] <- as.numeric(d[[nm]])[keep]
    }
    upd <- spec[["update"]](state, d, ch)
    if (mode == "trace") {
      rec <- c(state, upd[names(upd) != "state"], pj, list(p = exp(lp)))
      if (!length(tr)) tr <- lapply(rec, function(z) rep(NA_real_, n))
      keep <- m == 1
      for (nm in names(rec)) {
        z <- as.numeric(rec[[nm]])
        if (length(z) == 1L) z <- rep(z, ns)
        tr[[nm]][i[keep]] <- z[keep]
      }
    }
    state <- ln_blend(state, upd[["state"]], m)
  }
  switch(mode, loglik = ll, trace = tr,
         simulate = list(y = ysim, cols = sc))
}

#' The model half of a learning family.
#'
#' `init`, `choice` and `update` are the three tape-safe functions a
#' family supplies. Everything else about a family (its parameters,
#' their links, its addition terms, what it refuses) is declared to
#' frmtmb rather than computed.
#'
#' @param n_option Options per decision, one entry per decision within a
#'   trial. `2` is a two-armed choice; `c(2, 2)` is the two-step task,
#'   whose trial holds two choices.
#' @param init `function(ns, d1)` returning a named list of length-`ns`
#'   value stores, given the number of subjects and the first trial's
#'   per-row data. Reading the first trial is what lets an initial value
#'   be a fitted parameter.
#' @param choice `function(state, d, j)` returning `n_option[j]` utility
#'   vectors over subjects. The engine softmaxes them: a utility is a
#'   log unnormalized probability, so the inverse temperature
#'   multiplies inside this function.
#' @param update `function(state, d, ch)` returning `list(state = )`
#'   plus any per-subject quantity worth recording, such as a prediction
#'   error. `ch[[j]][[k]]` is 1 for the subjects that took option `k` at
#'   decision `j`.
#' @param resp For a trial with more than one decision, the name in `cd`
#'   of each later decision's observed choice.
#' @param draw_env Simulation only: `function(d, ch, j)` returning `d`
#'   with the environment's answer to the choice just drawn filled in.
#' @param sim_cols Names in `d` that a draw rewrites and a caller must
#'   read back to rebuild a whole simulated dataset. Empty for a task
#'   whose payoff schedule is fixed in advance, where the response
#'   vector is the entire draw.
#'
#' @noRd
ln_spec <- function(n_option, init, choice, update, resp = NULL,
                    draw_env = NULL, sim_cols = character(0)) {
  list(n_option = as.integer(n_option), init = init, choice = choice,
       update = update, resp = resp, draw_env = draw_env,
       sim_cols = sim_cols)
}
