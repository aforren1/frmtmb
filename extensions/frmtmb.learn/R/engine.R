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
#
# The walk yields per-trial FACTORS rather than a running total. That is
# what lets the same recursion answer all three of the protocol's
# likelihood questions, the whole response and one value per subject and
# one value per row, without a second pass and without three definitions
# of the same numbers that could disagree.

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

#' What the trace's summary column is called for this family.
#'
#' @noRd
ln_pcol <- function(spec) if (is.null(spec[["logp"]])) "p" else "dens"

#' Fold a set of option utilities into their log normalizing constant.
#'
#' @noRd
ln_lse <- function(u) {
  out <- u[[1L]]
  for (k in seq_along(u)[-1L]) out <- RTMB::logspace_add(out, u[[k]])
  out
}

#' Lay the recursion out over `nrep` copies of the same design.
#'
#' WHY THIS EXISTS. The importance correction evaluates the whole design
#' once per draw, stacked, so `y` and every distributional parameter
#' arrive `nrep` times as long while the frame block still describes one
#' copy. Row `i` of replicate `k` sits at position `i + (k - 1) * n`, so
#' the block's row NUMBERS shift by that and nothing else does: the walk
#' is the same
#' walk over a longer subject axis. This is the whole of what the
#' protocol's stacking contract costs a family whose loop is already
#' vectorized across subjects.
#'
#' @noRd
ln_stack <- function(block, nrep) {
  idx <- block[["idx"]]
  msk <- block[["mask"]]
  if (nrep == 1L) return(list(idx = idx, mask = msk))
  ns <- nrow(idx)
  take <- rep(seq_len(ns), nrep)
  list(idx = idx[take, , drop = FALSE] +
         rep((seq_len(nrep) - 1L) * block[["n"]], each = ns),
       mask = msk[take, , drop = FALSE])
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
#' @param mode `"terms"` returns the per-trial log-likelihood factors,
#'   one masked vector over subjects per trial, which the three
#'   likelihood slots are each one line of arithmetic away from;
#'   `"trace"` returns a named list of per-row vectors (choice
#'   probabilities, value estimates, prediction errors); `"simulate"`
#'   returns a drawn response vector.
#' @param nrep How many stacked copies of the design `y` and the
#'   distributional parameters carry. Only the likelihood path ever sees
#'   more than one; see `ln_stack()`.
#'
#' @noRd
ln_recurse <- function(block, cd, y, spec,
                       mode = c("terms", "trace", "simulate"),
                       nrep = 1L) {
  mode <- match.arg(mode)
  sk <- ln_stack(block, nrep)
  idx <- sk[["idx"]]
  msk <- sk[["mask"]]
  n <- block[["n"]] * nrep
  ns <- nrow(idx)
  nopt <- spec[["n_option"]]
  nd <- length(nopt)
  # The observed-choice indicators are DATA, so they are built once over
  # the whole response rather than compared at every trial. This is also
  # what keeps the chosen option out of a branch: selecting it is a
  # multiplication by a column of zeros and ones.
  cc <- spec[["choice_col"]]
  ind <- NULL
  if (mode != "simulate") {
    ind <- vector("list", nd)
    for (j in seq_len(nd)) {
      # `choice_col` is the seam a joint choice-and-response-time family
      # needs: its RESPONSE is the response time, so the option taken
      # arrives as an addition term instead of as `y`.
      v <- if (j == 1L) {
        if (is.null(cc)) y else spec[["choice_map"]](cd[[cc]])
      } else {
        cd[[spec[["resp"]][[j]]]]
      }
      ind[[j]] <- lapply(seq_len(nopt[j]), function(k) as.numeric(v == k))
    }
  }
  state <- spec[["init"]](ns, ln_at(cd, idx[, 1L]))
  tm <- if (mode == "terms") vector("list", ncol(idx))
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
    # the response at this trial's rows, which a density rule reads and
    # a softmax rule does not. In `simulate` mode it is the placeholder
    # the caller passed and no rule reads it.
    d[[".y"]] <- y[i]
    ch <- vector("list", nd)
    lp <- 0
    pj <- list()
    for (j in seq_len(nd)) {
      # A DENSITY RULE rather than a softmax over utilities. A family
      # whose response is a response time has no set of option utilities
      # to normalize: its trial contributes one log density, formed from
      # the value store and the response together, and the option taken
      # picks a boundary rather than a term of a sum. Supplying `logp`
      # replaces the choice-and-softmax pair with that; everything else
      # about the walk, the mask and the three likelihood slots is
      # unchanged, which is the whole reason it is a seam in the engine
      # rather than a second engine.
      if (!is.null(spec[["logp"]])) {
        if (mode == "simulate") {
          # the draw writes its own answer back into `d`, in the coding
          # the addition term uses, because only the family knows that
          # coding; the engine reads back the option code it returns
          dr <- spec[["draw"]](state, d)
          d <- dr[["d"]]
          ysim[i[m == 1]] <- dr[["y"]][m == 1]
          ch[[j]] <- lapply(seq_len(nopt[j]),
                            function(k) as.numeric(dr[["choice"]] == k))
          lpj <- 0 * m
        } else {
          ch[[j]] <- lapply(ind[[j]], function(v) v[i])
          lpj <- spec[["logp"]](state, d, ch[[j]])
        }
        lp <- lp + lpj
        next
      }
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
    if (mode == "terms") tm[[t]] <- m * lp
    if (length(sc)) {
      keep <- m == 1
      for (nm in names(sc)) sc[[nm]][i[keep]] <- as.numeric(d[[nm]])[keep]
    }
    upd <- spec[["update"]](state, d, ch)
    if (mode == "trace") {
      # `p` for a softmax family and `dens` for a density family, and
      # the difference is not cosmetic: exp(lp) is a probability in the
      # first case and a density in the second, which may exceed one.
      # Calling a density `p` would invite a reader to check it against
      # a probability.
      rec <- c(state, upd[names(upd) != "state"], pj,
               stats::setNames(list(exp(lp)), ln_pcol(spec)))
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
  switch(mode,
         terms = list(terms = tm, idx = idx, mask = msk),
         trace = tr,
         simulate = list(y = ysim, cols = sc))
}

# ------------------------------------------------------ the three slots
#
# WHY THEY ARE ALL ONE RECURSION. The structured protocol asks a family
# how finely its likelihood factorizes, and this one factorizes twice:
# over subjects, because two subjects share nothing but parameters, and
# again over trials, because
#
#     L(subject) = prod_t P(choice_t | history_t, theta)
#
# is a product, the value store at trial t being a function of that
# subject's data before t and of the parameters. Both statements are
# properties of the recursion rather than of any one rule, so both hold
# for every family in this package. Nothing extra is COMPUTED for them:
# the walk already forms one factor per trial per subject, and the total
# it used to accumulate was the sum of exactly these. Deriving all three
# from `ln_recurse(mode = "terms")` is what keeps them from drifting.

#' The whole response's log-likelihood, as the sum of its pieces.
#'
#' @noRd
ln_loglik <- function(block, cd, y, spec) {
  sum(Reduce(`+`, ln_recurse(block, cd, y, spec, "terms")[["terms"]]))
}

#' One value per subject per replicate, which is the coarsest honest
#' piece and the one `frm(importance =)` resamples.
#'
#' A subject's trials are a product, so its log-likelihood is the
#' elementwise sum of the per-trial vectors. The order is the one the
#' protocol fixes: replicate-major, and within a replicate the subject
#' order of `block[["group"]]`'s levels, which is the order `ln_pack()`
#' laid the rows out in.
#'
#' @noRd
ln_loglik_group <- function(block, cd, y, spec, nrep) {
  Reduce(`+`, ln_recurse(block, cd, y, spec, "terms", nrep)[["terms"]])
}

#' One value per row: the trial's choice probability CONDITIONAL on
#' everything its subject saw before it.
#'
#' These are the factors themselves, scattered back to the rows they
#' came from. ONE sub-assignment, not one per trial: `[<-` on a taped
#' vector copies the whole vector each time it is called, and this
#' package measured that copy costing an order of magnitude at tape
#' build. The padded cells are dropped rather than written, because a
#' pad repeats its subject's FIRST row number and would otherwise
#' overwrite that trial's own factor with a masked zero.
#'
#' @noRd
ln_loglik_row <- function(block, cd, y, spec, nrep, saturated) {
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  tm <- ln_recurse(block, cd, y, spec, "terms", nrep)
  keep <- as.vector(tm[["mask"]]) == 1
  out <- rep(0, length(y))
  out[as.vector(tm[["idx"]])[keep]] <- do.call(c, tm[["terms"]])[keep]
  # off the tape only, and only where the family has one: an advector
  # carries no attributes, and the numeric path is the one residuals()
  # would read. A saturated fit of a nominal choice puts probability one
  # on the option that was taken, so its log-density is zero; a family
  # whose row carries a DENSITY rather than a probability has no such
  # constant and declares none.
  if (saturated && !inherits(out, "advector")) {
    attr(out, "saturated") <- rep(0, length(out))
  }
  out
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
#' @param logp `function(state, d, ch)` returning one log density per
#'   subject, for a family whose response is not the option code.
#'   Supplying it replaces `choice` and the engine's softmax: there are
#'   no utilities to normalize when the trial's contribution is a
#'   first-passage density. `ch[[k]]` is 1 for the subjects that took
#'   option `k`, and `d[[".y"]]` is the response at this trial's rows.
#' @param draw Simulation only, and required with `logp`:
#'   `function(state, d)` returning `list(y = , choice = , d = )`, the
#'   drawn response, the option each subject took as a code in
#'   `1..n_option`, and `d` with the drawn choice written back under the
#'   name its addition term reads.
#' @param choice_col The name in `cd` of the option taken, for a family
#'   whose response is something else.
#' @param choice_map `function(v)` turning that column into option codes
#'   `1..n_option`. `dec()` codes a boundary 0 or 1, and the option
#'   codes are 1 and 2, so the map is what keeps the family's coding out
#'   of the engine.
#'
#' @noRd
ln_spec <- function(n_option, init, choice, update, resp = NULL,
                    draw_env = NULL, sim_cols = character(0),
                    logp = NULL, draw = NULL, choice_col = NULL,
                    choice_map = identity) {
  list(n_option = as.integer(n_option), init = init, choice = choice,
       update = update, resp = resp, draw_env = draw_env,
       sim_cols = sim_cols, logp = logp, draw = draw,
       choice_col = choice_col, choice_map = choice_map)
}
