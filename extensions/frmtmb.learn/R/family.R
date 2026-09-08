# Everything every family in this package shares: the frame block, the
# refusals, the family object, and the three post-fit slots. A family
# file is then the model and nothing else.

#' The frame block: from a long data frame to subject by trial.
#'
#' Runs once, at frame assembly. Everything it returns is plain data,
#' which is what the structured protocol requires and what lets
#' `refit()` rebuild it.
#'
#' `idx` holds row NUMBERS rather than values, so the recursion can
#' index any per-row quantity with it: the response, the payoffs and
#' every linear predictor. Subjects with fewer trials than the longest
#' are padded on the right with a repeat of their own first row, and
#' `mask` is 0 there.
#'
#' @noRd
ln_block <- function(resp, spec, av, mf, y, n) {
  fam <- resp[["family"]]
  nm <- fam[["family"]]
  lrn <- fam[["learn"]]
  gv <- eval(lrn[["subject_expr"]], mf, resp[["formula_env"]])
  if (anyNA(gv)) {
    stop(nm, "(): every row needs a subject, and ",
         deparse1(lrn[["subject_expr"]]), " has ", sum(is.na(gv)),
         " missing value(s)", call. = FALSE)
  }
  gv <- factor(gv)
  tv <- if (is.null(lrn[["trial_expr"]])) {
    # no trial column: the data frame's own row order within a subject
    out <- integer(n)
    for (g in split(seq_len(n), gv)) out[g] <- seq_along(g)
    out
  } else {
    v <- eval(lrn[["trial_expr"]], mf, resp[["formula_env"]])
    if (anyNA(v)) {
      stop(nm, "(): the trial variable '", deparse1(lrn[["trial_expr"]]),
           "' has missing values, so the order of the recursion is ",
           "undefined at those rows", call. = FALSE)
    }
    as.numeric(v)
  }
  ln_pack(gv, tv, n, nm)
}

#' Lay a subject column and a trial column out as the recursion walks
#' them.
#'
#' Shared by the fitted path, which reaches it through `ln_block()`, and
#' by the generative path in `frm_task_simulate()`, which has no frame.
#' One layout means a simulated dataset and the fit of it are indexed
#' the same way by construction rather than by agreement.
#'
#' @noRd
ln_pack <- function(gv, tv, n, nm) {
  rows <- lapply(split(seq_len(n), gv), function(r) r[order(tv[r])])
  key <- paste(as.integer(gv), tv, sep = "|")
  if (anyDuplicated(key)) {
    dup <- key[duplicated(key)][1L]
    stop(nm, "(): trial numbers must be unique within a subject; ",
         sum(key == dup), " rows share one. A learning rule updates ",
         "once per trial, so two rows at one trial have no order to ",
         "learn in", call. = FALSE)
  }
  len <- lengths(rows)
  nt <- max(len)
  # The matrix() before each transpose is load-bearing: vapply() drops
  # to a plain vector when nt == 1, t() of a vector is 1-by-n_subj, and
  # the block would come back with every subject read as one more trial
  # of a single learner, carrying values across subject boundaries.
  idx <- t(matrix(vapply(rows, function(r) c(r, rep(r[1L], nt - length(r))),
                         integer(nt)), nrow = nt))
  mask <- t(matrix(vapply(len, function(l) as.numeric(seq_len(nt) <= l),
                          numeric(nt)), nrow = nt))
  # `group` is the protocol's reserved name for the family's own
  # independent unit, and a subject IS that unit here: two subjects
  # share parameters and nothing else. The core reads it to align its
  # own grouping against this one before it admits the importance
  # correction, and to order what `loglik_group` returns. `subject` is
  # the same vector under this package's own name, which
  # frm_value_trace() reads; they are one column with two readers rather
  # than two columns that could drift, because both are `gv`.
  list(idx = idx, mask = mask, len = len, n_subj = length(rows),
       n_trial = nt, group = gv, subject = gv, trial = tv,
       levels = levels(gv), n = n)
}

#' The refusal a non-rowwise likelihood owes the rowwise contract.
#'
#' One template for every family in the package, filled in at run time:
#' the sentence is the same fact about all of them.
#'
#' @noRd
ln_no_rowwise <- function(nm) {
  stop("A ", nm, "() trial's probability depends on every earlier trial ",
       "of the same subject, so the family has no row-wise log-density. ",
       "Use logLik() for the total, or fitted() for the per-trial ",
       "probability of the choice that was made", call. = FALSE)
}

#' The response coding every family in this package uses.
#'
#' @noRd
ln_valid_y <- function(y, k, nm) {
  if (!all(y %in% seq_len(k))) {
    stop(nm, "(): the response is the option chosen on each trial, ",
         "coded 1 to ", k, ". Two-level codings that start at 0 are the ",
         "usual mistake; add one", call. = FALSE)
  }
}

#' Refusals that follow from the likelihood not factorizing over rows.
#'
#' @noRd
ln_check_spec <- function(nm, resp, spec, av) {
  if (length(spec[["responses"]]) > 1L || isTRUE(spec[["rescor"]])) {
    stop(nm, "() supports univariate models only: the recursion is a ",
         "likelihood over one response's trial sequences", call. = FALSE)
  }
  bad <- intersect(c("weights", "cens", "trunc_lb", "trunc_ub", "se"),
                   names(av))
  if (length(bad)) {
    stop(nm, "() cannot be combined with ", bad[1L], "(): that term ",
         "reshapes a per-row likelihood contribution, and a trial's ",
         "contribution here is conditional on every earlier trial of ",
         "the same subject", call. = FALSE)
  }
  invisible(NULL)
}

#' The four sentences every family in this package owes, in its own name.
#'
#' Three of them are consequences of the recursion reading a whole trial
#' history, and none of those is a switch that could be turned on: a
#' synthetic covariate grid has no history to replay, and the tape holds
#' no registered observation vector to differentiate. The fourth,
#' deviance, is a consequence of the response being NOMINAL rather than
#' of the likelihood: the magnitude is available now that the family
#' declares `loglik_row`, and only the sign is not.
#'
#' @noRd
ln_refusals <- function(nm, nominal = TRUE) {
  list(
    newdata_response = paste0(
      "predict(type = 'response') on a ", nm, "() fit is not available ",
      "for newdata: a trial's choice probability is conditional on the ",
      "outcomes and choices of every earlier trial of the same subject, ",
      "and newdata carries no block to replay. On new data, ask for ",
      "predict(type = 'link', dpar = ) for a learning parameter"),
    conditional_effects = paste0(
      "conditional_effects() is not available for a ", nm, "() fit: the ",
      "expected response is a choice probability that depends on a whole ",
      "trial history, which the synthetic grid this function builds does ",
      "not have. Plot a learning parameter itself with predict(dpar = )"),
    osa = paste0(
      "residuals(type = 'osa') is not available for a ", nm, "() fit: ",
      "the tape holds the recursion over each whole sequence with no ",
      "registered observation vector, and the response is a nominal ",
      "option code rather than a quantity the tape could differentiate. ",
      "frm_value_trace() gives the per-trial choice probability"),
    deviance = if (nominal) {
      paste0(
        "residuals(type = 'deviance') is not available for a ", nm,
        "() fit, and what is missing is the SIGN rather than the ",
        "magnitude. The family declares loglik_row(), so each trial's ",
        "log-likelihood and its saturated value are both available and ",
        "their difference is the unit deviance; a deviance RESIDUAL ",
        "then needs the direction of the row's departure from its ",
        "conditional mean, and the response is the option a subject ",
        "took, coded 1 to K, which is nominal and has no mean to ",
        "depart from. frm_value_trace() gives the per-trial log ",
        "likelihood as log(p)")
    } else {
      paste0(
        "residuals(type = 'deviance') is not available for a ", nm,
        "() fit. The family declares loglik_row(), so each trial's own ",
        "log-likelihood is available; what a unit deviance needs ",
        "beside it is the SATURATED value, and a row here contributes ",
        "a density rather than a probability. The supremum of a Wiener ",
        "log density at a fixed response time is unbounded over the ",
        "parameters, so there is no constant to compare against. ",
        "frm_value_trace() gives the per-trial log density as ",
        "log(dens)")
    })
}

#' The quantities the post-fit slots read, through exported accessors.
#'
#' @noRd
ln_parts <- function(fit, nm) {
  rspec <- single_response(fit, paste0("a ", nm, "() fit"))
  r <- rspec[["resp_name"]]
  list(dpars = lapply(eval_dpars(fit)[[r]], as.numeric),
       aterms = fit[["frame"]][["aterm_values"]][[r]],
       y = as.numeric(fit[["frame"]][["y"]][[r]]))
}

#' Replay the fitted recursion in plain numbers.
#'
#' @noRd
ln_trace_at <- function(fit, block, spec, nm) {
  pp <- ln_parts(fit, nm)
  ln_recurse(block, c(pp[["dpars"]], pp[["aterms"]]), pp[["y"]], spec,
             "trace")
}

#' The structure: this package's half of the protocol.
#'
#' Every capability flag starts FALSE for a likelihood that does not
#' factorize over rows, so the refusals are explanations rather than
#' switches.
#'
#' @noRd
ln_structure <- function(nm, spec, sim = TRUE, refusals = list(),
                         saturated = TRUE) {
  # A family whose choice rule is a density has a response that is not
  # an option code, and three of its refusals then say something else.
  # Read off the rule rather than declared again beside it.
  nominal <- is.null(spec[["logp"]])
  frmtmb_structure(
    frame_vars = function(fam) {
      list(fam[["learn"]][["subject_expr"]], fam[["learn"]][["trial_expr"]])
    },
    # An unanswered trial produces no prediction error, so dropping its
    # row IS the right recursion. This is where a learning rule differs
    # from a hidden Markov chain, which keeps the row because the state
    # still moves; see frmtmb.latent::hmm().
    keep_na = FALSE,
    check_spec = function(resp, spec_, av) ln_check_spec(nm, resp, spec_, av),
    frame_block = ln_block,
    # `weights` arrives at all THREE likelihood slots and all three
    # ignore it. The protocol's convention is that the family applies
    # the row weights exactly once, and ignoring them is only correct
    # because they cannot be anything but 1 here. Two routes could
    # supply them and both are shut: check_spec refuses `weights()` by
    # name, and cluster weights reach the objective only through
    # sandwich.R, which gates on the `cluster_robust` support flag that
    # this structure leaves FALSE. A trial's factor could not be
    # reweighted on its own in any case, because its value depends on
    # every earlier trial of the same subject.
    loglik = function(y, dpars, aterms, weights, block, extra) {
      ln_loglik(block, c(dpars, aterms), y, spec)
    },
    # THE TWO FACTORIZATION SLOTS, and each declares something that is
    # true of this recursion rather than something convenient.
    #
    # loglik_group: a subject. Two subjects share the parameters and
    # nothing else, so the response's likelihood is a product over them.
    #
    # loglik_row: a trial. Given the parameters, a subject's likelihood
    # is prod_t P(choice_t | history_t), and history_t is that subject's
    # own earlier rows. So a row HAS a conditional log-density, it
    # depends only on its own subject's random effects, and the row
    # values sum to the group values which sum to the total. All three
    # come off one call to the recursion (R/engine.R), which is what
    # stops them disagreeing.
    #
    # `unit` is a different question and keeps its own answer. These say
    # how finely the likelihood factorizes; `unit` says what may honestly
    # be left OUT, and dropping one trial changes every later trial's
    # value store, so the leave-one-out unit is still the whole sequence.
    loglik_row = function(y, dpars, aterms, weights, block, extra) {
      ln_loglik_row(block, c(dpars, aterms), y, spec,
                    NROW(y) %/% block[["n"]], saturated)
    },
    loglik_group = function(y, dpars, aterms, weights, block, extra) {
      ln_loglik_group(block, c(dpars, aterms), y, spec,
                      NROW(y) %/% block[["n"]])
    },
    unit = "one subject's trial sequence",
    # NO fitted_mean, and it is a decision rather than an omission.
    #
    # Core forms a residual as `y - fitted_mean(fit, block)` and a
    # fitted value as `fitted_mean` itself. The response here is the
    # option a subject took, coded 1 to K, and it is NOMINAL: arm 2 is
    # not twice arm 1 and the four decks of the Iowa gambling task have
    # no order at all. So there is no number `mu` for which `y - mu` is
    # a residual, and any `fitted_mean` this family supplied would make
    # fitted(), predict(type = "response") and residuals() return
    # arithmetic on a category code. core::cox() and
    # frmtmb.spline::royston_parmar() decline to invent a mean for the
    # same kind of reason.
    #
    # What replaces them is frm_value_trace(), which returns MORE than
    # a mean could: the value estimates each choice was made on, the
    # prediction error, and the fitted probability of the choice that
    # was made. The identity sum(log(trace$p)) == logLik(fit) is
    # asserted in test-engine.R.
    #
    # The cost is real and is stated here rather than discovered: the
    # rw_delta example in frmtmb DOES have fitted(), because it codes a
    # two-armed choice 0 and 1 and returns P(arm 1), so `y - p` is the
    # ordinary binary residual. That trick does not survive a fourth
    # arm, and one engine over two and four options is worth more here
    # than fitted() on the two-option half of it.
    sim_ctx = if (sim) {
      function(ctx) {
        blk <- ctx[["block"]]
        ln_recurse(blk,
                   c(lapply(ctx[["dpars"]], as.numeric), ctx[["aterms"]]),
                   rep(1, blk[["n"]]), spec, "simulate")[["y"]]
      }
    },
    refusals = utils::modifyList(ln_refusals(nm, nominal), refusals))
}

#' Assemble a learning family.
#'
#' The one place a family object is built, so that a family file holds
#' the model and nothing else.
#'
#' @noRd
ln_family <- function(nm, subject_expr, trial_expr, dpars, links, primary,
                      inits, aterms, spec, sim = TRUE, data_map = NULL,
                      constants = list(), refusals = list(),
                      sim_refusal = NULL, saturated = TRUE,
                      valid_y = NULL, finalize = NULL) {
  if (is.null(subject_expr)) {
    stop(nm, "(subject =) names the column that separates one learner's ",
         "trial sequence from the next", call. = FALSE)
  }
  k <- spec[["n_option"]][[1L]]
  fam <- frmtmb_family(
    nm,
    dpars = dpars,
    links = links,
    primary_dpars = primary,
    type = "discrete",
    required_aterms = aterms,
    # The likelihood factorizes over trials but is not ROWWISE: a
    # trial's factor reads the whole history of its subject, which this
    # signature cannot see. Returning something anyway is the silent lie
    # the structured protocol exists to remove.
    lpdf = function(y, dpars, aterms, extra = NULL) ln_no_rowwise(nm),
    # The six softmax families all take the same response, an option
    # code 1 to K, and share one validator. A family whose response is
    # something else brings its own.
    valid_y = if (is.null(valid_y)) {
      function(y, aterms) ln_valid_y(y, k, nm)
    } else {
      valid_y
    },
    family_finalize = finalize,
    init_dpars = inits,
    structure = ln_structure(nm, spec, sim, refusals, saturated))
  # Core reads this slot when a family has no simulator, and appends it
  # to whichever entry point refused. It is not a `supports` flag,
  # because having a simulator is not a capability the protocol tracks.
  if (!is.null(sim_refusal)) fam[["sim_refusal"]] <- sim_refusal
  fam[["learn"]] <- c(list(subject_expr = subject_expr,
                           trial_expr = trial_expr, spec = spec,
                           dpars = dpars, data_map = data_map),
                      constants)
  fam
}
