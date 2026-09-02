# Hidden Markov models: hmm(K, family) as a first-class family.
#
# WHAT THIS IS. An HMM makes the observation at time t come from one of K
# state-dependent distributions, and the state itself follows a Markov
# chain. The state sequence is DISCRETE, so the Laplace approximation
# neither applies to it nor is needed: the forward algorithm marginalizes
# it EXACTLY, as a per-sequence recursion of K x K products. That
# recursion tapes in RTMB, which is what makes this a family rather than
# a separate fitting engine (dev/hmm-feasibility.md, probe A).
#
# WHAT IS BORROWED FROM mixture(). Everything except the recursion:
# suffixed per-state dpars carrying the full formula grammar, the
# multinomial-logit block (K copies of it, one per source state), the
# logspace_add fold, quantile-spread initialization, and the
# sum-inside-the-integral objective shape (R/objective.R). Random effects
# need no special handling at all - they live in the state dpars' linear
# predictors, the forward algorithm is exact GIVEN them, and the Laplace
# integrates them outside. That composition is the v0.19 latent-class
# insight one rung up.
#
# WHY THE RECURSION IS SLICED BY TIME, NOT BY SEQUENCE. The obvious
# spelling loops over sequences and, within one, over time, extracting a
# K-vector of emission log-densities per step. That is O(n K^2) SCALAR
# extractions on the tape. Sorting the sequences by decreasing length
# makes the set of sequences still running at step s a PREFIX of that
# order, so one step can be taken for every live sequence at once: the
# alpha vector for state k is a vector over sequences, and a step is
# K^2 vectorized logspace_add calls of that length. The tape then holds
# O(Tmax K^2) nodes instead of O(n K^2), which is what keeps the tape
# build usable on many short sequences (probe A2: rows are the unit of
# cost, and 1000 sequences of 5 cost what 1 of 5000 does).
#
# CONVENTIONS, PINNED AGAINST depmixS4 (probe B1, every coefficient to
# four decimals):
#   - each ROW of the transition matrix is its own multinomial logit,
#     with STATE 1 as the reference cell in every row (depmixS4's
#     convention; hmmTMB instead references the diagonal, which spans the
#     same model with different coefficients);
#   - the covariate value at time t drives the transition from t to
#     t + 1, so the last row of a sequence never contributes one.
#
# `[[ ]]`, NEVER `$`, ON A FAMILY OR FRAME FIELD. `$` falls back to
# partial matching when no exact name matches, so a gate spelled
# `is.null(fam$hmm)` would silently start reading a future `fam$hmm_opts`
# and fire on families that are not HMMs at all. The same holds inside
# the sequence structure below, which carries both `m` and `mask`/`miss`
# and both `n` and `n_seq`: those reads are bracketed for the same
# reason, so that removing a field breaks loudly instead of quietly.

## ---- helpers ---------------------------------------------------------

#' The dpar name of transition cell `i -> j`. Column 1 is the reference
#' cell of row `i`, so only `j >= 2` has a name.
#'
#' @noRd
hmm_tr_name <- function(i, j) paste0("tr", i, j)

#' Every transition dpar name of a K-state chain, row-major.
#'
#' @noRd
hmm_tr_names <- function(K) {
  as.vector(t(outer(seq_len(K), seq_len(K)[-1L], hmm_tr_name)))
}

#' Whether one parsed dpar is an intercept (or a fixed constant): no
#' design columns beyond the intercept, no random effects, no smooth,
#' Gaussian-process, or spatial terms. Read off the SPEC, because
#' `init = "stationary"` has to be decided before the designs are built.
#'
#' @noRd
hmm_dpar_is_constant <- function(dp) {
  if (!is.null(dp$constant)) return(TRUE)
  if (!is.null(dp$nl_body)) return(FALSE)
  lists_empty <- vapply(
    c("re", "smooth", "mo", "miterms", "csterms", "gpterms",
      "carterms", "spdeterms"),
    function(nm) !length(dp[[nm]] %||% list()), TRUE
  )
  if (!all(lists_empty)) return(FALSE)
  f <- dp$fixed
  is.null(f) || identical(deparse1(f[[length(f)]]), "1")
}

#' A `time =` / `group =` argument as an expression. Both a bare symbol
#' (the `ar()` / `cosy()` spelling) and a one-sided formula (the
#' `mixture(groups = ~g)` spelling) are accepted, so neither habit is a
#' trap.
#'
#' @noRd
hmm_capture <- function(e, what) {
  if (is.null(e)) return(NULL)
  if (inherits(e, "formula")) e <- e[[length(e)]]
  if (is.call(e) && identical(e[[1L]], as.name("~"))) {
    return(e[[length(e)]])
  }
  if (!is.name(e) && !is.call(e)) {
    stop("hmm(): `", what, "` must name a variable in the data, as a ",
         "bare name (", what, " = subject) or a one-sided formula (",
         what, " = ~subject)", call. = FALSE)
  }
  e
}

#' Starting intercepts for state `k`'s location parameter, on the
#' RESPONSE scale.
#'
#' `mixture()`'s convention - the k-th of K + 1 evenly spaced response
#' quantiles - is the right default and the one probe F1 measured: a
#' start ON the label-symmetry axis (every state mean at the median) is a
#' fixed point of the objective's `K!` symmetry, and the fit stalls at
#' the one-state solution 262 log-likelihood units below the optimum.
#' Two cases need a fallback: a quantile outside the link's range (a zero
#' count under a log link), and a matrix response, which has no quantile
#' at all. Both spread the states around the pooled fit on the LINK
#' scale instead, which breaks the same symmetry.
#'
#' @noRd
hmm_mu_init <- function(k, K, link) {
  force(k); force(K); force(link)
  function(y, aterms) {
    spread <- function(centre) {
      link$linkinv(centre + seq(-1, 1, length.out = K)[k])
    }
    if (is.matrix(y) || !is.numeric(y)) return(spread(0))
    yv <- as.numeric(y)
    qs <- stats::quantile(yv, seq_len(K) / (K + 1), names = FALSE,
                          na.rm = TRUE)
    if (all(is.finite(link$linkfun(qs)))) return(qs[k])
    m <- mean(yv, na.rm = TRUE)
    centre <- link$linkfun(m)
    if (!is.finite(centre)) centre <- 0
    spread(centre)
  }
}

#' Starting transition logits: a sticky chain. Row `i`'s reference cell
#' is state 1, so staying is the +1.5 logit when `i > 1` and every free
#' logit is -1.5 when `i == 1`.
#'
#' @noRd
hmm_tr_init <- function(i, j) {
  v <- if (i == j) 1.5 else -1.5
  force(v)
  function(y, aterms) v
}

## ---- the family constructor -----------------------------------------

#' Hidden Markov models
#'
#' `hmm(K, family)` fits a `K`-state hidden Markov model: the response at
#' each time point is drawn from one of `K` state-dependent copies of
#' `family`, and the unobserved state follows a first-order Markov chain
#' along `time` within `group`. The state sequence is summed out exactly
#' by the forward algorithm, which is evaluated on the same AD tape as
#' everything else, so random effects, smooths and distributional
#' predictors compose with it unchanged.
#'
#' @section Parameters:
#' Each of the wrapped family's distributional parameters is copied once
#' per state and suffixed with the state index, exactly as [mixture()]
#' does: `hmm(2, gaussian())` has `mu1`, `mu2`, `sigma1`, `sigma2`. The
#' main model formula applies to every state's location parameter;
#' override one state with `bf(y ~ x, mu2 ~ x + (1 | id))`. Every dpar
#' takes the full formula grammar, random effects included.
#'
#' Transition probabilities are a row-wise multinomial logit with
#' **state 1 as the reference cell in every row**, named `tr{i}{j}` for
#' the move from state `i` to state `j` (`j >= 2`): a two-state chain has
#' `tr12` and `tr22`. This is `depmixS4`'s parameterization, and the
#' coefficients agree with it to four decimals (`dev/hmm-feasibility.md`,
#' probe B1). `hmmTMB` and `moveHMM` instead reference the diagonal;
#' the two span the same model and only the reported coefficients differ.
#' `trans` gives every transition cell the same default formula; a single
#' cell is overridden the ordinary way, `bf(y ~ 1, tr12 ~ x)`. The
#' covariate value at time `t` drives the transition from `t` to
#' `t + 1`, so a sequence's last row never contributes a transition -
#' again `depmixS4`'s convention.
#'
#' @section Initial distribution:
#' \describe{
#'   \item{`"stationary"`}{The stationary distribution of the transition
#'     matrix, solved on the tape. Costs no parameters and is the right
#'     default for long sequences and for many short ones. It needs a
#'     CONSTANT transition matrix, so it is refused when any transition
#'     dpar carries a predictor.}
#'   \item{`"estimated"`}{`K - 1` free logits (state 1 the reference),
#'     reported as `hmm_ldel_1 ...` and counted in `df`.}
#'   \item{`"uniform"`}{Fixed at `1 / K`.}
#' }
#'
#' @section Decoding and the fitted values:
#' `E[y_t]` under an HMM is `sum_k P(S_t = k | y) mu_k(x_t)`, which needs
#' the smoothed state probability and therefore a backward pass over the
#' whole sequence. [hmm_probs()] returns those probabilities and
#' [hmm_viterbi()] the maximum-a-posteriori state path; [fitted()],
#' `predict(type = "response")` and [residuals()] all route through
#' [hmm_probs()], so they report the occupancy-weighted mean rather than
#' any single state's.
#'
#' @section Label switching and local optima:
#' The likelihood is invariant to permuting the states, so the state
#' means start at spread response quantiles (see [mixture()]); a start
#' with every state mean equal sits on the symmetry axis and the
#' optimizer never leaves it. Relabeling between runs is expected and is
#' not fought. Multimodality is real and is not signalled by any
#' convergence diagnostic: on a random-effect model the default cold
#' start has been measured converging 8.1 log-likelihood units below the
#' global optimum with `convergence == 0` and a positive-definite
#' Hessian. Compare several starts ([frm_allfit()]) before reporting.
#'
#' @section Missing responses:
#' An `NA` response is a time point the chain passes THROUGH without
#' emitting: the row is kept, its emission factor is masked out, and the
#' transition into and out of it is unchanged. This departs from the
#' package-wide `na.action`, deliberately - dropping the row would
#' shorten the chain and make one transition stand in for several, which
#' biases the transition matrix (measurably: with 3 of 20 points missing
#' the fitted off-diagonal moved from 0.115 to 0.128 against a true
#' 0.10). `nobs()` therefore counts every row. `fitted()` and
#' `residuals()` are `NA` at a masked row, while [hmm_probs()] is not:
#' the neighbouring observations still say where the chain was.
#'
#' @section Cost:
#' Evaluating the likelihood is linear in the number of ROWS and free in
#' the number of sequences: a thousand sequences of five cost what one of
#' five thousand does. Both the tape and the decoding passes are
#' organized by time step, so what they actually cost is
#' `length of the LONGEST sequence` times `K^2`. At 20 000 rows a fit
#' takes about 3 s cut into a thousand short sequences and about the same
#' as one long chain, but [hmm_probs()] takes 0.06 s in the first shape
#' and 1.3 s in the second. Below 5 000 rows nothing is noticeable.
#' On a long chain the objective's own magnitude also makes the
#' optimizer's relative convergence test bite before its gradient test
#' does; judge `max|grad|` per observation rather than in absolute terms.
#'
#' @section Boundaries:
#' An unpenalized multinomial logit will send an emission or transition
#' probability to 0 whenever a category is rare inside a state, and the
#' optimizer then reports singular convergence at a perfectly good
#' optimum. [set_prior()] on the affected logit is the remedy.
#'
#' @section Random effects and the Laplace approximation:
#' A random effect in a state's linear predictor is integrated by the
#' Laplace approximation OUTSIDE the exact state sum. That integrand is a
#' mixture over state sequences and is not Gaussian, so the approximation
#' is genuinely approximate here even for a gaussian response. Measured
#' against adaptive Gauss-Hermite quadrature on a 40-sequence,
#' 25-observation model with one scalar random intercept, the bias was
#' **-0.126 in the log-likelihood** (8.9e-5 relative) and 4.4e-4 absolute
#' in the parameters (`dev/hmm-feasibility.md`, probe D1). `quadrature =
#' TRUE` is refused: its rule integrates a random effect against a
#' PRODUCT of per-row densities, and the forward algorithm is not one.
#'
#' @section What is refused:
#' `REML` (the restricted likelihood would integrate out only one state's
#' location coefficients, which matches no standard definition),
#' `quadrature`, `frmtmb_control(profile = TRUE)`, `weights()`, `cens()`,
#' `trunc()`, `se()` and `mi()` on the response, multivariate models and
#' `rescor`, `residuals(type = "osa")`, `predict(se.fit = TRUE)` on the
#' response scale, and `conditional_effects()`. A grouping in which every
#' sequence has length 1 is refused too: the chain is then unidentified
#' and the model is a [mixture()].
#'
#' @param K Number of hidden states (at least 2, at most 9 - beyond that
#'   the `tr{i}{j}` dpar names stop being unambiguous).
#' @param family The state-dependent (emission) family. Any univariate
#'   family with a `mu` parameter works, plus [multinomial()] for
#'   categorical emissions. Ordinal families, mixtures and nested `hmm()`
#'   components are refused.
#' @param time Variable giving the order of observations within a
#'   sequence, as a bare name or a one-sided formula. Omit it to use each
#'   row's position in the data.
#' @param group Variable identifying the sequences, as a bare name or a
#'   one-sided formula. Omit it to treat the whole data set as one
#'   sequence.
#' @param init Initial-state distribution: `"stationary"` (the default),
#'   `"estimated"`, or `"uniform"`.
#' @param trans Default one-sided formula for every transition cell.
#' @return A `frmtmb_family`.
#' @seealso [hmm_probs()], [hmm_viterbi()], [mixture()]
#' @examples
#' set.seed(11)
#' n_seq <- 20; len <- 25
#' G <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2, byrow = TRUE)
#' mu <- c(0, 3); sg <- c(0.6, 0.6)
#' dd <- do.call(rbind, lapply(seq_len(n_seq), function(id) {
#'   s <- integer(len); s[1] <- 1L
#'   for (t in 2:len) s[t] <- sample.int(2, 1, prob = G[s[t - 1], ])
#'   data.frame(id = id, t = seq_len(len),
#'              y = rnorm(len, mu[s], sg[s]), state = s)
#' }))
#' fit <- frm(bf(y ~ 1),
#'            family = hmm(K = 2, gaussian(), time = t, group = id),
#'            data = dd)
#' fixef(fit)
#'
#' # smoothed state probabilities and the MAP path
#' head(hmm_probs(fit))
#' mean(hmm_viterbi(fit) == dd$state)
#'
#' # fitted() is the occupancy-weighted mean, not state 1's
#' cor(fitted(fit), dd$y)
#'
#' \donttest{
#' # one state's mean takes its own predictor, random effects included
#' frm(bf(y ~ 1, mu2 ~ 1 + (1 | id)),
#'     family = hmm(K = 2, gaussian(), time = t, group = id),
#'     data = dd)
#'
#' # covariate-dependent transitions: trans = sets every cell's default
#' dd$x <- rnorm(nrow(dd))
#' frm(bf(y ~ 1),
#'     family = hmm(K = 2, gaussian(), time = t, group = id,
#'                  init = "estimated", trans = ~x),
#'     data = dd)
#' }
#' @export
hmm <- function(K, family = stats::gaussian(), time = NULL, group = NULL,
                init = c("stationary", "estimated", "uniform"),
                trans = ~1) {
  time_expr <- hmm_capture(substitute(time), "time")
  group_expr <- hmm_capture(substitute(group), "group")
  init <- match.arg(init)
  if (length(K) != 1L || !is.numeric(K) || is.na(K) ||
      K != round(K) || K < 2) {
    stop("hmm() needs the number of hidden states as a whole number ",
         "of at least 2: hmm(2, gaussian())", call. = FALSE)
  }
  K <- as.integer(K)
  if (K > 9L) {
    stop("hmm(): at most 9 states are supported, because the transition ",
         "dpar names concatenate the state indices (tr12, tr22, ...) and ",
         "stop being unambiguous at two digits; K = ", K, " was given",
         call. = FALSE)
  }
  if (!inherits(trans, "formula") || length(trans) != 2L) {
    stop("hmm(): `trans` must be a one-sided formula giving the default ",
         "predictor of every transition cell, for example trans = ~x",
         call. = FALSE)
  }
  comp <- as_frmtmb_family(family)
  if (!is.null(comp[["hmm"]])) {
    stop("hmm(): the state-dependent family cannot itself be an hmm(); ",
         "higher-order and hierarchical chains are a different model",
         call. = FALSE)
  }
  if (!is.null(comp[["mix"]])) {
    stop("hmm(): the state-dependent family cannot be a mixture(); a ",
         "mixture inside a state is not identified against the state ",
         "itself", call. = FALSE)
  }
  if (identical(comp$type, "ordinal") || isTRUE(comp$drop_intercept)) {
    stop("hmm(): ordinal families are not supported as state-dependent ",
         "distributions ('", comp$family, "'); their thresholds are ",
         "family-level extra parameters with no per-state copy",
         call. = FALSE)
  }
  if (!is.null(comp$extra_pars)) {
    stop("hmm(): the state-dependent family '", comp$family, "' carries ",
         "family-level extra parameters, which have no per-state copy",
         call. = FALSE)
  }
  primaries <- comp$primary_dpars %||% "mu"
  if (!"mu" %in% comp$dpars && !identical(comp$family, "multinomial")) {
    stop("hmm(): the state-dependent family needs a 'mu' parameter; '",
         comp$family, "' has ", paste(comp$dpars, collapse = ", "),
         call. = FALSE)
  }

  dpars <- character(0)
  links <- list()
  init_fns <- list()
  for (k in seq_len(K)) {
    for (dp in comp$dpars) {
      nm <- paste0(dp, k)
      dpars <- c(dpars, nm)
      links[[nm]] <- comp$links[[dp]]
      init_fns[[nm]] <- if (dp %in% primaries) {
        hmm_mu_init(k, K, comp$links[[dp]])
      } else {
        comp$init_dpars[[dp]]
      }
    }
  }
  for (i in seq_len(K)) {
    for (j in seq_len(K)[-1L]) {
      nm <- hmm_tr_name(i, j)
      dpars <- c(dpars, nm)
      links[[nm]] <- "identity"
      init_fns[[nm]] <- hmm_tr_init(i, j)
    }
  }
  init_fns <- Filter(Negate(is.null), init_fns)

  state_dpars <- function(dp_all, k) {
    stats::setNames(
      lapply(comp$dpars, function(dp) dp_all[[paste0(dp, k)]]),
      comp$dpars
    )
  }
  extra_pars <- if (init == "estimated") {
    function(y, aterms) list(hmm_ldel = rep(0, K - 1L))
  }

  fam <- frmtmb_family(
    paste0("hmm(", K, ", ", comp$family, ")"),
    dpars = dpars,
    links = links,
    # The likelihood does NOT factorize over rows, so there is no
    # per-row log-density to hand back. Returning one (state 1's, say)
    # is exactly the silent lie rung 2 exists to remove, so the slot
    # refuses instead; the objective takes its own branch and never
    # calls this.
    lpdf = function(y, dpars, aterms, extra = NULL) {
      stop("The likelihood of an hmm() family is a per-SEQUENCE forward ",
           "recursion, not a product of per-row densities, so it has no ",
           "row-wise log-density. Use logLik() for the total, or ",
           "hmm_probs() for the state probabilities", call. = FALSE)
    },
    valid_y = comp$valid_y,
    init_dpars = init_fns,
    type = comp$type,
    primary_dpars = as.vector(t(outer(primaries, seq_len(K), paste0))),
    extra_pars = extra_pars
  )
  fam[["hmm"]] <- list(
    K = K,
    init = init,
    comp = comp,
    time_expr = time_expr,
    group_expr = group_expr,
    trans = trans,
    tr_names = hmm_tr_names(K),
    state_dpars = state_dpars,
    state_lpdf = function(y, dp_all, aterms, k) {
      comp$lpdf(y, state_dpars(dp_all, k), aterms)
    },
    state_mean = function(dp_all, aterms, k) {
      response_mean(comp, state_dpars(dp_all, k), aterms)
    },
    state_sim = function(dp_all, aterms, n, k, rows) {
      if (is.null(comp$sim)) {
        stop("simulate(): the state-dependent family '", comp$family,
             "' of this hmm() fit has no simulator yet", call. = FALSE)
      }
      dk <- lapply(state_dpars(dp_all, k), function(v) {
        rep(v, length.out = n)[rows]
      })
      comp$sim(dk, aterms, length(rows))
    }
  )
  # The transition cells share one default predictor, which the parser
  # copies in wherever the user wrote no formula of their own.
  if (!identical(deparse1(trans[[2L]]), "1")) {
    fam[["default_forms"]] <- stats::setNames(
      rep(list(trans), length(fam[["hmm"]][["tr_names"]])),
      fam[["hmm"]][["tr_names"]]
    )
  }
  fam
}

#' Whether a spec holds an `hmm()` response.
#'
#' @noRd
has_hmm <- function(spec) {
  any(vapply(spec$responses,
             function(r) !is.null(r$family[["hmm"]]), TRUE))
}

## ---- frame: the sequence structure ----------------------------------

#' Time-ordered row indices per sequence, sorted so the sequences still
#' running at step `s` are a PREFIX of the order.
#'
#' `slice[[s]]` holds the row index of the `s`-th observation of every
#' sequence that has one, and `rslice[[s]]` the `s`-th FROM THE END,
#' which is what the backward pass walks. `m[s]` is how many sequences
#' those are. See the file header for why the recursion is organized this
#' way instead of one sequence at a time.
#'
#' @noRd
hmm_seq_structure <- function(gidx, tidx, n) {
  by_g <- split(seq_len(n), gidx)
  by_g <- lapply(by_g, function(r) r[order(tidx[r])])
  len <- lengths(by_g)
  ord <- order(-len, seq_along(len))
  by_g <- by_g[ord]
  len <- len[ord]
  tmax <- max(len)
  m <- vapply(seq_len(tmax), function(s) sum(len >= s), integer(1))
  slice <- lapply(seq_len(tmax), function(s) {
    vapply(by_g[seq_len(m[s])], function(r) r[s], integer(1),
           USE.NAMES = FALSE)
  })
  rslice <- lapply(seq_len(tmax), function(s) {
    vapply(by_g[seq_len(m[s])], function(r) r[length(r) - s + 1L],
           integer(1), USE.NAMES = FALSE)
  })
  list(rows = by_g, len = len, m = m, slice = slice, rslice = rslice,
       order = ord, n_seq = length(by_g))
}

#' Sequence structure, checks and the NA mask for one `hmm()` response.
#'
#' Everything here is DATA, resolved once: nothing below ever branches on
#' a parameter. The checks that must happen before the tape exists live
#' here rather than in the objective, because a refusal is only useful
#' while the user can still see their own call.
#'
#' @noRd
#' Addition terms and model shapes an `hmm()` response cannot carry.
#'
#' Called EARLY in the frame loop, before the generic per-aterm guards,
#' so a user who writes `cens()` on an HMM is told about the HMM rather
#' than about the family's missing CDF.
#'
#' @noRd
hmm_check_aterms <- function(resp, spec, av) {
  if (is.null(resp$family[["hmm"]])) return(invisible(NULL))
  if (length(spec$responses) > 1L || isTRUE(spec$rescor)) {
    stop("hmm() supports univariate models only: the forward recursion ",
         "is a likelihood over one response's sequences, and mvbf() / ",
         "rescor = TRUE would need a joint state process across ",
         "responses", call. = FALSE)
  }
  bad <- intersect(c("weights", "cens", "trunc_lb", "trunc_ub", "se"),
                   names(av))
  if (isTRUE(resp$aterms$mi)) bad <- c(bad, "mi")
  if (length(bad)) {
    lab <- c(weights = "weights()", cens = "cens()",
             trunc_lb = "trunc()", trunc_ub = "trunc()", se = "se()",
             mi = "mi()")[[bad[1L]]]
    stop("hmm() cannot be combined with ", lab, ": that term reshapes a ",
         "PER-ROW likelihood contribution, and an HMM's contribution is ",
         "per SEQUENCE - the forward recursion couples the rows of a ",
         "sequence and leaves no row-wise factor to weight, censor or ",
         "truncate. A missing response needs none of this: an NA is ",
         "kept and its emission masked, so the chain keeps its length",
         call. = FALSE)
  }
  invisible(NULL)
}

hmm_frame_block <- function(resp, spec, av, mf, y, n) {
  hs <- resp$family[["hmm"]]
  gv <- if (is.null(hs$group_expr)) {
    factor(rep(1L, n))
  } else {
    v <- eval(hs$group_expr, mf, resp$formula_env)
    if (anyNA(v)) {
      stop("hmm(): the sequences are defined by group = ",
           deparse1(hs$group_expr),
           ", so every row needs a group; that variable has ",
           sum(is.na(v)), " missing value(s)", call. = FALSE)
    }
    factor(v)
  }
  gidx <- as.integer(gv)
  tv <- if (is.null(hs$time_expr)) {
    idx <- integer(n)
    for (g in split(seq_len(n), gidx)) idx[g] <- seq_along(g)
    idx
  } else {
    v <- eval(hs$time_expr, mf, resp$formula_env)
    if (anyNA(v)) {
      stop("hmm(): the time variable '", deparse1(hs$time_expr),
           "' has missing values, so the order of the chain is ",
           "undefined at those rows", call. = FALSE)
    }
    if (is.factor(v)) match(as.character(v), levels(v)) else as.numeric(v)
  }
  key <- paste(gidx, tv, sep = "\r")
  if (anyDuplicated(key)) {
    dup <- key[duplicated(key)][1L]
    parts <- strsplit(dup, "\r", fixed = TRUE)[[1L]]
    stop("hmm(): time points must be unique within a sequence; group '",
         levels(gv)[as.integer(parts[1L])], "' has ", sum(key == dup),
         " rows at time '", parts[2L], "'. A Markov chain has one state ",
         "per time point, so repeated measurements at one time need ",
         "either a finer time variable or a grouping that separates ",
         "them", call. = FALSE)
  }
  st <- hmm_seq_structure(gidx, tv, n)
  tr_fixed <- vapply(hs$tr_names,
                     function(nm) !is.null(resp$dpars[[nm]]$constant), TRUE)
  if (max(st$len) < 2L && !all(tr_fixed)) {
    # With one observation per group no transition is ever taken, so the
    # transition logits are flat directions of the likelihood: the fit
    # converges, and reports a df (and an AIC) counting parameters the
    # data never touched (probe F2 measured df 7 against a mixture's 5).
    # The model IS a finite mixture at that point.
    stop("hmm(): every sequence has length 1, so no transition is ever ",
         "taken and the ", length(hs$tr_names), " transition ",
         "parameter(s) are flat directions of the likelihood - the fit ",
         "would report a df counting them. A model with one observation ",
         "per group is a finite mixture: use mixture(",
         paste(rep(hs$comp$family, min(hs$K, 2L)), collapse = ", "),
         if (hs$K > 2L) ", ..." else "",
         "), or hold every transition dpar at a constant (bf(..., ",
         hs$tr_names[1L], " = 0, ...)) if the degenerate chain is ",
         "deliberate", call. = FALSE)
  }

  const_trans <- all(vapply(hs$tr_names,
                            function(nm) hmm_dpar_is_constant(resp$dpars[[nm]]),
                            TRUE))
  if (hs$init == "stationary" && !const_trans) {
    varying <- hs$tr_names[!vapply(
      hs$tr_names,
      function(nm) hmm_dpar_is_constant(resp$dpars[[nm]]), TRUE)]
    stop("hmm(init = \"stationary\") needs a constant transition ",
         "matrix, and ", paste0("'", varying, "'", collapse = ", "),
         " carries a predictor. A chain whose transition matrix changes ",
         "from row to row has no single stationary distribution; use ",
         "init = \"estimated\" or init = \"uniform\"", call. = FALSE)
  }

  yv <- y
  miss <- if (is.matrix(yv)) {
    rowSums(is.na(yv)) > 0
  } else {
    is.na(yv)
  }
  mask <- NULL
  if (any(miss)) {
    if (all(miss)) {
      stop("hmm(): every response value is missing", call. = FALSE)
    }
    # An NA response is a time point the chain PASSES THROUGH without
    # emitting: the correct likelihood drops that step's emission factor
    # and keeps its transition. Dropping the row instead (which is what
    # na.omit does) shortens the chain and lets one transition stand in
    # for several - measurably biased (probe F4). The mask is data, so
    # multiplying the emission log-density by it puts no branch on the
    # tape; the placeholder value only has to be inside the family's
    # support, since its density is multiplied by zero.
    fill <- if (is.matrix(yv)) yv[which(!miss)[1L], ] else yv[!miss][1L]
    if (is.matrix(yv)) yv[miss, ] <- rep(fill, each = sum(miss))
    else yv[miss] <- fill
    mask <- as.numeric(!miss)
  }

  list(K = hs$K, init = hs$init, m = st$m, slice = st$slice,
       rslice = st$rslice, rows = st$rows, len = st$len,
       n_seq = st$n_seq, order = st$order,
       levels = levels(gv)[st$order], gindex = gidx, time = tv,
       const_trans = const_trans, miss = miss, mask = mask,
       n = n, y = yv)
}

## ---- the taped forward recursion ------------------------------------

#' Row-wise multinomial-logit transition probabilities, on the log scale.
#'
#' `out[[i]][[j]]` is `log P(state i -> state j)` at every row, with
#' state 1 the reference cell of every row (see the file header). The
#' normalizer is folded with `logspace_add` starting from the reference
#' cell's zero, so a large logit cannot overflow the way `log(1 +
#' sum(exp(eta)))` would.
#'
#' @noRd
hmm_log_tpm <- function(dp, K) {
  out <- vector("list", K)
  for (i in seq_len(K)) {
    lse <- 0 * dp[[hmm_tr_name(i, 2L)]]
    for (j in seq_len(K)[-1L]) {
      lse <- RTMB::logspace_add(lse, dp[[hmm_tr_name(i, j)]])
    }
    row <- vector("list", K)
    row[[1L]] <- -lse
    for (j in seq_len(K)[-1L]) row[[j]] <- dp[[hmm_tr_name(i, j)]] - lse
    out[[i]] <- row
  }
  out
}

#' The initial log-distribution as a length-K list of scalars.
#'
#' `"stationary"` solves `delta (I - Gamma + 1 1') = 1'` ON THE TAPE,
#' with `RTMB::solve` (the S4 advector method; `base::solve` would
#' silently take the numeric path). It is only reachable when the
#' transition matrix is constant, which the frame has already checked, so
#' reading row 1 of the per-row transition vectors is reading the whole
#' matrix.
#'
#' @noRd
hmm_log_delta <- function(hs, lg, K, extra) {
  "c" <- RTMB::ADoverload("c")
  switch(
    hs$init,
    uniform = as.list(rep(-log(K), K)),
    estimated = {
      ld <- extra[["hmm_ldel"]]
      lse <- RTMB::logspace_add(0 * ld[1L], ld[1L])
      for (k in seq_len(K - 1L)[-1L]) {
        lse <- RTMB::logspace_add(lse, ld[k])
      }
      out <- vector("list", K)
      out[[1L]] <- -lse
      for (k in seq_len(K - 1L)) out[[k + 1L]] <- ld[k] - lse
      out
    },
    stationary = {
      vals <- NULL
      for (j in seq_len(K)) {
        for (i in seq_len(K)) {
          e <- exp(lg[[i]][[j]][1L])
          vals <- if (is.null(vals)) e else c(vals, e)
        }
      }
      G <- RTMB::matrix(vals, K, K)
      A <- RTMB::matrix(rep(1, K * K), K, K) + diag(K) - G
      d <- as.vector(RTMB::solve(t(A), rep(1, K)))
      lapply(seq_len(K), function(k) log(d[k]))
    }
  )
}

#' The total log-likelihood of every sequence, on the tape.
#'
#' One step of the recursion is `K^2` vectorized `logspace_add` calls
#' over the sequences still running, not `K^2` scalar operations per row.
#'
#' @noRd
hmm_forward_ad <- function(lp, lg, ld, hg, K) {
  m <- hg[["m"]]
  S <- length(m)
  la <- lapply(seq_len(K), function(k) ld[[k]] + lp[[k]][hg$slice[[1L]]])
  total <- 0
  mprev <- m[1L]
  for (s in seq_len(S)) {
    if (s > 1L) {
      mm <- m[s]
      rp <- hg$slice[[s - 1L]][seq_len(mm)]
      if (mm < mprev) la <- lapply(la, function(v) v[seq_len(mm)])
      acc <- vector("list", K)
      for (j in seq_len(K)) {
        v <- lg[[1L]][[j]][rp] + la[[1L]]
        for (i in seq_len(K)[-1L]) {
          v <- RTMB::logspace_add(v, lg[[i]][[j]][rp] + la[[i]])
        }
        acc[[j]] <- v
      }
      rc <- hg$slice[[s]]
      la <- lapply(seq_len(K), function(k) acc[[k]] + lp[[k]][rc])
      mprev <- mm
    }
    mnext <- if (s < S) m[s + 1L] else 0L
    if (m[s] > mnext) {
      # the sequences that end at step s: sum out their final state
      idx <- (mnext + 1L):m[s]
      tv <- la[[1L]][idx]
      for (k in seq_len(K)[-1L]) tv <- RTMB::logspace_add(tv, la[[k]][idx])
      total <- total + sum(tv)
    }
  }
  total
}

#' The HMM contribution to the objective: the marginal log-likelihood of
#' every sequence, with the discrete states summed out exactly.
#'
#' @noRd
hmm_loglik_ad <- function(fam, dp, hg, av, extra) {
  hs <- fam[["hmm"]]
  K <- hs$K
  lp <- lapply(seq_len(K), function(k) {
    v <- hs$state_lpdf(hg[["y"]], dp, av, k)
    if (!is.null(hg[["mask"]])) v <- v * hg[["mask"]]
    v
  })
  lg <- hmm_log_tpm(dp, K)
  ld <- hmm_log_delta(hs, lg, K, extra)
  hmm_forward_ad(lp, lg, ld, hg, K)
}

## ---- numeric post-processing ----------------------------------------

#' Row-wise logsumexp of a matrix, overflow-safe.
#'
#' @noRd
hmm_lse_rows <- function(M) {
  # the row maximum by a K-step pmax fold rather than apply() or
  # as.data.frame(): this runs once per state per time slice, and on a
  # single long sequence every slice is one row, so the per-call
  # overhead is the whole cost of a decoding pass
  mx <- M[, 1L]
  nc <- ncol(M)
  if (nc > 1L) for (j in seq_len(nc)[-1L]) mx <- pmax(mx, M[, j])
  mx[!is.finite(mx)] <- 0
  mx + log(.rowSums(exp(M - mx), nrow(M), nc))
}

#' Numeric emission log-densities, log transition matrices and initial
#' log-distribution of a fitted hmm, all at the estimates.
#'
#' The same quantities the tape holds, recomputed off it: the tape is not
#' interrogable row by row, and the decoding passes are plain numeric
#' work that belongs outside the objective.
#'
#' @noRd
hmm_parts <- function(fit) {
  rspec <- uni_resp(fit, "hmm_probs()")
  fam <- rspec$family
  if (is.null(fam[["hmm"]])) {
    stop("This fit does not use an hmm() family, so it has no hidden ",
         "states to decode", call. = FALSE)
  }
  hg <- fit$frame[["hmm_g"]][[rspec$resp_name]]
  hs <- fam[["hmm"]]
  K <- hs$K
  dp <- eval_dpars(fit)[[rspec$resp_name]]
  av <- fit$frame$aterm_values[[rspec$resp_name]]
  n <- hg[["n"]]
  lpmat <- matrix(vapply(seq_len(K), function(k) {
    v <- as.numeric(hs$state_lpdf(hg[["y"]], dp, av, k))
    v <- rep(v, length.out = n)
    if (!is.null(hg[["mask"]])) v <- v * hg[["mask"]]
    v
  }, numeric(n)), n, K)
  # log transition probabilities as a K x K list of length-n vectors
  lg <- vector("list", K)
  for (i in seq_len(K)) {
    E <- matrix(vapply(seq_len(K), function(j) {
      if (j == 1L) rep(0, n) else {
        rep(as.numeric(dp[[hmm_tr_name(i, j)]]), length.out = n)
      }
    }, numeric(n)), n, K)
    lg[[i]] <- E - hmm_lse_rows(E)
  }
  ld <- switch(
    hs$init,
    uniform = rep(-log(K), K),
    estimated = {
      v <- c(0, as.numeric(fit$estimates[["hmm_ldel"]]))
      v - (max(v) + log(sum(exp(v - max(v)))))
    },
    stationary = {
      G <- t(vapply(lg, function(Li) exp(Li[1L, ]), numeric(K)))
      log(as.vector(solve(t(matrix(1, K, K) + diag(K) - G), rep(1, K))))
    }
  )
  list(rspec = rspec, fam = fam, hs = hs, K = K, hg = hg, dp = dp,
       av = av, n = n, lpmat = lpmat, lg = lg, ld = ld)
}

#' Smoothed state probabilities by forward-backward.
#'
#' Both passes are sliced by time exactly as the taped forward pass is,
#' so a data set of many short sequences costs `O(Tmax K^2)` R-level
#' operations rather than `O(n K^2)`. `alpha` and `beta` are scattered
#' back into `n x K` log matrices; the row-wise normalization at the end
#' is what divides by the sequence's likelihood.
#'
#' @noRd
hmm_fb <- function(p) {
  hg <- p$hg
  K <- p$K
  m <- hg[["m"]]
  S <- length(m)
  LA <- matrix(-Inf, p[["n"]], K)
  LB <- matrix(-Inf, p[["n"]], K)

  la <- matrix(p$ld, m[1L], K, byrow = TRUE) + p$lpmat[hg$slice[[1L]], ,
                                                       drop = FALSE]
  LA[hg$slice[[1L]], ] <- la
  for (s in seq_len(S)[-1L]) {
    mm <- m[s]
    rp <- hg$slice[[s - 1L]][seq_len(mm)]
    lap <- la[seq_len(mm), , drop = FALSE]
    acc <- matrix(0, mm, K)
    for (j in seq_len(K)) {
      # matrix(), not the vapply result: with ONE live sequence vapply
      # gives a length-K vector and every row-wise step below silently
      # changes meaning
      M <- matrix(vapply(seq_len(K),
                         function(i) p$lg[[i]][rp, j] + lap[, i],
                         numeric(mm)), mm, K)
      acc[, j] <- hmm_lse_rows(M)
    }
    la <- acc + p$lpmat[hg$slice[[s]], , drop = FALSE]
    LA[hg$slice[[s]], ] <- la
  }

  lb <- matrix(0, m[1L], K)
  LB[hg$rslice[[1L]], ] <- lb
  for (u in seq_len(S)[-1L]) {
    mm <- m[u]
    rt <- hg$rslice[[u]][seq_len(mm)]
    rt1 <- hg$rslice[[u - 1L]][seq_len(mm)]
    lbp <- lb[seq_len(mm), , drop = FALSE]
    em <- p$lpmat[rt1, , drop = FALSE]
    new <- matrix(0, mm, K)
    for (i in seq_len(K)) {
      M <- matrix(vapply(seq_len(K),
                         function(j) p$lg[[i]][rt, j] + em[, j] + lbp[, j],
                         numeric(mm)), mm, K)
      new[, i] <- hmm_lse_rows(M)
    }
    lb <- new
    LB[rt, ] <- lb
  }

  L <- LA + LB
  mx <- L[, 1L]
  if (K > 1L) for (j in seq_len(K)[-1L]) mx <- pmax(mx, L[, j])
  P <- exp(L - mx)
  P / .rowSums(P, p[["n"]], K)
}

#' Posterior state probabilities of an hmm fit
#'
#' The smoothed occupancy `P(S_t = k | y)` for every row of the data,
#' from a forward-backward pass at the estimates. This is the [mixture()]
#' analog [mixture_probs()], one rung up: an HMM's per-row state
#' probability conditions on the WHOLE sequence, not on that row alone,
#' which is why it needs a backward pass and cannot come out of the
#' family's density.
#'
#' Under a random-effect model the probabilities are conditional on the
#' random-effect modes, as every other post-fit quantity in the package
#' is. Rows whose response is `NA` still get a probability: the chain
#' passes through them and the neighboring observations inform the state.
#'
#' @param fit A `frmtmb_fit` with an [hmm()] family.
#' @return An `n x K` matrix of probabilities whose rows sum to one, with
#'   columns named `state1 ... stateK`.
#' @seealso [hmm_viterbi()], [hmm()]
#' @examples
#' set.seed(12)
#' dd <- data.frame(id = 1, t = 1:120)
#' s <- integer(120); s[1] <- 1L
#' G <- matrix(c(0.92, 0.08, 0.15, 0.85), 2, 2, byrow = TRUE)
#' for (t in 2:120) s[t] <- sample.int(2, 1, prob = G[s[t - 1], ])
#' dd$y <- rnorm(120, c(0, 3)[s], 0.5)
#' fit <- frm(bf(y ~ 1),
#'            family = hmm(K = 2, gaussian(), time = t, group = id),
#'            data = dd)
#' p <- hmm_probs(fit)
#' head(p)
#' rowSums(p)[1:3]
#' mean(max.col(p) == s)
#' @export
hmm_probs <- function(fit) {
  p <- hmm_parts(fit)
  P <- hmm_fb(p)
  colnames(P) <- paste0("state", seq_len(p$K))
  P
}

#' Most likely state path of an hmm fit (Viterbi)
#'
#' Global decoding: the single state sequence with the highest posterior
#' probability, per sequence, by the Viterbi recursion. This is not the
#' same as taking the most likely state at each row separately
#' (`max.col(hmm_probs(fit))`, local decoding), and the two can disagree
#' - only the Viterbi path is guaranteed to be a possible path under the
#' transition matrix.
#'
#' State LABELS are arbitrary: the likelihood is invariant to permuting
#' the states, so which run calls a state "1" is not meaningful across
#' fits. Compare paths after matching states by their fitted means.
#'
#' @param fit A `frmtmb_fit` with an [hmm()] family.
#' @return An integer vector, one state index per row of the data.
#' @seealso [hmm_probs()], [hmm()]
#' @examples
#' set.seed(12)
#' dd <- data.frame(id = 1, t = 1:120)
#' s <- integer(120); s[1] <- 1L
#' G <- matrix(c(0.92, 0.08, 0.15, 0.85), 2, 2, byrow = TRUE)
#' for (t in 2:120) s[t] <- sample.int(2, 1, prob = G[s[t - 1], ])
#' dd$y <- rnorm(120, c(0, 3)[s], 0.5)
#' fit <- frm(bf(y ~ 1),
#'            family = hmm(K = 2, gaussian(), time = t, group = id),
#'            data = dd)
#' v <- hmm_viterbi(fit)
#' table(v, truth = s)
#' @export
hmm_viterbi <- function(fit) {
  p <- hmm_parts(fit)
  hg <- p$hg
  K <- p$K
  m <- hg[["m"]]
  S <- length(m)
  delta <- vector("list", S)
  bp <- vector("list", S)
  d <- matrix(p$ld, m[1L], K, byrow = TRUE) +
    p$lpmat[hg$slice[[1L]], , drop = FALSE]
  delta[[1L]] <- d
  for (s in seq_len(S)[-1L]) {
    mm <- m[s]
    rp <- hg$slice[[s - 1L]][seq_len(mm)]
    dp_ <- d[seq_len(mm), , drop = FALSE]
    acc <- matrix(0, mm, K)
    ptr <- matrix(0L, mm, K)
    for (j in seq_len(K)) {
      M <- matrix(vapply(seq_len(K),
                         function(i) p$lg[[i]][rp, j] + dp_[, i],
                         numeric(mm)), mm, K)
      ptr[, j] <- max.col(M, ties.method = "first")
      acc[, j] <- M[cbind(seq_len(mm), ptr[, j])]
    }
    d <- acc + p$lpmat[hg$slice[[s]], , drop = FALSE]
    delta[[s]] <- d
    bp[[s]] <- ptr
  }
  out <- integer(p[["n"]])
  cur <- integer(m[1L])
  for (s in rev(seq_len(S))) {
    mm <- m[s]
    mnext <- if (s < S) m[s + 1L] else 0L
    if (mm > mnext) {
      # sequences that END at step s start their backtrack here
      idx <- (mnext + 1L):mm
      cur[idx] <- max.col(delta[[s]][idx, , drop = FALSE],
                          ties.method = "first")
    }
    out[hg$slice[[s]]] <- cur[seq_len(mm)]
    if (s > 1L) {
      cur[seq_len(mm)] <- bp[[s]][cbind(seq_len(mm), cur[seq_len(mm)])]
    }
  }
  out
}

#' The occupancy-weighted response mean of an hmm fit,
#' `E[y_t | y] = sum_k P(S_t = k | y) mu_k(x_t)`.
#'
#' This is the quantity `fitted()`, `predict(type = "response")` and the
#' response/pearson residuals all report. There is no per-row `mean_fn`
#' that could produce it, which is exactly why rung 1 of this feature
#' silently reported state 1's mean everywhere.
#'
#' @noRd
hmm_mean_response <- function(fit) {
  p <- hmm_parts(fit)
  P <- hmm_fb(p)
  out <- 0
  for (k in seq_len(p$K)) {
    mk <- p$hs$state_mean(p$dp, p$av, k)
    if (is.null(mk)) {
      stop("The state-dependent family '", p$hs$comp$family, "' of this ",
           "hmm() fit has no mean, so the occupancy-weighted expected ",
           "response is undefined", call. = FALSE)
    }
    out <- out + P[, k] * rep(as.numeric(mk), length.out = p[["n"]])
  }
  out
}

#' Per-row state variance of an hmm fit, for pearson residuals: the law
#' of total variance over the occupancy distribution, so the between-state
#' spread counts as well as the within-state one.
#'
#' @noRd
hmm_var_response <- function(fit) {
  p <- hmm_parts(fit)
  if (is.null(p$hs$comp$post$var_fn)) return(NULL)
  P <- hmm_fb(p)
  mu <- matrix(0, p[["n"]], p$K)
  vv <- matrix(0, p[["n"]], p$K)
  for (k in seq_len(p$K)) {
    dk <- p$hs$state_dpars(p$dp, k)
    mu[, k] <- rep(as.numeric(response_mean(p$hs$comp, dk, p$av)),
                   length.out = p[["n"]])
    vv[, k] <- rep(as.numeric(p$hs$comp$post$var_fn(dk, p$av)),
                   length.out = p[["n"]])
  }
  m <- rowSums(P * mu)
  rowSums(P * (vv + mu^2)) - m^2
}

#' Forward-simulate a state path per sequence, then one emission draw per
#' row from that row's state.
#'
#' @noRd
hmm_simulate_rows <- function(fit) {
  p <- hmm_parts(fit)
  hg <- p$hg
  K <- p$K
  m <- hg[["m"]]
  S <- length(m)
  st <- integer(p[["n"]])
  cur <- integer(m[1L])
  d0 <- exp(p$ld)
  cur <- sample.int(K, m[1L], replace = TRUE, prob = d0)
  st[hg$slice[[1L]]] <- cur
  for (s in seq_len(S)[-1L]) {
    mm <- m[s]
    rp <- hg$slice[[s - 1L]][seq_len(mm)]
    prev <- cur[seq_len(mm)]
    nxt <- integer(mm)
    for (i in seq_len(K)) {
      sel <- which(prev == i)
      if (!length(sel)) next
      Pi <- exp(p$lg[[i]][rp[sel], , drop = FALSE])
      nxt[sel] <- vapply(seq_along(sel), function(z) {
        sample.int(K, 1L, prob = Pi[z, ])
      }, integer(1))
    }
    cur <- nxt
    st[hg$slice[[s]]] <- cur
  }
  ys <- numeric(p[["n"]])
  for (k in seq_len(K)) {
    rows <- which(st == k)
    if (length(rows)) {
      ys[rows] <- p$hs$state_sim(p$dp, p$av, p[["n"]], k, rows)
    }
  }
  ys
}

## ---- fit-time guards -------------------------------------------------

#' The refusals an `hmm()` fit needs from the FITTING options, and the
#' one warning about the starting values.
#'
#' REML is the important one. `frm(REML = TRUE)` integrates the
#' `primary_dpars` fixed effects, which for an HMM is the state means
#' only: the transition blocks and the state dispersions stay in the
#' outer problem. That is a PARTIAL restricted likelihood corresponding
#' to no standard definition, and it ran silently before this guard
#' existed (probe F5).
#'
#' @noRd
hmm_check_fit <- function(spec, frame, template, REML, quadrature,
                          control) {
  if (!has_hmm(spec)) return(invisible(NULL))
  if (isTRUE(REML)) {
    stop("REML = TRUE cannot be combined with hmm(): REML integrates ",
         "out the location coefficients, which here are one set per ",
         "hidden state, and leaves the transition and dispersion ",
         "parameters outside. The result is a partial restricted ",
         "likelihood matching no standard definition. Use REML = FALSE",
         call. = FALSE)
  }
  if (isTRUE(quadrature)) {
    stop("quadrature = TRUE cannot be combined with hmm(): the rule ",
         "integrates a random effect against a PRODUCT of per-row ",
         "densities, and an HMM's likelihood is a forward recursion ",
         "over each sequence, not such a product. Use quadrature = ",
         "FALSE (Laplace) and check_laplace() to judge the ",
         "approximation", call. = FALSE)
  }
  if (isTRUE(control$profile)) {
    stop("frmtmb_control(profile = TRUE) cannot be combined with hmm(): ",
         "profiling moves the state means into an inner Laplace ",
         "problem, and an HMM likelihood is multimodal in them (the ",
         "states are exchangeable). Use profile = FALSE", call. = FALSE)
  }
  hmm_warn_symmetric_start(spec, frame, template)
  invisible(NULL)
}

#' A start with every state's location predictor equal sits ON the
#' label-symmetry axis: the gradient is invariant under a label swap, the
#' optimizer never leaves the diagonal, and the fit converges to the
#' one-state solution (probe F1 measured 262 log-likelihood units below
#' the optimum, with no diagnostic firing). The default inits are spread
#' quantiles and never do this; a user-supplied `start` can.
#'
#' @noRd
hmm_warn_symmetric_start <- function(spec, frame, template) {
  for (rn in names(spec$responses)) {
    hs <- spec$responses[[rn]]$family[["hmm"]]
    if (is.null(hs)) next
    prim <- spec$responses[[rn]]$family[["primary_dpars"]]
    vals <- vapply(prim, function(dnm) {
      lp <- frame$linpreds[[linpred_key(rn, dnm)]]
      if (is.null(lp) || !ncol(lp$X)) return(NA_real_)
      sum(lp$X[1L, ] * template[[lp$par]][lp$idx])
    }, numeric(1))
    vals <- vals[is.finite(vals)]
    if (length(vals) > 1L && diff(range(vals)) < 1e-8) {
      warning("hmm(): every state's location predictor starts at the ",
              "same value (", format(vals[1L], digits = 4),
              "). That start is a fixed point of the label symmetry: ",
              "the optimizer cannot separate the states from it and ",
              "the fit will collapse to a one-state solution. Spread ",
              "the starting intercepts, or drop `start` and let the ",
              "response-quantile defaults be used", call. = FALSE)
    }
  }
}
