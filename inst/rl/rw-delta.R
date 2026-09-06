# Rescorla-Wagner delta learning on a two-armed bandit, written as a
# structured frmtmb family. vignette("reinforcement-learning") shows
# every chunk below in place; tests/testthat/test-rl-example.R fits the
# same code. Load it with:
#
#   source(system.file("rl", "rw-delta.R", package = "frmtmb"))
#
# The model is the one hBayesDM calls bandit2arm_delta (Ahn, Haines and
# Zhang 2017). Q starts at 0 for both arms. After each choice
# Q[chosen] = Q[chosen] + alpha * (reward - Q[chosen]), and the choice
# rule is P(arm 1) = plogis(beta * (Q1 - Q2)). No code here comes from
# hBayesDM: the equations are the published model.

## ---- rl-aterm ----
# reward() carries what each arm WOULD have paid on this trial, in arm
# order. The likelihood only ever reads the chosen arm's entry, so data
# that records the received outcome alone passes it twice. The second
# column is what makes the simulator coherent, because a simulated
# choice needs the payoff of the arm the subject did not take in the
# data. vreal(pay1, pay2) would carry the same two columns. The registry
# is what lets a family name its data in the words of its own
# literature, and `reward` is the word this literature uses.
frmtmb_register_aterm("reward", arity = 2L)

## ---- rl-block ----
# The frame block is DATA, resolved once, and it is what turns a long
# data frame into the subject-by-trial layout the recursion walks. Every
# subject gets a row of `idx` holding its row numbers in trial order,
# padded on the right to the longest subject with a repeat of its own
# first row. `mask` is 1 on a real trial and 0 on a pad. It multiplies
# both the log-likelihood and the learning rate, so a padded cell
# contributes nothing and updates nothing.
rw_block <- function(resp, spec, av, mf, y, n) {
  rw <- resp$family[["rw"]]
  gv <- eval(rw[["subject_expr"]], mf, resp$formula_env)
  if (anyNA(gv)) {
    stop("rw_delta(): every row needs a subject, and ",
         deparse1(rw[["subject_expr"]]), " has ", sum(is.na(gv)),
         " missing value(s)", call. = FALSE)
  }
  gv <- factor(gv)
  tv <- if (is.null(rw[["trial_expr"]])) {
    # no trial column: the data frame's own row order within a subject
    out <- integer(n)
    for (g in split(seq_len(n), gv)) out[g] <- seq_along(g)
    out
  } else {
    v <- eval(rw[["trial_expr"]], mf, resp$formula_env)
    if (anyNA(v)) {
      stop("rw_delta(): the trial variable '", deparse1(rw[["trial_expr"]]),
           "' has missing values, so the order of the recursion is ",
           "undefined at those rows", call. = FALSE)
    }
    as.numeric(v)
  }
  rows <- lapply(split(seq_len(n), gv), function(r) r[order(tv[r])])
  key <- paste(as.integer(gv), tv, sep = "|")
  if (anyDuplicated(key)) {
    dup <- key[duplicated(key)][1L]
    stop("rw_delta(): trial numbers must be unique within a subject; ",
         sum(key == dup), " rows share one. A delta rule updates once ",
         "per trial, so two rows at one trial have no order to learn in",
         call. = FALSE)
  }
  len <- lengths(rows)
  nt <- max(len)
  # The matrix() before each transpose is load-bearing: vapply() drops to
  # a plain vector when nt == 1, and t() of a vector is 1-by-n_subj, so
  # the block would come back with every subject read as one more trial
  # of a single learner and Q carried across subject boundaries.
  idx <- t(matrix(vapply(rows, function(r) c(r, rep(r[1L], nt - length(r))),
                         integer(nt)), nrow = nt))
  mask <- t(matrix(vapply(len, function(l) as.numeric(seq_len(nt) <= l),
                          numeric(nt)), nrow = nt))
  # `group` is a reserved name: it is what the per-subject
  # log-likelihood returns one value per, in level order, and what the
  # core checks a model's own grouping factor against before it lets
  # frm(importance =) sum this family's pieces into its groups.
  list(idx = idx, mask = mask, len = len, n_subj = length(rows),
       n_trial = nt, group = gv, subject = gv, trial = tv,
       levels = levels(gv), n = n)
}

## ---- rl-loglik ----
# Broadcast a dpar the core handed back as a single value. rep(v, each =)
# strips the advector class (see the RTMB notes in dev/); multiplying by
# a vector of ones does not.
rw_bcast <- function(v, n) if (length(v) == 1L) v * rep(1, n) else v

# The subject-by-trial layout, repeated once per replicate of the
# design. The importance correction evaluates the whole design once per
# draw, stacked, so the recursion has to run over subject-crossed-with-
# draw rather than over subject. Because the loop below is already
# vectorized across subjects, that is the same loop over a longer
# vector: `nrep` blocks of subjects, each block's row numbers shifted by
# a whole design.
rw_stack <- function(block, nrep) {
  idx <- block[["idx"]]
  msk <- block[["mask"]]
  if (nrep == 1L) return(list(idx = idx, mask = msk))
  take <- rep(seq_len(nrow(idx)), nrep)
  list(idx = idx[take, , drop = FALSE] +
         rep((seq_len(nrep) - 1L) * block[["n"]], each = nrow(idx)),
       mask = msk[take, , drop = FALSE])
}

# The trial-by-trial factors of the likelihood, as one vector per trial
# over all (subject, replicate) pairs. Every quantity the family reports
# is built from these, so there is one recursion and no second
# definition of the same arithmetic to drift from it.
#
# One iteration per TRIAL, each one a handful of vector operations over
# all subjects at once. Writing it the other way round, one iteration
# per row with Q[s] <- ... inside, is the shape RTMB punishes: an
# element-wise assignment into a taped vector costs about three orders
# of magnitude more than the vector operation it replaces.
rw_terms <- function(y, dpars, aterms, block) {
  n <- block[["n"]]
  nrep <- NROW(y) %/% n
  sk <- rw_stack(block, nrep)
  idx <- sk[["idx"]]
  msk <- sk[["mask"]]
  nn <- n * nrep
  alpha <- rw_bcast(dpars[["alpha"]], nn)
  beta <- rw_bcast(dpars[["beta"]], nn)
  pay1 <- aterms[["reward1"]]
  pay2 <- aterms[["reward2"]]
  q1 <- rep(0, nrow(idx))
  q2 <- q1
  out <- vector("list", ncol(idx))
  for (t in seq_len(ncol(idx))) {
    i <- idx[, t]
    m <- msk[, t]
    c1 <- y[i]
    eta <- beta[i] * (q1 - q2)
    # log P(choice) = c1 * eta - log(1 + exp(eta)), folded with
    # logspace_add so a decisive subject does not overflow
    out[[t]] <- m * (c1 * eta - RTMB::logspace_add(0 * eta, eta))
    # the mask enters through the learning rate, so a padded cell learns
    # nothing and no branch reaches the tape
    a <- alpha[i] * m
    q1 <- q1 + (a * c1) * (pay1[i] - q1)
    q2 <- q2 + (a * (1 - c1)) * (pay2[i] - q2)
  }
  list(terms = out, idx = idx, mask = msk)
}

# One value per subject, which is this family's coarsest honest piece
# and the one frm(importance =) resamples. A subject's trials are a
# product, so its log-likelihood is the elementwise sum of the per-trial
# vectors: the piece is free, given the recursion.
#
# weights() is refused in check_spec, so `weights` is 1 in all three
# slots. A trial's factor cannot be reweighted on its own: its value
# depends on every earlier trial of the same subject.
rw_loglik_group <- function(y, dpars, aterms, weights, block, extra) {
  Reduce(`+`, rw_terms(y, dpars, aterms, block)[["terms"]])
}

# One value per trial: the choice probability CONDITIONAL on everything
# the subject saw before it. These are the factors themselves, scattered
# back to the rows they came from, and a Bernoulli factor's saturated
# value is zero, which is what makes a deviance residual definable here.
rw_loglik_row <- function(y, dpars, aterms, weights, block, extra) {
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  tm <- rw_terms(y, dpars, aterms, block)
  # ONE sub-assignment, not one per trial: `[<-` copies the whole vector
  # each time it is called on a taped one. The padded cells are dropped
  # rather than written, because a pad repeats its subject's FIRST row
  # number and would otherwise overwrite that trial's own factor with a
  # masked zero.
  keep <- as.vector(tm[["mask"]]) == 1
  out <- rep(0, length(y))
  out[as.vector(tm[["idx"]])[keep]] <-
    do.call(c, tm[["terms"]])[keep]
  # off the tape only: the numeric path is the one residuals() reads,
  # and an advector carries no attributes
  if (!inherits(out, "advector")) {
    attr(out, "saturated") <- rep(0, length(out))
  }
  out
}

# The whole response's log-likelihood, defined as the sum of the pieces
# rather than accumulated separately, so the three slots cannot disagree.
rw_loglik <- function(y, dpars, aterms, weights, block, extra) {
  sum(rw_loglik_group(y, dpars, aterms, weights, block, extra))
}

## ---- rl-replay ----
# The same recursion off the tape, in plain numbers, at the estimates.
# It answers two questions: the per-trial choice probability, which is
# the fitted value, and a forward draw, which is the simulator. They
# share one function because a draw is a replay whose choices come from
# rbinom() instead of from the data.
rw_replay <- function(dpars, aterms, block, y, draw = FALSE) {
  idx <- block[["idx"]]
  msk <- block[["mask"]]
  n <- block[["n"]]
  alpha <- rep(as.numeric(dpars[["alpha"]]), length.out = n)
  beta <- rep(as.numeric(dpars[["beta"]]), length.out = n)
  pay1 <- aterms[["reward1"]]
  pay2 <- aterms[["reward2"]]
  p <- rep(NA_real_, n)
  q1 <- rep(0, nrow(idx))
  q2 <- q1
  for (t in seq_len(ncol(idx))) {
    i <- idx[, t]
    m <- msk[, t]
    keep <- m == 1
    pt <- stats::plogis(beta[i] * (q1 - q2))
    p[i[keep]] <- pt[keep]
    c1 <- if (draw) stats::rbinom(length(i), 1L, pt) else y[i]
    if (draw) y[i[keep]] <- c1[keep]
    a <- alpha[i] * m
    q1 <- q1 + (a * c1) * (pay1[i] - q1)
    q2 <- q2 + (a * (1 - c1)) * (pay2[i] - q2)
  }
  list(p = p, y = y)
}

# The three post-fit slots, each one line of work on top of the replay.
rw_parts <- function(fit) {
  rspec <- single_response(fit, "a rw_delta() fit")
  list(dpars = eval_dpars(fit)[[rspec$resp_name]],
       aterms = fit$frame[["aterm_values"]][[rspec$resp_name]],
       y = fit$frame[["y"]][[rspec$resp_name]])
}

rw_fitted <- function(fit, block) {
  pp <- rw_parts(fit)
  rw_replay(pp$dpars, pp$aterms, block, pp$y)$p
}

## ---- rl-sim ----
# The structured simulator. simulate(), posterior_predict() and
# frm_simulate() all reach it through the same context, so the draw is
# written once. It walks each subject forward, drawing a choice from the
# current Q pair and learning from the payoff that choice earns.
rw_sim <- function(ctx) {
  blk <- ctx[["block"]]
  rw_replay(ctx[["dpars"]], ctx[["aterms"]], blk,
            y = rep(0, blk[["n"]]), draw = TRUE)$y
}

## ---- rl-family ----
# The family. Two distributional parameters, each with its own linear
# predictor: alpha on (0, 1) through a logit link, beta on (0, Inf)
# through a log link. `primary_dpars = "alpha"` sends the main
# right-hand side of the formula to the learning rate.
rw_delta <- function(subject, trial = NULL) {
  subject_expr <- substitute(subject)
  trial_expr <- substitute(trial)
  if (is.null(subject_expr)) {
    stop("rw_delta(subject =) names the column that separates one ",
         "learner's trial sequence from the next", call. = FALSE)
  }
  fam <- frmtmb_family(
    "rw_delta",
    dpars = c("alpha", "beta"),
    links = list(alpha = "logit", beta = "log"),
    primary_dpars = "alpha",
    type = "discrete",
    required_aterms = c("reward1", "reward2"),
    # The likelihood factorizes over trials but is not ROWWISE: a
    # trial's factor reads the whole history of its subject, which this
    # signature cannot see. Returning something anyway is the silent lie
    # the structured protocol exists to remove.
    lpdf = function(y, dpars, aterms, extra = NULL) {
      stop("A rw_delta() trial's probability depends on every earlier ",
           "trial of the same subject, so the family has no row-wise ",
           "log-density. Use logLik() for the total, or fitted() for ",
           "the per-trial choice probabilities", call. = FALSE)
    },
    valid_y = function(y, aterms) {
      if (!all(y %in% c(0, 1))) {
        stop("rw_delta(): the response is the arm chosen on each trial, ",
             "coded 1 for the first arm and 0 for the second",
             call. = FALSE)
      }
    },
    init_dpars = list(alpha = function(y, aterms) 0.3,
                      beta = function(y, aterms) 1),
    structure = rw_structure()
  )
  fam[["rw"]] <- list(subject_expr = subject_expr, trial_expr = trial_expr)
  fam
}

## ---- rl-structure ----
# The structure is the family's half of the protocol: where its data
# block comes from, what its likelihood is, what it can answer
# afterwards, and what it refuses. Every capability flag starts FALSE
# for a likelihood that does not factorize over rows, so the refusals
# below are explanations, not switches.
rw_structure <- function() {
  frmtmb_structure(
    frame_vars = function(fam) {
      list(fam[["rw"]][["subject_expr"]], fam[["rw"]][["trial_expr"]])
    },
    # An unanswered trial produces no prediction error, so dropping its
    # row IS the right recursion. This is where a delta rule differs
    # from a hidden Markov chain, which keeps the row because the state
    # still moves; see frmtmb.latent::hmm().
    keep_na = FALSE,
    check_spec = function(resp, spec, av) {
      if (length(spec$responses) > 1L || isTRUE(spec$rescor)) {
        stop("rw_delta() supports univariate models only: the recursion ",
             "is a likelihood over one response's trial sequences",
             call. = FALSE)
      }
      bad <- intersect(c("weights", "cens", "trunc_lb", "trunc_ub", "se"),
                       names(av))
      if (length(bad)) {
        stop("rw_delta() cannot be combined with ", bad[1L], "(): that ",
             "term reshapes a per-row likelihood contribution, and a ",
             "trial's contribution here is conditional on every earlier ",
             "trial of the same subject", call. = FALSE)
      }
    },
    frame_block = rw_block,
    loglik = rw_loglik,
    # The two factorizations. `unit` is a different question and keeps
    # its own answer: the finest factorization here is the TRIAL, while
    # the smallest thing that can honestly be left out is a whole
    # subject, because dropping one trial changes every later trial's Q.
    loglik_row = rw_loglik_row,
    loglik_group = rw_loglik_group,
    unit = "one subject's trial sequence",
    fitted_mean = rw_fitted,
    fitted_var = function(fit, block) {
      p <- rw_fitted(fit, block)
      p * (1 - p)
    },
    sim_ctx = rw_sim,
    # The one capability a factorized likelihood buys back here. Every
    # other flag stays FALSE: a per-trial factor is not a per-trial
    # MODEL, so newdata, conditional effects and one-step-ahead
    # residuals are refused for the reasons below and not for want of
    # the factors.
    supports = list(deviance = TRUE),
    refusals = list(
      newdata_response = paste0(
        "predict(type = 'response') on a rw_delta() fit is not ",
        "available for newdata: a trial's choice probability is ",
        "conditional on the payoffs and choices of every earlier trial ",
        "of the same subject, and newdata carries no block to replay. ",
        "Use fitted() on the training data, or predict(dpar = 'alpha') ",
        "for a learning rate"),
      conditional_effects = paste0(
        "conditional_effects() is not available for a rw_delta() fit: ",
        "the expected response is a choice probability that depends on ",
        "a whole trial history, which the synthetic grid this function ",
        "builds does not have. Plot the learning rate itself with ",
        "predict(dpar = 'alpha')"),
      osa = paste0(
        "residuals(type = 'osa') is not available for a rw_delta() fit: ",
        "the tape holds the recursion over each whole sequence with no ",
        "registered observation vector. Use type = 'pearson', which ",
        "divides by the binomial variance of the choice probability")
    )
  )
}

## ---- rl-design ----
# Scaffolding for the worked example, not part of the family. It builds
# a two-armed bandit with the payoff schedule of BOTH arms fixed in
# advance, one row per subject and trial, and a placeholder response the
# simulator overwrites. Arm 1 pays with probability p1 and arm 2 with
# probability p2, so the learnable contrast is p1 - p2.
rl_bandit_design <- function(n_subj, n_trial, p1 = 0.7, p2 = 0.3) {
  dd <- expand.grid(trial = seq_len(n_trial), id = seq_len(n_subj))
  dd$condition <- factor(rep(0:1, length.out = n_subj)[dd$id], 0:1,
                         c("ctl", "trt"))
  dd$id <- factor(dd$id)
  dd$pay1 <- stats::rbinom(nrow(dd), 1L, p1)
  dd$pay2 <- stats::rbinom(nrow(dd), 1L, p2)
  dd$choice <- 0
  # expand.grid() leaves an out.attrs attribute that str() prints at
  # length and nothing reads
  attr(dd, "out.attrs") <- NULL
  dd
}

# The learning rate carries the condition effect; the softmax
# temperature is a per-subject intercept correlated with it.
rl_bform <- bf(choice | reward(pay1, pay2) ~ condition + (1 | p | id),
               beta ~ 1 + (1 | p | id))

rl_family <- function() rw_delta(subject = id, trial = trial)

# Parameter values in frm_simulate()'s natural spelling: the same names
# fixef() and VarCorr() report back.
rl_truth <- list(
  alpha_Intercept = stats::qlogis(0.35),
  alpha_conditiontrt = 0.8,
  beta_Intercept = log(3),
  sd_id__choice.alphaIntercept = 0.5,
  sd_id__choice.betaIntercept = 0.4,
  cor_id__choice.alphaIntercept__choice.betaIntercept = 0.3
)

# Draws come from the family's own simulator: frm_simulate() draws a new
# set of subject effects from `truth`, then walks each subject forward.
rl_simulate <- function(design, truth = rl_truth, nsim = 1, seed = NULL) {
  sims <- frm_simulate(rl_bform, data = design, family = rl_family(),
                       newparams = truth, nsim = nsim, seed = seed)
  lapply(seq_len(nsim), function(s) {
    d <- design
    d$choice <- sims[[s]]
    d
  })
}
