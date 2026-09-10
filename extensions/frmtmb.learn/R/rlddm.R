#' Reinforcement learning with a drift-diffusion choice rule
#'
#' A delta learning rule feeding the drift rate of a two-boundary Wiener
#' diffusion, so that the choice and the response time are one
#' likelihood rather than two. The learning half is
#' [bandit2arm_delta()]'s exactly:
#'
#' ```
#' Q[chosen] <- Q[chosen] + alpha * (reward - Q[chosen])
#' ```
#'
#' and the choice half replaces the softmax with a first-passage
#' density. Evidence accumulates at a rate proportional to how far apart
#' the two value estimates have grown,
#'
#' ```
#' v <- drift * (Q[upper] - Q[lower])
#' ```
#'
#' between boundaries `bs` apart, starting a fraction `bias` of the way
#' up, and the response time is `ndt` plus the time the accumulator
#' takes to touch a boundary. A trial's contribution is then the joint
#' density of WHICH boundary was touched and WHEN, which is what makes
#' this a model of response times rather than a choice model with a
#' timing footnote: a learning rate estimated from choices alone cannot
#' tell a subject who deliberated from one who guessed, and this one
#' can.
#'
#' The model is Pedersen, Frank and Biele (2017). hBayesDM does not
#' carry it; see the section below for what it does carry.
#'
#' @section What it corresponds to in hBayesDM:
#' Nothing, and that is worth stating plainly because every other family
#' in this package has a counterpart. hBayesDM ships `choiceRT_ddm`,
#' which is this family's choice rule with NO learning: one drift rate
#' per subject and no value store. It also ships `bandit2arm_delta`,
#' which is this family's learning rule with a softmax instead of a
#' diffusion. `rlddm()` is the two joined, and the join is the model
#' rather than a renaming of either.
#'
#' Where the names do map, they map as hBayesDM's `choiceRT_ddm` spells
#' the diffusion: its `alpha` is the boundary separation, which is `bs`
#' here, its `beta` is the start point, which is `bias`, its `delta` is
#' the drift rate, which is `drift` times the value difference, and its
#' `tau` is the non-decision time, which is `ndt`. `alpha` here is the
#' LEARNING RATE, as it is in every other family in this package. That
#' collision is why the diffusion parameters carry `frmtmb.eam`'s names
#' rather than hBayesDM's: one convention across this package's families
#' is worth more than an inherited one that would make `alpha` mean two
#' different things two families apart.
#'
#' @section The data a trial carries:
#' The response is the response TIME. The boundary reached travels in
#' `dec()`, coded 0 for the lower boundary and 1 for the upper, which is
#' the term and the coding [frmtmb.eam::wiener()] uses. Arm 1 is the
#' lower boundary and arm 2 the upper, so `reward(pay1, pay2)` is in the
#' same arm order as every other family here.
#'
#' ```
#' frm(bf(rt | dec(choice) + reward(pay1, pay2) ~ 1 + (1 | id),
#'        drift ~ 1, bs ~ 1, ndt ~ 1, bias ~ 1),
#'     family = rlddm(subject = id, trial = trial), data = d)
#' ```
#'
#' @section Holding the start point at a half:
#' Most of the reinforcement-learning literature fixes `bias` at 0.5,
#' because with two ARMS rather than a correct and an error response
#' there is no reason for the accumulator to start nearer one boundary.
#' It is an ordinary distributional parameter here, so it is estimated
#' by default and held with `bf()`'s constant form:
#'
#' ```
#' frm(bf(rt | dec(choice) + reward(pay1, pay2) ~ 1 + (1 | id),
#'        drift ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5), ...)
#' ```
#'
#' Estimating it is the more general model and costs one parameter;
#' fixing it is what makes `drift` and `bs` easier to separate on a
#' short session.
#'
#' @section The non-decision time is bounded, not logged:
#' The density is zero at and below `ndt`, so the likelihood has a hard
#' edge at the fastest response and a log link would let the optimizer
#' walk over it. `ndt` therefore gets a logit scaled onto `(0, ub)`,
#' which makes the constraint structural. The bound comes from
#' [frmtmb.eam::ndt_bound()], so this family and [frmtmb.eam::wiener()]
#' derive it the same way and refuse the same things.
#'
#' @section One bound, or one per learner:
#' `ub` is one number for the whole data set unless the model says
#' otherwise, and the whole data set's fastest response is the wrong
#' ceiling for a hierarchical fit. The information about one subject's
#' non-decision time is that SUBJECT's own fastest response, and with
#' 100 learners the global minimum is the fastest of all of them, so a
#' subject deviation on `ndt` is a deviation on a fraction of somebody
#' else's floor. The scale tier measured what that costs: the 100 by 200
#' design reached a maximum gradient of 6.36e+09, a Hessian that was not
#' positive definite and four `NaN` standard errors.
#'
#' `ndt_group()` is the fix, and it is opt-in:
#'
#' ```
#' frm(bf(rt | dec(choice) + reward(pay1, pay2) + ndt_group(id) ~
#'          1 + (1 | p | id),
#'        drift ~ 1 + (1 | p | id), bs ~ 1 + (1 | p | id),
#'        ndt ~ 1 + (1 | p | id), bias = 0.5),
#'     family = rlddm(subject = id, trial = trial), data = d)
#' ```
#'
#' With it, each row's bound is its own group's fastest response, `ndt`
#' is a FRACTION of that bound on a plain logit, and the density
#' multiplies it back out. `predict(dpar = "ndt", type = "response")`
#' then reports the fraction, and [frmtmb.eam::ndt_time()] reports the
#' non-decision time in seconds whichever parameterization a fit is in.
#' Without it nothing about the family moved.
#'
#' The grouping need not be the learner: it is whatever shares a floor.
#' `max_ndt` and `ndt_group()` together are refused, because they set
#' one bound to two different things.
#'
#' @section What the per-learner spread of ndt does and does not say:
#' Under `ndt_group()` a learner's non-decision time is a FRACTION times
#' that learner's own floor, so the fitted times differ between learners
#' even when the random effect on `ndt` is exactly zero, purely because
#' the floors differ. A spread across learners is therefore not evidence
#' that the variance component recovered anything.
#'
#' Measured on the 100 by 200 design of `dev/rlddm-findings.md`, where
#' the drawn truths have a spread of 0.0393041: the full model returns
#' **0.0333713** and the same model with the random effect on `ndt`
#' switched OFF returns **0.0380049**. On this statistic the random
#' effect is not neutral, it is worse.
#'
#' What the random effect IS worth on the same fit is the per-learner
#' accuracy, and there it is not close: 13.996 ms of error against the
#' drawn truths and a correlation of 0.9404, where no estimator without
#' it can beat 21.075 ms or 0.8478, because the correlation of any
#' constant multiple of the floors with the truth is fixed at 0.847780.
#' So read the per-learner error and the correlation, not the spread.
#' [frmtmb.eam::wiener()] carries the same warning for the same reason.
#'
#' @section Pinning the non-decision time needs max_ndt:
#' `bf(ndt = 0.2)` on a bare `rlddm()` is refused, and the message says
#' to pass `max_ndt`. `bf()` transforms a pinned constant at PARSE time,
#' before the response has been seen, so a family with no bound yet has
#' no scale to transform it on. Write the bound down instead:
#'
#' ```
#' frm(bf(rt | dec(choice) + reward(pay1, pay2) ~ 1,
#'        drift ~ 1, bs ~ 1, ndt = 0.2, bias = 0.5),
#'     family = rlddm(subject = id, trial = trial, max_ndt = 0.26),
#'     data = d)
#' ```
#'
#' **This is a change, and it is one thing that stops working.** Through
#' 0.3.0 `bf(ndt = 0.2)` on a bare family was accepted and fitted 0.2
#' exactly, so far as the constant was below the fastest response. What
#' it did NOT do is check the constant against the bound the fit would
#' actually use: `bf(ndt = 0.3)` with `max_ndt = 0.26` was accepted too,
#' and reached the objective as `NaN`, surfacing as an optimizer failure
#' that named nothing. Both cases are now settled before the formula is
#' parsed, which is how [frmtmb.eam::wiener()] and its siblings have
#' behaved since frmtmb.eam 0.7.0.
#'
#' @section Recovery, measured:
#' `dev/learn-recovery-round2.R` runs it and `?frmtmb.learn` carries the
#' table. All five parameters recover at 30 subjects by 100 trials with
#' a random intercept on the learning rate: every bias is inside its
#' Monte Carlo error, and 60 of 60 replicates produced a usable
#' interval. Coverage is at the nominal rate for `drift`, `bs`, `ndt`
#' and `bias`, and 0.87 for `alpha`, so a learning rate's Wald interval
#' is slightly optimistic here and should be read as such.
#'
#' THE PAIR THAT TRADES OFF IS NOT THE ONE TO EXPECT, and this section
#' said otherwise before the study was run. A first-passage density is
#' driven largely by the RATIO of the drift to the boundary, so `drift`
#' and `bs` look like the pair at risk; measured, their estimates
#' correlate at 0.17 across replicates, which is nothing. What does
#' co-vary is `bs` with `ndt`, at -0.57, and `drift` with `bias`, at
#' +0.56. Both make sense after the fact: the boundary separation and
#' the non-decision time are the two ways to make responses slower, and
#' the drift and the start point are the two ways to favor a boundary.
#' Read those two pairs together.
#'
#' @section Drawing data needs RWiener, and fitting does not:
#' [frm_task_simulate()] on this family draws a boundary and a time
#' JOINTLY from the diffusion, through
#' [frmtmb.eam::ddm_simulate()], and that route is exact only when
#' `RWiener` is installed. Without it the draw is refused rather than
#' approximated: `frmtmb.eam`'s own fallback draws the boundary first
#' and then rejects paths until one reaches it, which fails outright on
#' the rows where one boundary is strongly favored, and a simulator
#' that works on most rows and errors on the rest is worse than one
#' that says what it needs. Nothing about FITTING needs `RWiener`; the
#' likelihood is [frmtmb.eam::wiener_lpdf()] and is exact either way.
#'
#' @section Where the density and the bound come from:
#' [frmtmb.eam::wiener_lpdf()] and [frmtmb.eam::ndt_bound()], through
#' the exports rather than through a colon. The Navarro and Fuss (2009)
#' series pair, its smooth blend and its tape safety belong to that
#' package, and so does the bound, its `ndt_group()` behavior and every
#' refusal either of them owes. This family writes the learning rule and
#' one multiplication.
#'
#' @param subject The column separating one learner's trial sequence
#'   from the next, given unquoted.
#' @param trial The column giving trial order within a subject, given
#'   unquoted. `NULL` uses the order the rows appear in.
#' @param max_ndt Upper bound for the non-decision time. `NULL`, the
#'   default, uses the fastest response in the data. A value above it is
#'   refused, because it admits parameters at which the fastest trial
#'   has no likelihood. It cannot be combined with
#'   `ndt_group()`, which sets the same bound per group. It is also what
#'   `bf(ndt = )` now needs: see below.
#'
#' @return A `frmtmb_family` object, for `frm(family = )`.
#'
#' @references
#' Pedersen, M. L., Frank, M. J. and Biele, G. (2017). The drift
#' diffusion model as the choice rule in reinforcement learning.
#' *Psychonomic Bulletin and Review* 24, 1234-1251.
#'
#' Navarro, D. J. and Fuss, I. G. (2009). Fast and accurate calculations
#' for first-passage times in Wiener diffusion models. *Journal of
#' Mathematical Psychology* 53, 222-230.
#'
#' @seealso [bandit2arm_delta()] for the same learning rule under a
#'   softmax, [frmtmb.eam::wiener()] for the diffusion without learning,
#'   [frm_value_trace()] for the per-trial drift rates.
#'
#' @examples
#' d <- frm_task_design("bandit2arm", n_subject = 6, n_trial = 40,
#'                      seed = 1)
#' s <- frm_task_simulate(rlddm(subject = id, trial = trial), d,
#'                        pars = list(alpha = 0.4, drift = 3, bs = 1.6,
#'                                    ndt = 0.2, bias = 0.5), seed = 1)
#' fit <- frmtmb::frm(
#'   frmtmb::bf(rt | dec(choice) + reward(pay1, pay2) ~ 1,
#'              drift ~ 1, bs ~ 1, ndt ~ 1, bias ~ 1),
#'   family = rlddm(subject = id, trial = trial), data = s[[1]])
#' frmtmb::fixef(fit)
#' head(frm_value_trace(fit))
#' @export
rlddm <- function(subject, trial = NULL, max_ndt = NULL) {
  spec <- ln_spec(
    n_option = 2L,
    init = function(ns, d1) list(q1 = rep(0, ns), q2 = rep(0, ns)),
    # no softmax: the trial's contribution is a density over the
    # response time, and the option taken selects a boundary rather
    # than one term of a normalizing sum
    choice = NULL,
    choice_col = "dec",
    # dec() codes the lower boundary 0 and the upper 1; the engine
    # counts options from 1, and arm 2 is the upper boundary
    choice_map = function(v) v + 1,
    logp = function(state, d, ch) {
      v <- d[["drift"]] * (state[["q2"]] - state[["q1"]])
      wiener_lpdf(d[[".y"]] - ln_ndt_at(d), v, d[["bs"]], d[["bias"]],
                  ch[[2L]])
    },
    # The draw reads the non-decision time the same way the density
    # does, and today that is a no-op on the only path that reaches it:
    # frm_task_simulate() builds its per-row values out of `pars` and
    # the data map, so no ndt_floor and no ndt_group is ever in them and
    # `ndt` arrives as a time already. It goes through ln_ndt_at()
    # anyway, because the two halves of a family reading one parameter
    # two different ways is how a fraction gets drawn as seconds the day
    # a simulator is turned on here.
    draw = function(state, d) {
      ns <- length(state[["q1"]])
      v <- as.numeric(d[["drift"]] * (state[["q2"]] - state[["q1"]]))
      sim <- ddm_simulate(
        ns, v, as.numeric(rep_len(d[["bs"]], ns)),
        as.numeric(rep_len(ln_ndt_at(d), ns)),
        as.numeric(rep_len(d[["bias"]], ns)))
      d[["dec"]] <- sim[["upper"]]
      list(y = sim[["rt"]], choice = sim[["upper"]] + 1, d = d)
    },
    update = function(state, d, ch) {
      q1 <- state[["q1"]]
      q2 <- state[["q2"]]
      c1 <- ch[[1L]][[1L]]
      c2 <- ch[[1L]][[2L]]
      pe1 <- d[["reward1"]] - q1
      pe2 <- d[["reward2"]] - q2
      list(state = list(q1 = q1 + c1 * d[["alpha"]] * pe1,
                        q2 = q2 + c2 * d[["alpha"]] * pe2),
           drift_t = d[["drift"]] * (q2 - q1),
           pe = c1 * pe1 + c2 * pe2)
    },
    sim_cols = "dec")
  fam <- ln_family("rlddm", substitute(subject), substitute(trial),
            dpars = c("alpha", "drift", "bs", "ndt", "bias"),
            links = list(alpha = "logit", drift = "identity",
                         bs = "log", ndt = "log", bias = "logit"),
            primary = "alpha",
            inits = list(alpha = function(y, aterms) 0.3,
                         drift = function(y, aterms) 1,
                         bs = function(y, aterms) 1.5,
                         ndt = function(y, aterms) 0.5 * min(y),
                         bias = function(y, aterms) 0.5),
            aterms = c("dec", "reward1", "reward2"), spec = spec,
            data_map = c(reward1 = "pay1", reward2 = "pay2",
                         dec = "choice"),
            sim = FALSE,
            constants = list(sim_out = list(dec = "choice"),
                             sim_response = "rt",
                             # frm_task_simulate() checks this once per
                             # call, before the walk starts. See
                             # ?rlddm's "Drawing data" section for why
                             # the draw needs a package the likelihood
                             # does not.
                             sim_needs = "RWiener"),
            valid_y = ln_valid_rt,
            # A row's contribution is a DENSITY, and a density has no
            # saturated value: the supremum of a Wiener log density at a
            # fixed response time is unbounded over the parameters, so
            # there is no constant to compare a fitted row against. The
            # six softmax families do have one, and it is zero.
            saturated = FALSE,
            finalize = function(fam, y, aterms) {
              # The refusals come FIRST, before the carried bound is
              # read, so that a max_ndt above the fastest response is
              # still refused on the second finalize of a family whose
              # bound is already settled. frmtmb.eam's own families put
              # them in the same order for the same reason.
              bd <- ndt_bound(y, aterms, max_ndt, "rlddm")
              if (!is.null(ndt_bound_of(fam))) return(fam)
              ndt_bound_attach(fam, bd)
            },
            sim_refusal = paste0(
              "One rlddm() trial's draw is two numbers, the boundary ",
              "reached and the time it took, and a response vector ",
              "holds one. Drawing the time alone would keep the ",
              "observed choice against a freshly drawn response time, ",
              "which is a draw from a different model. Use ",
              "frmtmb.learn::frm_task_simulate(), which returns whole ",
              "data frames"))
  # THE BOUND BEFORE THERE IS DATA. `ndt`'s declared link above is a
  # placeholder, and through 0.3.0 it was what a family object carried
  # until frm() replaced it.
  #
  # WHAT THAT DID AND DID NOT COST, measured, because a first version of
  # this comment claimed the wrong one. frmtmb transforms a constant
  # dpar twice: parse.R only RANGE-CHECKS it against whatever link the
  # family is carrying, and frame.R does the transform that reaches the
  # parameter, after family_finalize() has settled the link. So an
  # in-range `bf(ndt = 0.2)` was fitted at exactly 0.2 and nothing was
  # wrong with it. What was wrong is the range check: at 0.3.0 it ran
  # against the placeholder `log`, which accepts any positive constant,
  # so `bf(ndt = 0.3)` with max_ndt = 0.26 passed it, produced NaN in
  # frame.R with no finiteness check, and surfaced as "NA/NaN gradient
  # evaluation" from the optimizer rather than as anything naming the
  # constant.
  #
  # A pending bound makes the check honest, at the price of `bf(ndt = )`
  # needing `max_ndt` on a bare family. That price is the one
  # wiener(), lba(), rdm() and wiener_gng() have paid since
  # frmtmb.eam 0.7.0, and paying it here is what makes one contract
  # across the five families rather than four and an exception.
  ndt_bound_attach(fam, ndt_bound_pending(max_ndt, "rlddm"))
}

#' The response of a joint choice-and-time family.
#'
#' @noRd
ln_valid_rt <- function(y, aterms) {
  if (!is.numeric(y) || anyNA(y) || any(y <= 0)) {
    stop("rlddm(): the response is the response TIME, a positive ",
         "number, and the option taken travels in dec(). A response ",
         "coded 1 or 2 is the usual mistake: that is what the six ",
         "softmax families here take, and this one takes the time",
         call. = FALSE)
  }
  invisible(NULL)
}

#' The non-decision time at one trial's rows, as a TIME.
#'
#' The whole of what a per-group bound costs this likelihood, and it is
#' one call. With `ndt_group()` the `ndt` link is a plain logit on a
#' FRACTION of the row's own bound, and that bound arrives beside the
#' addition-term values as `ndt_floor`, one entry per row, so the
#' recursion indexes it at each trial exactly as it indexes the payoffs.
#' Without the grouping there is no `ndt_floor`, `ndt` is already a
#' time, and it comes back untouched: a model that does not opt in is
#' the arithmetic 0.3.0 shipped.
#'
#' THE ARITHMETIC IS NOT THE PART THAT MATTERS. A grouped model whose
#' density is reached without `ndt_floor` would read the fraction as
#' seconds and report a converged fit at a non-decision time several
#' times too small. That refusal is longer than the multiplication it
#' guards, which is why it is [frmtmb.eam::ndt_apply()]'s and not a
#' copy here: the review of this item found that a consumer copying the
#' documented arithmetic and not the undocumented refusal gets exactly
#' the defect the item exists to remove.
#'
#' The engine merges the distributional parameters and the addition-term
#' values into one list, so `d` is passed once and serves as both.
#' `ndt_apply()` is imported by name rather than reached with `::`,
#' because this runs once per trial per objective evaluation and a
#' namespace lookup does not belong in that loop.
#'
#' @noRd
ln_ndt_at <- function(d) ndt_apply(d, d, "rlddm()")
