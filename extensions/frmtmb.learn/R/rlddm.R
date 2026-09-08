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
#' edge at `ndt = min(rt)` and a log link would let the optimizer walk
#' over it. `ndt` therefore gets a logit scaled onto `(0, max_ndt)`,
#' with `max_ndt` defaulting to the fastest response in the data, which
#' makes the constraint structural. This is the same construction
#' [frmtmb.eam::wiener()] uses and it is written again here rather than
#' borrowed, because a link is not what that package exports.
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
#' @section Where the density comes from:
#' [frmtmb.eam::wiener_lpdf()], through the export rather than through a
#' colon. The Navarro and Fuss (2009) series pair, its smooth blend and
#' its tape safety all belong to that package and are not repeated here.
#' This is the only dependency one extension of frmtmb has on another,
#' and it is one function.
#'
#' @param subject The column separating one learner's trial sequence
#'   from the next, given unquoted.
#' @param trial The column giving trial order within a subject, given
#'   unquoted. `NULL` uses the order the rows appear in.
#' @param max_ndt Upper bound for the non-decision time. `NULL`, the
#'   default, uses the fastest response in the data. A value above it is
#'   refused, because it admits parameters at which the fastest trial
#'   has no likelihood.
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
      wiener_lpdf(d[[".y"]] - d[["ndt"]], v, d[["bs"]], d[["bias"]],
                  ch[[2L]])
    },
    draw = function(state, d) {
      ns <- length(state[["q1"]])
      v <- as.numeric(d[["drift"]] * (state[["q2"]] - state[["q1"]]))
      sim <- ddm_simulate(
        ns, v, as.numeric(rep_len(d[["bs"]], ns)),
        as.numeric(rep_len(d[["ndt"]], ns)),
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
  ln_family("rlddm", substitute(subject), substitute(trial),
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
            finalize = function(fam, y, aterms) ln_ndt_link(fam, y, max_ndt),
            sim_refusal = paste0(
              "One rlddm() trial's draw is two numbers, the boundary ",
              "reached and the time it took, and a response vector ",
              "holds one. Drawing the time alone would keep the ",
              "observed choice against a freshly drawn response time, ",
              "which is a draw from a different model. Use ",
              "frmtmb.learn::frm_task_simulate(), which returns whole ",
              "data frames"))
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

#' A logit scaled onto `(0, ub)` for the non-decision time.
#'
#' WHY IT IS WRITTEN HERE. The density is zero at and below `ndt`, so
#' the likelihood has a hard edge at `ndt = min(rt)`; a log link lets an
#' optimizer step over it into a region where the fastest trial has no
#' likelihood at all. Making the bound part of the link makes the
#' constraint structural instead.
#'
#' `frmtmb.eam` builds the same link for the same reason, and this is
#' the one thing this family reproduces rather than imports. What that
#' package exports is a DENSITY; a link is not part of the promise it
#' makes, and asking for a second export whose subject is link
#' construction would widen its public surface for ten lines of
#' arithmetic. The duplication is recorded in dev/learn2-findings.md
#' rather than hidden.
#'
#' @noRd
ln_ndt_link <- function(fam, y, max_ndt) {
  ub <- if (is.null(max_ndt)) min(y) else max_ndt
  if (!is.null(max_ndt) && ub > min(y)) {
    stop("rlddm(max_ndt = ", format(ub), ") is above the fastest ",
         "response (", format(min(y)), "). Nothing is observed before ",
         "the non-decision time, so a bound above the fastest response ",
         "admits values at which that trial has no likelihood",
         call. = FALSE)
  }
  fam[["links"]][["ndt"]] <- list(
    name = paste0("scaled_logit(0, ", signif(ub, 4), ")"),
    linkfun = function(mu) log(mu / (ub - mu)),
    linkinv = function(eta) ub / (1 + exp(-eta)),
    mu_eta = function(eta) {
      p <- 1 / (1 + exp(-eta))
      ub * p * (1 - p)
    })
  fam
}
