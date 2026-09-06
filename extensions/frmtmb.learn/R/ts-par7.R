#' The two-step task: a model-based and model-free hybrid
#'
#' The task of Daw and others (2011), and the model that made it famous.
#' A trial has two choices. The first leads, usually but not always, to
#' one of two second-stage states; the second is a choice within that
#' state, and it pays. A purely model-free learner repeats what was
#' rewarded. A model-based learner knows the transition structure and
#' works backward through it, so a reward that arrived after a RARE
#' transition makes it switch rather than stay. `w` is how much of each.
#'
#' ```
#' QMB[a] <- P(common) * max(Q2[state a leads to]) +
#'           P(rare)   * max(Q2[the other state])
#' Q1[a]  <- w * QMB[a] + (1 - w) * QMF[a]
#' P(a)   <- softmax(tau1 * (Q1 + pers * repeated))
#' # after the outcome
#' delta1 <- Q2[s2, a2] - QMF[a1];  QMF[a1] <- QMF[a1] + alpha1 * delta1
#' delta2 <- reward   - Q2[s2, a2]; Q2[s2, a2] <- Q2[s2, a2] + alpha2 * delta2
#' QMF[a1] <- QMF[a1] + alpha1 * lambda * delta2
#' ```
#'
#' `lambda` is the eligibility trace: how much of the second-stage
#' prediction error reaches the first-stage value directly. `pers` is
#' perseveration, a constant bonus for repeating the previous trial's
#' first-stage choice. hBayesDM calls the seven-parameter version
#' `ts_par7` and spells them `a1`, `beta1`, `a2`, `beta2`, `lambda`,
#' `w` and `pi`.
#'
#' @section Two ways this differs from hBayesDM's `ts_par7`:
#' Five of the seven parameters are a plain rename. Two things are not,
#' and both are read off hBayesDM 2.0.0's published Stan source.
#'
#' **`pers` does not port, because perseveration sits on the other side
#' of the temperature.** hBayesDM adds it OUTSIDE `beta1`:
#' `inv_logit(beta1 * (v2 - v1) + pi * (...))`. This family adds it
#' INSIDE the bracket `tau1` multiplies, so its contribution to the same
#' logit is `tau1 * pers * (...)`. Matching term by term,
#'
#' ```
#' pers = pi_hBayesDM / tau1
#' ```
#'
#' `tau1` is estimated, so as with [prl_fictitious()]'s `bias` there is
#' no constant transform: convert through a fit's own `tau1`, or refit.
#'
#' **The eligibility trace uses a different prediction error, and this
#' family follows the paper.** The trace term here is formed from the
#' PRE-update `delta2`, which is the equation Daw and others (2011)
#' print. hBayesDM's code updates the stage-two value first and then
#' forms the trace from the already-updated value, so its effective
#' trace is `lambda * a1 * (1 - a2) * delta2`. The two agree only at
#' `a2 = 0`. This is a difference in the reference implementation rather
#' than an error here, and it is recorded because a user comparing fits
#' will otherwise meet it with no stated cause.
#'
#' @section The data a trial carries:
#' The response is the FIRST-stage choice, 1 or 2. `stage2(state,
#' choice)` carries the second-stage state (1 or 2, for the two states)
#' and the choice made in it (1 or 2). `payoff(p1, p2, p3, p4)` carries
#' what each of the four second-stage options would have paid, indexed
#' `2 * (state - 1) + choice`.
#'
#' The transition probability is a known constant of the TASK rather
#' than a parameter, so `P(state | first choice)` contributes a constant
#' to the log likelihood and is dropped. It still enters the model-based
#' value, which is what `p_common` is for.
#'
#' @section Why `simulate()` is refused:
#' One trial's draw is three numbers: the first-stage choice, the state
#' the environment answers with, and the second-stage choice.
#' [stats::simulate()] and [frmtmb::frm_simulate()] return one response
#' vector, and there is nowhere to put the other two, so both refuse by
#' name. [frm_task_simulate()] returns whole data frames and does the
#' draw properly; it is also what every other family in this package
#' uses for a recovery study.
#'
#' @inheritSection bandit2arm_delta The Laplace caveat
#' @inheritParams bandit2arm_delta
#' @param p_common The probability that a first-stage choice leads to
#'   its common second-stage state. 0.7 in the original task. Known,
#'   not estimated: the subject's BELIEF about it is not separable from
#'   `w` on the data these studies collect.
#'
#' @return A `frmtmb_family` object, for `frm(family = )`.
#'
#' @references
#' Daw, N. D., Gershman, S. J., Seymour, B., Dayan, P. and Dolan, R. J.
#' (2011). Model-based influences on humans' choices and striatal
#' prediction errors. *Neuron* 69, 1204-1215.
#'
#' @seealso [frm_task_design()], [frm_task_simulate()]
#'
#' @examples
#' d <- frm_task_design("twostep", n_subject = 5, n_trial = 40, seed = 8)
#' d <- frm_task_simulate(
#'   ts_par7(subject = id, trial = trial), d,
#'   pars = list(alpha1 = 0.4, tau1 = 3, alpha2 = 0.4, tau2 = 3,
#'               lambda = 0.6, w = 0.5, pers = 0.2), seed = 8)[[1]]
#' head(d[, c("id", "trial", "choice", "state2", "choice2")])
#' @export
ts_par7 <- function(subject, trial = NULL, p_common = 0.7) {
  if (!is.numeric(p_common) || length(p_common) != 1L || is.na(p_common) ||
        p_common <= 0.5 || p_common >= 1) {
    stop("ts_par7(p_common =) is the probability of the COMMON ",
         "transition and must be one number strictly between 0.5 and 1; ",
         "at 0.5 the two states are indistinguishable and there is no ",
         "model-based value to compute", call. = FALSE)
  }
  pc <- p_common
  q2_of <- function(state, s2ind, a) {
    # the stage-two value of action `a` in the OBSERVED state, selected
    # with a data indicator rather than a branch
    s2ind * state[[paste0("q2_", a)]] +
      (1 - s2ind) * state[[paste0("q2_", a + 2L)]]
  }
  spec <- ln_spec(
    n_option = c(2L, 2L),
    resp = list(NULL, "stage22"),
    init = function(ns, d1) {
      z <- rep(0, ns)
      list(qmf1 = z, qmf2 = z, q2_1 = z, q2_2 = z, q2_3 = z, q2_4 = z,
           rep1 = z, rep2 = z)
    },
    choice = function(state, d, j) {
      if (j == 1L) {
        # model-based: work backward through the known transition
        best1 <- ln_max2(state[["q2_1"]], state[["q2_2"]])
        best2 <- ln_max2(state[["q2_3"]], state[["q2_4"]])
        qmb1 <- pc * best1 + (1 - pc) * best2
        qmb2 <- (1 - pc) * best1 + pc * best2
        w <- d[["w"]]
        p <- d[["pers"]]
        u1 <- w * qmb1 + (1 - w) * state[["qmf1"]] + p * state[["rep1"]]
        u2 <- w * qmb2 + (1 - w) * state[["qmf2"]] + p * state[["rep2"]]
        list(d[["tau1"]] * u1, d[["tau1"]] * u2)
      } else {
        s2ind <- as.numeric(d[["stage21"]] == 1)
        list(d[["tau2"]] * q2_of(state, s2ind, 1L),
             d[["tau2"]] * q2_of(state, s2ind, 2L))
      }
    },
    update = function(state, d, ch) {
      s2ind <- as.numeric(d[["stage21"]] == 1)
      c1 <- ch[[1L]][[1L]]
      c2 <- ch[[1L]][[2L]]
      k1 <- ch[[2L]][[1L]]
      k2 <- ch[[2L]][[2L]]
      # a one-of-four indicator over the stage-two options, all data
      sel <- list(s2ind * k1, s2ind * k2, (1 - s2ind) * k1,
                  (1 - s2ind) * k2)
      r <- 0
      qs2 <- 0
      for (k in 1:4) {
        r <- r + sel[[k]] * d[[paste0("payoff", k)]]
        qs2 <- qs2 + sel[[k]] * state[[paste0("q2_", k)]]
      }
      qmf <- c1 * state[["qmf1"]] + c2 * state[["qmf2"]]
      delta1 <- qs2 - qmf
      delta2 <- r - qs2
      # the eligibility trace: the stage-two prediction error reaches
      # the stage-one value directly, discounted by lambda
      step1 <- d[["alpha1"]] * (delta1 + d[["lambda"]] * delta2)
      step2 <- d[["alpha2"]] * delta2
      out <- list(qmf1 = state[["qmf1"]] + c1 * step1,
                  qmf2 = state[["qmf2"]] + c2 * step1)
      for (k in 1:4) {
        out[[paste0("q2_", k)]] <- state[[paste0("q2_", k)]] +
          sel[[k]] * step2
      }
      out[["rep1"]] <- c1
      out[["rep2"]] <- c2
      list(state = out, pe = delta1, pe2 = delta2, qmf = qmf)
    },
    draw_env = function(d, ch, j) {
      if (j == 1L) {
        # the environment answers the first choice with a state: the
        # common one with probability p_common
        a1 <- ch[[1L]][[2L]] + 1
        common <- stats::runif(length(a1)) < pc
        d[["stage21"]] <- ifelse(common, a1, 3 - a1)
      }
      d
    },
    sim_cols = c("stage21", "stage22"))
  ln_family("ts_par7", substitute(subject), substitute(trial),
            dpars = c("w", "alpha1", "tau1", "alpha2", "tau2", "lambda",
                      "pers"),
            links = list(w = "logit", alpha1 = "logit", tau1 = "log",
                         alpha2 = "logit", tau2 = "log",
                         lambda = "logit", pers = "identity"),
            # w is primary, so the main right-hand side reaches the
            # model-based weight. It is the quantity the two-step
            # literature compares between groups.
            primary = "w",
            inits = list(w = function(y, aterms) 0.5,
                         alpha1 = function(y, aterms) 0.3,
                         tau1 = function(y, aterms) 1,
                         alpha2 = function(y, aterms) 0.3,
                         tau2 = function(y, aterms) 1,
                         lambda = function(y, aterms) 0.5,
                         pers = function(y, aterms) 0),
            aterms = c("stage21", "stage22", paste0("payoff", 1:4)),
            spec = spec, sim = FALSE,
            data_map = c(stage21 = "state2", stage22 = "choice2",
                         stats::setNames(paste0("pay", 1:4),
                                         paste0("payoff", 1:4))),
            constants = list(p_common = p_common,
                             sim_out = list(stage21 = "state2",
                                            stage22 = "choice2")),
            sim_refusal = paste0(
              "One two-step trial's draw is three numbers, the ",
              "stage-one choice, the stage-two state the environment ",
              "answers with and the stage-two choice, and a response ",
              "vector holds one. Use frmtmb.learn::frm_task_simulate(), ",
              "which returns whole data frames"))
}
