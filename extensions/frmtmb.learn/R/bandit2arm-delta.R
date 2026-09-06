#' Rescorla-Wagner delta learning on a two-armed bandit
#'
#' The simplest value-learning model, and the one every other family
#' here is a variation on. A subject keeps a value estimate `Q` for each
#' arm, starts both at zero, and after each choice moves the chosen
#' arm's estimate toward what it just paid:
#'
#' ```
#' Q[chosen] <- Q[chosen] + alpha * (reward - Q[chosen])
#' ```
#'
#' The choice is a softmax over the two values, which for two arms is a
#' logistic function of their difference:
#'
#' ```
#' P(arm 1) <- plogis(tau * (Q1 - Q2))
#' ```
#'
#' This is the model hBayesDM calls `bandit2arm_delta`. Its two
#' parameters are named `A` and `tau` there and `alpha` and `tau` here;
#' see [frm_learn_families()] for the whole naming map.
#'
#' @section The data a trial carries:
#' The response is the arm chosen, coded 1 or 2. `reward(pay1, pay2)`
#' carries what each arm WOULD have paid on this trial, in arm order.
#' The likelihood only ever reads the chosen arm's entry, so data that
#' records the received outcome alone can pass it twice. The second
#' column is what makes the simulator coherent: a simulated choice needs
#' the payoff of the arm the subject did not take in the data.
#'
#' @section Reversal learning, and what the grammar buys:
#' Nothing here is specific to a stationary bandit. A probabilistic
#' reversal-learning task is this family with a covariate on the
#' learning rate, because `alpha` is an ordinary distributional
#' parameter with its own linear predictor:
#'
#' ```
#' frm(bf(choice | reward(pay1, pay2) ~ after_reversal + (1 | id),
#'        tau ~ 1 + (1 | id)),
#'     family = bandit2arm_delta(subject = id, trial = trial), data = d)
#' ```
#'
#' A separate `prl` family would fit one learning rate before the
#' reversal and one after. This fits the DIFFERENCE, with a standard
#' error, and it takes a factor with any number of levels, a smooth term
#' or a random slope in the same place. `vignette("learning")` works the
#' example through.
#'
#' @section The Laplace caveat:
#' frmtmb integrates the subject effects out with a Laplace
#' approximation, which is exact only when the conditional log-density
#' is quadratic. For binary choices it is not, and the fewer trials a
#' subject has the less quadratic it is. Measured on this family, the
#' fixed effects survive short sessions and the variance components do
#' not: see the Laplace section of `vignette("learning")` and
#' `?frmtmb.learn` for the numbers. `frm(importance =)`, which is the
#' usual way to price that error, is REFUSED for every family here; the
#' refusal names the seam.
#'
#' @param subject The column separating one learner's trial sequence
#'   from the next, given unquoted. It is carried by the family rather
#'   than by the formula, the way `frmtmb.latent::hmm()` carries its
#'   sequence grouping, because it orders the likelihood rather than
#'   entering a linear predictor.
#' @param trial The column giving trial order within a subject, given
#'   unquoted. `NULL` uses the order the rows appear in.
#'
#' @return A `frmtmb_family` object, for `frm(family = )`.
#'
#' @references
#' Ahn, W.-Y., Haines, N. and Zhang, L. (2017). Revealing
#' neurocomputational mechanisms of reinforcement learning and
#' decision-making with the hBayesDM package. *Computational
#' Psychiatry* 1, 24-57.
#'
#' @seealso [bandit2arm_dual()] for separate rates on gains and losses,
#'   [prl_fictitious()] for counterfactual updating,
#'   [frm_value_trace()] for the fitted values and prediction errors.
#'
#' @examples
#' d <- frm_task_design("bandit2arm", n_subject = 8, n_trial = 40,
#'                      seed = 1)
#' d$choice <- frm_task_simulate(
#'   bandit2arm_delta(subject = id, trial = trial), d,
#'   pars = list(alpha = 0.4, tau = 3), seed = 1)[[1]]$choice
#' fit <- frmtmb::frm(frmtmb::bf(choice | reward(pay1, pay2) ~ 1,
#'                               tau ~ 1),
#'                    family = bandit2arm_delta(subject = id,
#'                                              trial = trial), data = d)
#' frmtmb::fixef(fit)
#' head(frm_value_trace(fit))
#' @export
bandit2arm_delta <- function(subject, trial = NULL) {
  spec <- ln_spec(
    n_option = 2L,
    init = function(ns, d1) list(q1 = rep(0, ns), q2 = rep(0, ns)),
    choice = function(state, d, j) {
      list(d[["tau"]] * state[["q1"]], d[["tau"]] * state[["q2"]])
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
           pe = c1 * pe1 + c2 * pe2)
    })
  ln_family("bandit2arm_delta", substitute(subject), substitute(trial),
            dpars = c("alpha", "tau"),
            links = list(alpha = "logit", tau = "log"),
            primary = "alpha",
            inits = list(alpha = function(y, aterms) 0.3,
                         tau = function(y, aterms) 1),
            aterms = c("reward1", "reward2"), spec = spec,
            data_map = c(reward1 = "pay1", reward2 = "pay2"))
}
