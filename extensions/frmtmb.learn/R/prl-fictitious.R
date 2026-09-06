#' Counterfactual (fictitious) updating on two options
#'
#' The unchosen option is updated too, in the opposite direction. In a
#' two-armed task where one option pays when the other does not, that is
#' what a subject who understands the task should do, and it is the
#' model of Glascher, Hampton and O'Doherty (2009) for probabilistic
#' reversal learning. hBayesDM calls it `prl_fictitious`.
#'
#' ```
#' EV[chosen]   <- EV[chosen]   + alpha * ( outcome - EV[chosen])
#' EV[unchosen] <- EV[unchosen] + alpha * (-outcome - EV[unchosen])
#' P(option 1)  <- plogis(tau * (EV1 - EV2) + bias)
#' ```
#'
#' The counterfactual outcome is the NEGATIVE of the realized one,
#' which is the model's assumption of an anticorrelated task rather
#' than a reading of the second `reward()` column. The second column is
#' still what the simulator pays a counterfactual choice with.
#'
#' @section What `bias` is, and why it does not port:
#' `bias` is the indecision point: a constant added to the utility
#' difference, so a subject with `bias` above zero prefers option 1 at
#' equal value. It has an identity link.
#'
#' **It is not hBayesDM's `alpha` renamed, and there is no fixed
#' transform between them.** hBayesDM's `prl_fictitious` writes the
#' choice probability as `inv_logit(beta * (alpha - (ev1 - ev2)))`,
#' which is `inv_logit(beta * (ev1 - ev2) - beta * alpha)` with the
#' option labels swapped. This family writes
#' `inv_logit(tau * (ev1 - ev2) + bias)`. Matching the two term by term
#' gives
#'
#' ```
#' bias = -tau * alpha_hBayesDM
#' ```
#'
#' The factor is `tau`, which is ESTIMATED, so the map depends on the
#' fit and no constant relates the two parameters. At the `tau = 3` this
#' page's example uses, an indecision point carried across unchanged is
#' wrong by a factor of three, and the resulting model converges quietly.
#' Convert through the equation above, using the same fit's `tau`, or
#' refit.
#'
#' What is implemented here is the equation printed at the top of this
#' page, and it is that equation the Stan identity in
#' `tests/testthat/test-stan-identity.R` checks to 8.5e-14. The relation
#' above is read off hBayesDM 2.0.0's published Stan source rather than
#' measured against a fit: `dev/hbayesdm-crosscheck.R` would measure it
#' and has not been run, for the toolchain reason recorded in its
#' header.
#'
#' @inheritSection bandit2arm_delta The Laplace caveat
#' @inheritParams bandit2arm_delta
#'
#' @return A `frmtmb_family` object, for `frm(family = )`.
#'
#' @references
#' Glascher, J., Hampton, A. N. and O'Doherty, J. P. (2009). Determining
#' a role for ventromedial prefrontal cortex in encoding action-based
#' value signals during reward-related decision making. *Cerebral
#' Cortex* 19, 483-495.
#'
#' @seealso [bandit2arm_delta()], [bandit2arm_dual()]
#'
#' @examples
#' d <- frm_task_design("reversal", n_subject = 6, n_trial = 40, seed = 4)
#' # a reversal task pays -1 as often as +1
#' d$pay1 <- 2 * d$pay1 - 1
#' d$pay2 <- 2 * d$pay2 - 1
#' d$choice <- frm_task_simulate(
#'   prl_fictitious(subject = id, trial = trial), d,
#'   pars = list(alpha = 0.3, bias = 0, tau = 3), seed = 4)[[1]]$choice
#' fit <- frmtmb::frm(
#'   frmtmb::bf(choice | reward(pay1, pay2) ~ 1, bias ~ 1, tau ~ 1),
#'   family = prl_fictitious(subject = id, trial = trial), data = d)
#' frmtmb::fixef(fit)
#' @export
prl_fictitious <- function(subject, trial = NULL) {
  spec <- ln_spec(
    n_option = 2L,
    init = function(ns, d1) list(ev1 = rep(0, ns), ev2 = rep(0, ns)),
    choice = function(state, d, j) {
      list(d[["tau"]] * state[["ev1"]] + d[["bias"]],
           d[["tau"]] * state[["ev2"]])
    },
    update = function(state, d, ch) {
      ev1 <- state[["ev1"]]
      ev2 <- state[["ev2"]]
      c1 <- ch[[1L]][[1L]]
      c2 <- ch[[1L]][[2L]]
      # the outcome the subject saw, and its counterfactual
      oc <- c1 * d[["reward1"]] + c2 * d[["reward2"]]
      a <- d[["alpha"]]
      # option 1 is updated with `oc` when it was chosen and with -oc
      # when it was not, and option 2 the other way round
      t1 <- c1 * oc - c2 * oc
      t2 <- c2 * oc - c1 * oc
      list(state = list(ev1 = ev1 + a * (t1 - ev1),
                        ev2 = ev2 + a * (t2 - ev2)),
           pe = c1 * (t1 - ev1) + c2 * (t2 - ev2))
    })
  ln_family("prl_fictitious", substitute(subject), substitute(trial),
            dpars = c("alpha", "bias", "tau"),
            links = list(alpha = "logit", bias = "identity", tau = "log"),
            primary = "alpha",
            inits = list(alpha = function(y, aterms) 0.3,
                         bias = function(y, aterms) 0,
                         tau = function(y, aterms) 1),
            aterms = c("reward1", "reward2"), spec = spec,
            data_map = c(reward1 = "pay1", reward2 = "pay2"))
}
