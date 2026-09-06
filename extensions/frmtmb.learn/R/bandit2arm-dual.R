#' Two-armed delta learning with separate rates for gains and losses
#'
#' [bandit2arm_delta()] with one learning rate replaced by two, so that
#' good news and bad news are absorbed at different speeds. This
#' asymmetry is the parameter most of the clinical literature is about,
#' and it is the model hBayesDM calls `prl_rp` when the split is on the
#' sign of the OUTCOME.
#'
#' ```
#' pe  <- reward - Q[chosen]
#' Q[chosen] <- Q[chosen] + (pe > 0 ? Arew : Apun) * pe
#' ```
#'
#' @section Which sign, and why it matters here more than elsewhere:
#' Two models are in circulation and they are not the same model. One
#' splits on the sign of the PREDICTION ERROR: an outcome better than
#' expected is a gain even when it is a small reward. The other splits
#' on the sign of the OUTCOME, which is what `prl_rp` does. `split`
#' chooses, and the default is `"pe"`.
#'
#' They also cost different things to tape. The outcome's sign is DATA,
#' so `split = "outcome"` selects the rate with a column of zeros and
#' ones computed once and adds nothing to the tape. A prediction error
#' is an AD quantity, and RTMB refuses a comparison on one outright, so
#' `split = "pe"` selects with `sign()`. That is exact rather than
#' smoothed, and the kink it introduces at `pe == 0` is the model's own:
#' the update is zero there under both rates, so the two branches agree
#' at the crossing and the derivative from each side is the right one.
#'
#' With `split = "outcome"`, an outcome of exactly zero counts as
#' punishment, which is what `prl_rp` does and what the usual `0`/`1`
#' payoff coding needs. Code punishment as `-1` if you want the
#' symmetric reading.
#'
#' One consequence of the kink is worth knowing before fitting `"pe"`
#' hierarchically on GRADED payoffs. The joint log density is then
#' non-smooth in the random effects wherever a prediction error crosses
#' zero, and the Laplace inner solve stops short of the conditional
#' mode: measured on the identity fixture, the largest gradient on the
#' subject effects is 4.0e-02 under `"pe"` against 1.4e-15 under
#' `"outcome"` on the same data. The fit and its log likelihood are
#' unaffected to machine precision, but a convergence warning on a
#' `"pe"` fit with graded payoffs is expected rather than alarming. On
#' binary payoffs neither split has a live kink, and the two are the
#' same model anyway.
#'
#' @inheritSection bandit2arm_delta The Laplace caveat
#' @inheritParams bandit2arm_delta
#' @param split `"pe"` splits on the sign of the prediction error,
#'   `"outcome"` on the sign of the outcome itself.
#'
#' @return A `frmtmb_family` object, for `frm(family = )`.
#'
#' @references
#' Ahn, W.-Y., Haines, N. and Zhang, L. (2017). Revealing
#' neurocomputational mechanisms of reinforcement learning and
#' decision-making with the hBayesDM package. *Computational
#' Psychiatry* 1, 24-57.
#'
#' @seealso [bandit2arm_delta()], [frm_learn_families()]
#'
#' @examples
#' d <- frm_task_design("bandit2arm", n_subject = 6, n_trial = 40, seed = 2)
#' d$choice <- frm_task_simulate(
#'   bandit2arm_dual(subject = id, trial = trial), d,
#'   pars = list(Arew = 0.5, Apun = 0.15, tau = 3),
#'   seed = 2)[[1]]$choice
#' fit <- frmtmb::frm(
#'   frmtmb::bf(choice | reward(pay1, pay2) ~ 1, Apun ~ 1, tau ~ 1),
#'   family = bandit2arm_dual(subject = id, trial = trial), data = d)
#' frmtmb::fixef(fit)
#' @export
bandit2arm_dual <- function(subject, trial = NULL,
                            split = c("pe", "outcome")) {
  split <- match.arg(split)
  upd <- function(state, d, ch) {
    q1 <- state[["q1"]]
    q2 <- state[["q2"]]
    c1 <- ch[[1L]][[1L]]
    c2 <- ch[[1L]][[2L]]
    pe1 <- d[["reward1"]] - q1
    pe2 <- d[["reward2"]] - q2
    pec <- c1 * pe1 + c2 * pe2
    w <- if (split == "outcome") {
      as.numeric(c1 * d[["reward1"]] + c2 * d[["reward2"]] > 0)
    } else {
      0.5 * (1 + sign(pec))
    }
    a <- w * d[["Arew"]] + (1 - w) * d[["Apun"]]
    list(state = list(q1 = q1 + c1 * a * pe1, q2 = q2 + c2 * a * pe2),
         pe = pec, rate = a)
  }
  spec <- ln_spec(
    n_option = 2L,
    init = function(ns, d1) list(q1 = rep(0, ns), q2 = rep(0, ns)),
    choice = function(state, d, j) {
      list(d[["tau"]] * state[["q1"]], d[["tau"]] * state[["q2"]])
    },
    update = upd)
  ln_family("bandit2arm_dual", substitute(subject), substitute(trial),
            dpars = c("Arew", "Apun", "tau"),
            links = list(Arew = "logit", Apun = "logit",
                         tau = "log"),
            primary = "Arew",
            inits = list(Arew = function(y, aterms) 0.3,
                         Apun = function(y, aterms) 0.3,
                         tau = function(y, aterms) 1),
            aterms = c("reward1", "reward2"), spec = spec,
            data_map = c(reward1 = "pay1", reward2 = "pay2"),
            constants = list(split = split))
}
