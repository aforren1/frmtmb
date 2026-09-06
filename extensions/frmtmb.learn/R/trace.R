#' Per-trial value estimates, prediction errors and choice probabilities
#'
#' The quantities these papers plot. A learning model is fitted for its
#' latent trajectory as much as for its parameters, and the trajectory
#' is not in the fitted object: it is recomputed by replaying the
#' recursion at the estimates.
#'
#' The columns after `trial` are the family's own. Every family returns
#' `p`, the fitted probability of the choice that was actually made, and
#' `pe`, the prediction error the outcome produced. The value stores
#' vary: `q1` and `q2` for a two-armed family, `mu1` to `mu4` and `s1`
#' to `s4` for the Kalman filter, `qmf1` and `qmb1` for the two-step
#' learner's model-free and model-based stage-one values. Each value is
#' the one the choice on that trial was made ON, before the outcome
#' moved it, which is the ordering a plot of learning needs.
#'
#' @section What `p` is, and what it is not:
#' `p` is the per-trial factor of the likelihood, so `sum(log(p))` is
#' the CONDITIONAL data log-likelihood given the fitted subject effects.
#' In a fit with no random effects that is `logLik(fit)` exactly, to
#' machine precision, and the test suite asserts it.
#'
#' In a HIERARCHICAL fit the two differ, and by more than rounding:
#' `logLik()` reports the Laplace-approximated MARGINAL likelihood,
#' which adds the random-effect density and the curvature term that a
#' conditional quantity does not carry. Measured on a 30-subject fit
#' with `sd(id) = 1.07`, `sum(log(p))` is -1507.8 against a `logLik()`
#' of -1542.5. Neither is wrong; they are different quantities, and only
#' the conditional one exists per trial.
#'
#' This function REPLACES [stats::fitted()] for these families rather
#' than supplementing it. They decline to declare a mean on the response
#' scale, because the response is a nominal option code and `y - mu`
#' would be arithmetic on a category, so `fitted()` and
#' `predict(type = "response")` refuse. Everything a mean would have
#' carried is in this table, with more beside it.
#'
#' @param fit A fit whose family came from this package.
#'
#' @return A data frame with one row per row of the model frame, in the
#'   data's own order.
#'
#' @seealso [bandit2arm_delta()], [frm_task_simulate()]
#'
#' @examples
#' d <- frm_task_design("bandit2arm", n_subject = 6, n_trial = 30,
#'                      seed = 3)
#' d$choice <- frm_task_simulate(
#'   bandit2arm_delta(subject = id, trial = trial), d,
#'   pars = list(alpha = 0.4, tau = 3), seed = 3)[[1]]$choice
#' fit <- frmtmb::frm(frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
#'                    family = bandit2arm_delta(subject = id,
#'                                              trial = trial), data = d)
#' tr <- frm_value_trace(fit)
#' head(tr)
#' # the identity that ties the trace to the likelihood
#' c(from_trace = sum(log(tr$p)), logLik = as.numeric(stats::logLik(fit)))
#' @export
frm_value_trace <- function(fit) {
  rspec <- single_response(fit, "a frmtmb.learn fit")
  fam <- rspec[["family"]]
  lrn <- fam[["learn"]]
  if (is.null(lrn)) {
    stop("frm_value_trace() reads a fit whose family came from ",
         "frmtmb.learn, and this one is a '", fam[["family"]],
         "' fit. There is no value recursion to replay", call. = FALSE)
  }
  blk <- frame_block_of(fit[["frame"]], rspec[["resp_name"]])
  tr <- ln_trace_at(fit, blk, lrn[["spec"]], fam[["family"]])
  out <- data.frame(subject = blk[["subject"]], trial = blk[["trial"]])
  # p last, because it is the summary of the row rather than part of the
  # trajectory, and the value stores read left to right in trial order
  for (nm in c(setdiff(names(tr), "p"), "p")) out[[nm]] <- tr[[nm]]
  out
}

#' The families this package supplies, and what each one is
#'
#' A reference table rather than a computation. It carries the name each
#' family has here, the name hBayesDM gives the same model where it has
#' one, the task, the options per trial, the addition terms and the
#' parameters with their links.
#'
#' The parameter names are this package's, not hBayesDM's, and the
#' `hbayesdm_pars` column is the map, in the same order as `pars`. A
#' cell that is a name is a RENAME; a cell that is an expression is a
#' TRANSFORM, and three of the six rows carry one. Two of those three
#' involve an estimated parameter (`-tau * alpha` for
#' [prl_fictitious()]'s `bias`, `pi / tau1` for [ts_par7()]'s `pers`),
#' so no constant relates them and a value carried across unchanged
#' fits a different model without complaint. `?prl_fictitious` and
#' `?ts_par7` derive both. One
#' convention across six families is worth more than six inherited ones:
#' `alpha` is always a learning rate on `(0, 1)` and `tau` is always the
#' softmax sensitivity on `(0, Inf)`, so what a formula does to `alpha`
#' means the same thing whichever family it is written against.
#' hBayesDM spells the two-armed bandit's learning rate `A`, the
#' reversal models' sensitivity `beta`, and the Iowa gambling task's
#' UTILITY EXPONENT `alpha`.
#'
#' Two cautions the map does not remove. Core refuses a distributional
#' parameter name containing a dot or an underscore, which is why the
#' dual-rate family spells its rates `Arew` and `Apun` rather than
#' anything more readable. And `lambda` is used by the literature for
#' three unrelated quantities: the arm decay in the Kalman filter, the
#' eligibility trace in the two-step model, and loss aversion in the
#' Iowa gambling task. Each family's help says which.
#'
#' @return A data frame, one row per family.
#'
#' @examples
#' frm_learn_families()[, c("family", "hbayesdm", "pars")]
#' @export
frm_learn_families <- function() {
  d <- data.frame(
    family = c("bandit2arm_delta", "bandit2arm_dual", "prl_fictitious",
               "bandit4arm2_kalman_filter", "ts_par7", "igt_pvl_delta"),
    hbayesdm = c("bandit2arm_delta", "prl_rp (split = 'outcome')",
                 "prl_fictitious", "bandit4arm2_kalman_filter", "ts_par7",
                 "igt_pvl_delta"),
    task = c("two-armed bandit", "two-armed bandit, reward and punishment",
             "probabilistic reversal, counterfactual updating",
             "restless four-armed bandit",
             "two-stage Markov decision task", "Iowa gambling task"),
    options = c("2", "2", "2", "4", "2 then 2", "4"),
    aterms = c("reward(pay1, pay2)", "reward(pay1, pay2)",
               "reward(pay1, pay2)", "payoff(pay1, ..., pay4)",
               "stage2(state, choice) + payoff(pay1, ..., pay4)",
               "payoff(pay1, ..., pay4)"),
    pars = c("alpha, tau", "Arew, Apun, tau", "alpha, bias, tau",
             "tau, lambda, center, mu0, sigma0, sigmaD",
             "w, alpha1, tau1, alpha2, tau2, lambda, pers",
             "alpha, shape, lambda, tau"),
    # Cells show a TRANSFORM wherever the relation is not a rename.
    # Two of the six are related to hBayesDM by a factor of an ESTIMATED
    # parameter, so no constant maps them and a ported estimate is a
    # different model that converges quietly. Each family's help derives
    # the expression.
    hbayesdm_pars = c("A, tau", "Arew, Apun, beta",
                      "eta, -tau * alpha, beta",
                      "beta, lambda, theta, mu0, sigma0, sigmaD",
                      "w, a1, beta1, a2, beta2, lambda, pi / tau1",
                      "A, alpha, lambda, 3^cons - 1"),
    reference = c("Ahn, Haines and Zhang (2017)",
                  "Ahn, Haines and Zhang (2017)",
                  "Glascher, Hampton and O'Doherty (2009)",
                  "Daw, O'Doherty, Dayan, Seymour and Dolan (2006)",
                  "Daw, Gershman, Seymour, Dayan and Dolan (2011)",
                  "Ahn, Busemeyer, Wagenmakers and Stout (2008)"),
    stringsAsFactors = FALSE)
  d
}
