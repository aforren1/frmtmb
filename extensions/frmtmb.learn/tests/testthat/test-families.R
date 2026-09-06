## Each family: that it fits, that its own mechanism does what its help
## says, and that the number it returns is the number the model defines.
## The identity against Stan is in test-stan-identity.R and is gated;
## these run everywhere.

ln_fit <- function(fam, d, ...) {
  frmtmb::frm(do.call(frmtmb::bf, list(...)), family = fam, data = d)
}

test_that("bandit2arm_delta learns the better arm", {
  d <- frm_task_design("bandit2arm", n_subject = 10L, n_trial = 60L,
                       seed = 51L)
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = id, trial = trial), d,
    pars = list(alpha = 0.4, tau = 4), seed = 51)[[1L]]$choice
  # arm 1 pays with probability 0.7 and arm 2 with 0.3, so a learner
  # takes arm 1 most of the time. This is a property of the SIMULATOR,
  # and it is asserted because a recursion that never learned would
  # still fit.
  expect_gt(mean(d$choice == 1), 0.6)
  # and late trials more than early ones
  expect_gt(mean(d$choice[d$trial > 40] == 1),
            mean(d$choice[d$trial <= 10] == 1))
  fit <- ln_fit(bandit2arm_delta(subject = id, trial = trial), d,
                choice | reward(pay1, pay2) ~ 1, tau ~ 1)
  fx <- unlist(frmtmb::fixef(fit))
  expect_gt(stats::plogis(fx[["alpha.(Intercept)"]]), 0.15)
  expect_lt(stats::plogis(fx[["alpha.(Intercept)"]]), 0.75)
  expect_gt(exp(fx[["tau.(Intercept)"]]), 1.5)
})

test_that("a covariate on the learning rate is the reversal model", {
  # The worked example of what the grammar buys: one fit, one interval,
  # for the CHANGE in learning rate after a reversal.
  d <- frm_task_design("reversal", n_subject = 14L, n_trial = 60L,
                       seed = 52L)
  ab <- as.numeric(d$after_reversal == "after")
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = id, trial = trial), d,
    pars = list(alpha = 0.3, tau = 3), seed = 52)[[1L]]$choice
  fit <- ln_fit(bandit2arm_delta(subject = id, trial = trial), d,
                choice | reward(pay1, pay2) ~ after_reversal, tau ~ 1)
  fx <- unlist(frmtmb::fixef(fit))
  expect_true("alpha.after_reversalafter" %in% names(fx))
  # the effect has a standard error, which is the point
  ci <- stats::confint(fit)
  expect_true(any(grepl("after_reversalafter", rownames(ci))))
  expect_true(all(is.finite(ci[grepl("after_reversalafter", rownames(ci)),
                               , drop = FALSE])))
  # a smooth reaches the same parameter
  fs <- ln_fit(bandit2arm_delta(subject = id, trial = trial), d,
               choice | reward(pay1, pay2) ~ s(trial, k = 5),
               tau ~ 1)
  expect_true(is.finite(as.numeric(stats::logLik(fs))))
})

test_that("bandit2arm_dual: the two splits coincide on binary payoffs", {
  # Not a coincidence and worth pinning. With payoffs in {0, 1} the
  # value store stays in [0, 1), so the prediction error is positive
  # exactly when the payoff is 1, and the two splits select the same
  # rate on every trial. They are different models only when the payoff
  # is GRADED.
  d <- frm_task_design("bandit2arm", n_subject = 8L, n_trial = 50L,
                       seed = 53L)
  d$choice <- frm_task_simulate(
    bandit2arm_dual(subject = id, trial = trial), d,
    pars = list(Arew = 0.5, Apun = 0.15, tau = 3), seed = 53)[[1L]]$choice
  fpe <- ln_fit(bandit2arm_dual(subject = id, trial = trial), d,
                choice | reward(pay1, pay2) ~ 1, Apun ~ 1, tau ~ 1)
  fout <- ln_fit(bandit2arm_dual(subject = id, trial = trial,
                                 split = "outcome"), d,
                 choice | reward(pay1, pay2) ~ 1, Apun ~ 1, tau ~ 1)
  expect_equal(as.numeric(stats::logLik(fpe)),
               as.numeric(stats::logLik(fout)), tolerance = 1e-8)

  # graded payoffs separate them
  set.seed(53)
  dg <- d
  dg$pay1 <- stats::runif(nrow(dg), -1, 1)
  dg$pay2 <- stats::runif(nrow(dg), -1, 1)
  gpe <- ln_fit(bandit2arm_dual(subject = id, trial = trial), dg,
                choice | reward(pay1, pay2) ~ 1, Apun ~ 1, tau ~ 1)
  gout <- ln_fit(bandit2arm_dual(subject = id, trial = trial,
                                 split = "outcome"), dg,
                 choice | reward(pay1, pay2) ~ 1, Apun ~ 1, tau ~ 1)
  expect_false(isTRUE(all.equal(as.numeric(stats::logLik(gpe)),
                                as.numeric(stats::logLik(gout)))))
})

test_that("bandit2arm_dual recovers an asymmetry that is there", {
  d <- frm_task_design("bandit2arm", n_subject = 25L, n_trial = 120L,
                       seed = 54L)
  d$pay1 <- 2 * d$pay1 - 1
  d$pay2 <- 2 * d$pay2 - 1
  d$choice <- frm_task_simulate(
    bandit2arm_dual(subject = id, trial = trial), d,
    pars = list(Arew = 0.6, Apun = 0.1, tau = 3), seed = 54)[[1L]]$choice
  fit <- ln_fit(bandit2arm_dual(subject = id, trial = trial), d,
                choice | reward(pay1, pay2) ~ 1, Apun ~ 1, tau ~ 1)
  fx <- unlist(frmtmb::fixef(fit))
  expect_gt(stats::plogis(fx[["Arew.(Intercept)"]]),
            stats::plogis(fx[["Apun.(Intercept)"]]))
})

test_that("prl_fictitious moves the option that was not chosen", {
  d <- frm_task_design("reversal", n_subject = 8L, n_trial = 40L,
                       seed = 55L)
  d$pay1 <- 2 * d$pay1 - 1
  d$pay2 <- 2 * d$pay2 - 1
  d$choice <- frm_task_simulate(
    prl_fictitious(subject = id, trial = trial), d,
    pars = list(alpha = 0.35, bias = 0, tau = 3), seed = 55)[[1L]]$choice
  fit <- ln_fit(prl_fictitious(subject = id, trial = trial), d,
                choice | reward(pay1, pay2) ~ 1, bias ~ 1, tau ~ 1)
  tr <- frm_value_trace(fit)
  # after trial 1 both values have moved, which is what counterfactual
  # updating means: a delta rule would leave the unchosen one at zero
  t2 <- tr[tr$trial == 2, ]
  expect_true(all(t2$ev1 != 0))
  expect_true(all(t2$ev2 != 0))
  # and they move in opposite directions
  expect_true(all(sign(t2$ev1) == -sign(t2$ev2)))
})

test_that("the Kalman gain falls as a subject learns", {
  d <- frm_task_design("bandit4arm_restless", n_subject = 6L, n_trial = 50L,
                       seed = 56L)
  fam <- bandit4arm2_kalman_filter(subject = id, trial = trial)
  d$choice <- frm_task_simulate(
    fam, d, pars = list(tau = 0.15, lambda = 0.98, center = 50, mu0 = 50,
                        sigma0 = 12, sigmaD = 2), seed = 56)[[1L]]$choice
  fit <- ln_fit(fam, d, choice | payoff(pay1, pay2, pay3, pay4) ~ 1,
                lambda ~ 1, center ~ 1, mu0 ~ 1, sigma0 ~ 1, sigmaD ~ 1)
  tr <- frm_value_trace(fit)
  expect_true(all(c("mu1", "s1", "s4") %in% names(tr)))
  # the posterior variance of an arm starts at sigma0^2 and is pulled
  # down every time that arm is pulled. This is the whole point of the
  # family: the learning rate is not a parameter, it is the variance.
  v1 <- tr$s1[tr$trial == 1]
  vlast <- tr$s1[tr$trial == max(tr$trial)]
  expect_true(all(v1 > 0))
  expect_lt(mean(vlast), mean(v1))
  # and it never goes negative, which the update guarantees rather than
  # a floor: the gain is in (0, 1) by construction
  expect_true(all(tr$s1 > 0, na.rm = TRUE))
})

test_that("igt_pvl_delta learns away from the frequent-loss deck", {
  d <- frm_task_design("igt", n_subject = 10L, n_trial = 80L, seed = 57L)
  fam <- igt_pvl_delta(subject = id, trial = trial)
  d$choice <- frm_task_simulate(
    fam, d, pars = list(alpha = 0.3, shape = 0.4, lambda = 2, tau = 1.5),
    seed = 57)[[1L]]$choice
  # Deck 1 is the one this MODEL learns away from, and it is not simply
  # the deck with the worst objective expected value. Under the
  # prospect utility at shape 0.4 and loss aversion 2, deck 2's rare
  # large loss is compressed to a fraction of its size, so decks 2, 3
  # and 4 come out close together in subjective value while deck 1,
  # whose loss is frequent, sits far below all three. That is the
  # well-known deck-B effect: it is a property of the model rather than
  # a fault in the design, so the assertion follows the model.
  # Measured on this fixture: deck 1 takes 0.270 of the first ten
  # trials, which is chance, and 0.100 of the last twenty. Comparing
  # against ALL earlier trials would compare an asymptote with itself,
  # because most of the learning is done by trial 20.
  expect_gt(mean(d$choice[d$trial <= 10] == 1L),
            mean(d$choice[d$trial > 60] == 1L))
  # and deck 1 ends up the least chosen of the four
  expect_equal(unname(which.min(table(d$choice))), 1L)
  fit <- ln_fit(fam, d, choice | payoff(pay1, pay2, pay3, pay4) ~ 1,
                shape ~ 1, lambda ~ 1, tau ~ 1)
  tr <- frm_value_trace(fit)
  expect_true("utility" %in% names(tr))
  # the outcome transform: a loss is scaled by loss aversion and a gain
  # is not, so a negative outcome maps to a utility more negative than
  # its magnitude^shape
  x <- rowSums(as.matrix(d[, paste0("pay", 1:4)]) *
                 stats::model.matrix(~ 0 + factor(d$choice)))
  neg <- x < 0
  expect_true(all(tr$utility[neg] < 0))
  expect_true(all(tr$utility[!neg & x > 0] > 0))
  expect_equal(sign(tr$utility), sign(x))
})

test_that("ts_par7 is model-based when w says so", {
  d <- frm_task_design("twostep", n_subject = 10L, n_trial = 80L, seed = 58L)
  fam <- ts_par7(subject = id, trial = trial)
  base <- list(alpha1 = 0.5, tau1 = 5, alpha2 = 0.5, tau2 = 5,
               lambda = 0.5, pers = 0)
  dmb <- frm_task_simulate(fam, d, pars = c(list(w = 0.99), base),
                           seed = 58)[[1L]]
  dmf <- frm_task_simulate(fam, d, pars = c(list(w = 0.01), base),
                           seed = 58)[[1L]]
  # The signature of the two-step task: a model-free learner stays after
  # ANY reward, a model-based one stays after a common-transition reward
  # and switches after a rare-transition reward. Measure the stay rate
  # split by transition type on the previous trial.
  stay_rate <- function(dd) {
    out <- c(common = NA_real_, rare = NA_real_)
    common <- (dd$choice == 1 & dd$state2 == 1) |
      (dd$choice == 2 & dd$state2 == 2)
    rew <- vapply(seq_len(nrow(dd)), function(i) {
      dd[[paste0("pay", 2 * (dd$state2[i] - 1) + dd$choice2[i])]][i]
    }, numeric(1))
    nxt <- c(dd$choice[-1], NA)
    same_subj <- c(dd$id[-1] == dd$id[-nrow(dd)], FALSE)
    stay <- nxt == dd$choice & same_subj
    ok <- rew == 1 & same_subj & !is.na(stay)
    out[["common"]] <- mean(stay[ok & common])
    out[["rare"]] <- mean(stay[ok & !common])
    out
  }
  smb <- stay_rate(dmb)
  smf <- stay_rate(dmf)
  # a model-based learner stays more after a common reward than a rare
  # one; a model-free learner does not care which
  expect_gt(smb[["common"]] - smb[["rare"]],
            smf[["common"]] - smf[["rare"]])
  fit <- ln_fit(fam, dmb,
                choice | stage2(state2, choice2) +
                  payoff(pay1, pay2, pay3, pay4) ~ 1,
                alpha1 ~ 1, tau1 ~ 1, alpha2 ~ 1, tau2 ~ 1, lambda ~ 1,
                pers ~ 1)
  # and the fit recovers a high w from model-based data
  expect_gt(stats::plogis(unlist(frmtmb::fixef(fit))[["w.(Intercept)"]]),
            0.5)
})

test_that("family constructors refuse what they cannot do", {
  expect_error(bandit4arm2_kalman_filter(subject = id, sigma_o = -1),
               "one positive number")
  expect_error(ts_par7(subject = id, p_common = 0.5),
               "strictly between 0.5 and 1")
  expect_error(frm_task_design("bandit2arm", n_subject = 0),
               "at least 1")
  fam <- bandit2arm_delta(subject = id, trial = trial)
  d <- frm_task_design("bandit2arm", n_subject = 3L, n_trial = 5L, seed = 59)
  expect_error(frm_task_simulate(fam, d, pars = list(alpha = 0.3)),
               "no value for tau")
  expect_error(frm_task_simulate(fam, d, pars = list(alpha = rep(0.3, 7),
                                                     tau = 1)),
               "7 values")
  expect_error(frm_task_simulate(stats::gaussian(), d, pars = list()),
               "carries no learning recursion")
})


test_that("a parameter may vary within a subject, one value per row", {
  # The case the package exists for: a reversal changes the learning
  # rate PART WAY THROUGH a session, so the simulator has to accept a
  # value per row and not only per subject.
  d <- frm_task_design("reversal", n_subject = 8L, n_trial = 60L,
                       seed = 60L)
  a <- ifelse(d$after_reversal == "after", 0.65, 0.15)
  fam <- bandit2arm_delta(subject = id, trial = trial)
  d$choice <- frm_task_simulate(fam, d, pars = list(alpha = a, tau = 3),
                                seed = 60)[[1L]]$choice
  fit <- ln_fit(fam, d, choice | reward(pay1, pay2) ~ after_reversal,
                tau ~ 1)
  fx <- unlist(frmtmb::fixef(fit))
  # the effect is in the right direction and reaches significance
  expect_gt(fx[["alpha.after_reversalafter"]], 0)
  ci <- stats::confint(fit)
  row <- ci[grepl("after_reversalafter", rownames(ci)), , drop = FALSE]
  expect_gt(row[1L, 1L], 0)
  # a per-row vector of the wrong length is still refused
  expect_error(frm_task_simulate(fam, d, pars = list(alpha = a[-1], tau = 1)),
               "one per row")
})