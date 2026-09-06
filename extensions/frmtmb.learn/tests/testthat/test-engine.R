## The engine, checked against things that are not the engine.
##
## A recursion vectorized across subjects, masked for unequal trial
## counts and selecting the chosen option by multiplying by an indicator
## is exactly the kind of code that can be wrong in a way that still
## converges. Every test here compares it with the model written out the
## way the model is STATED: one subject, one trial, plain arithmetic.

ln_ref_delta <- function(d, alpha, tau) {
  ll <- 0
  for (s in levels(d$id)) {
    rows <- which(d$id == s)
    rows <- rows[order(d$trial[rows])]
    q <- c(0, 0)
    for (i in rows) {
      p <- exp(tau * q) / sum(exp(tau * q))
      ll <- ll + log(p[d$choice[i]])
      k <- d$choice[i]
      pay <- if (k == 1) d$pay1[i] else d$pay2[i]
      q[k] <- q[k] + alpha * (pay - q[k])
    }
  }
  ll
}

ln_toy <- function(ns = 8L, nt = 30L, seed = 41L, alpha = 0.35, tau = 3) {
  d <- frm_task_design("bandit2arm", n_subject = ns, n_trial = nt,
                       seed = seed)
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = id, trial = trial), d,
    pars = list(alpha = alpha, tau = tau), seed = seed)[[1L]]$choice
  d
}

test_that("the recursion equals the model written out longhand", {
  d <- ln_toy()
  fit <- frmtmb::frm(
    frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
    family = bandit2arm_delta(subject = id, trial = trial), data = d)
  fx <- unlist(frmtmb::fixef(fit))
  ref <- ln_ref_delta(d, stats::plogis(fx[["alpha.(Intercept)"]]),
                      exp(fx[["tau.(Intercept)"]]))
  expect_equal(ref, as.numeric(stats::logLik(fit)), tolerance = 1e-10)
})

test_that("the trace's p is the probability of the choice that was made", {
  d <- ln_toy()
  fit <- frmtmb::frm(
    frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
    family = bandit2arm_delta(subject = id, trial = trial), data = d)
  p <- frm_value_trace(fit)$p
  expect_length(p, nrow(d))
  expect_true(all(p > 0 & p < 1))
  # The identity that ties the per-row quantity to the total, and the
  # reason this package can refuse fitted() without losing anything:
  # the per-trial factors are all here and they multiply to the
  # likelihood.
  expect_equal(sum(log(p)), as.numeric(stats::logLik(fit)),
               tolerance = 1e-10)
  # fitted() refuses rather than returning arithmetic on a nominal
  # option code; see the fitted row of frm_compat()
  expect_error(stats::fitted(fit), "bandit2arm_delta")
})

test_that("padding is inert: an unbalanced design scores its own rows", {
  d <- ln_toy(ns = 6L, nt = 40L, seed = 42L)
  # drop a different number of late trials from each subject
  keep <- rep(TRUE, nrow(d))
  for (k in seq_along(levels(d$id))) {
    rows <- which(d$id == levels(d$id)[k])
    keep[rows[d$trial[rows] > 40 - 5 * (k - 1)]] <- FALSE
  }
  du <- d[keep, ]
  fit <- frmtmb::frm(
    frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
    family = bandit2arm_delta(subject = id, trial = trial), data = du)
  fx <- unlist(frmtmb::fixef(fit))
  ref <- ln_ref_delta(droplevels(du), stats::plogis(fx[["alpha.(Intercept)"]]),
                      exp(fx[["tau.(Intercept)"]]))
  expect_equal(ref, as.numeric(stats::logLik(fit)), tolerance = 1e-10)
  expect_equal(stats::nobs(fit), nrow(du))
})

test_that("one trial per subject is not read as one subject with n trials", {
  # The transposition trap: vapply() drops to a plain vector when the
  # longest subject has one trial, and t() of a vector is 1-by-n_subj,
  # so six subjects come back as one learner with six trials and the
  # value store is carried across subject boundaries.
  d <- frm_task_design("bandit2arm", n_subject = 6L, n_trial = 1L,
                       seed = 43L)
  d$choice <- c(1L, 2L, 1L, 1L, 2L, 1L)
  fit <- frmtmb::frm(
    frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
    family = bandit2arm_delta(subject = id, trial = trial), data = d)
  # every subject's FIRST trial starts from Q1 = Q2 = 0, so the choice
  # probability is 0.5 whatever the parameters are. A test on the block
  # dimensions alone would pass on a block that was right by accident;
  # this one fails on the arithmetic.
  expect_equal(as.numeric(stats::logLik(fit)), 6 * log(0.5),
               tolerance = 1e-8)
  blk <- frmtmb::frame_block_of(fit$frame, "choice")
  expect_equal(dim(blk$idx), c(6L, 1L))
})

test_that("a subject's value store does not leak into the next subject", {
  # Two subjects with identical data must score exactly twice one
  # subject's log-likelihood. They do not if the recursion carries Q
  # across the boundary.
  d1 <- ln_toy(ns = 1L, nt = 25L, seed = 44L)
  d2 <- rbind(d1, transform(d1, id = factor(2)))
  d2$id <- factor(d2$id)
  f1 <- frmtmb::frm(
    frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
    family = bandit2arm_delta(subject = id, trial = trial), data = d1)
  f2 <- frmtmb::frm(
    frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
    family = bandit2arm_delta(subject = id, trial = trial), data = d2)
  expect_equal(2 * as.numeric(stats::logLik(f1)),
               as.numeric(stats::logLik(f2)), tolerance = 1e-6)
})

test_that("the trace records the value the choice was made on", {
  d <- ln_toy(ns = 4L, nt = 20L, seed = 45L)
  fit <- frmtmb::frm(
    frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
    family = bandit2arm_delta(subject = id, trial = trial), data = d)
  tr <- frm_value_trace(fit)
  expect_true(all(c("subject", "trial", "q1", "q2", "pe", "p") %in%
                    names(tr)))
  first <- tr[tr$trial == 1, ]
  # both values start at zero, so the first trial's choice probability
  # is one over the number of options whatever the parameters are
  expect_true(all(first$q1 == 0 & first$q2 == 0))
  expect_true(all(abs(first$p - 0.5) < 1e-12))
  # and the recorded probability reproduces the recorded values
  tau <- exp(unlist(frmtmb::fixef(fit))[["tau.(Intercept)"]])
  pk <- ifelse(d$choice == 1,
               stats::plogis(tau * (tr$q1 - tr$q2)),
               stats::plogis(tau * (tr$q2 - tr$q1)))
  expect_equal(pk, tr$p, tolerance = 1e-10)
})

test_that("the two simulation routes agree", {
  # frm_task_simulate() takes parameters directly and builds its own
  # block; frm_simulate() goes through the fitted grammar. They share
  # the recursion, so at one seed and one parameter set they must give
  # the same draw.
  d <- ln_toy(ns = 5L, nt = 20L, seed = 46L)
  fam <- bandit2arm_delta(subject = id, trial = trial)
  set.seed(99)
  a <- frm_task_simulate(fam, d, pars = list(alpha = 0.4, tau = 2),
                         seed = 7)[[1L]]$choice
  set.seed(99)
  b <- frmtmb::frm_simulate(
    frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1), data = d,
    family = fam,
    newparams = list(alpha_Intercept = stats::qlogis(0.4),
                     tau_Intercept = log(2)), nsim = 1, seed = 7)[[1L]]
  expect_identical(as.integer(a), as.integer(b))
})

test_that("the block refuses data it cannot order", {
  d <- ln_toy(ns = 3L, nt = 10L, seed = 47L)
  bad <- d
  bad$trial[2] <- bad$trial[1]
  expect_error(
    frmtmb::frm(frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
                family = bandit2arm_delta(subject = id, trial = trial),
                data = bad),
    "unique within a subject")
  # An NA in the trial column never reaches the block's own guard: the
  # trial variable is in the model frame (the family puts it there
  # through frame_vars), so na.action drops the row first and says so.
  # That is the right outcome and the guard behind it is unreachable
  # through this path, which is recorded rather than removed: a trial
  # expression that COMPUTES an NA from non-missing columns would still
  # meet it.
  bad2 <- d
  bad2$trial[3] <- NA
  expect_message(
    f2 <- frmtmb::frm(frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
                      family = bandit2arm_delta(subject = id,
                                                trial = trial),
                      data = bad2),
    "removed because of missing values")
  expect_equal(stats::nobs(f2), nrow(d) - 1L)
})

test_that("the response coding is checked and named", {
  d <- ln_toy(ns = 3L, nt = 10L, seed = 48L)
  d$choice <- d$choice - 1L
  expect_error(
    frmtmb::frm(frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
                family = bandit2arm_delta(subject = id, trial = trial),
                data = d),
    "coded 1 to 2")
})
