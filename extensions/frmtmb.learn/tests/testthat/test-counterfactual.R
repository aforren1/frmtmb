## What a fit needs and what a DRAW needs are not the same data, and the
## difference is a whole column. Every family here reads only the chosen
## option's payoff, so a record of the received outcome alone is fitted
## EXACTLY right by passing that column once per option. A simulated
## subject chooses for itself, and the option it picks has no recorded
## payoff unless the schedule of every option is there. These tests pin
## both halves: the fit that is exact, and the draw that is refused, on
## all eight families rather than on a chosen five.

ct_design <- function(seed = 1L, n_subject = 6L, n_trial = 60L) {
  d <- frm_task_design("bandit2arm", n_subject = n_subject,
                       n_trial = n_trial, seed = seed)
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = id, trial = trial), d,
    pars = list(alpha = 0.4, tau = 3), seed = seed)[[1L]]$choice
  # what the subject actually received, which is all a one-column
  # record holds, and a pair of columns that agree with the record on
  # the chosen entry and hold noise everywhere else
  d$rec <- ifelse(d$choice == 1L, d$pay1, d$pay2)
  set.seed(seed + 500L)
  d$junk <- stats::rnorm(nrow(d), 100, 50)
  d$p1 <- ifelse(d$choice == 1L, d$pay1, d$junk)
  d$p2 <- ifelse(d$choice == 2L, d$pay2, d$junk)
  d
}

# One plausible value for each of a family's parameters, taken from its
# own starting values rather than written down here, so that a family
# renaming a parameter cannot make this file assert nothing.
ct_pars <- function(fam) {
  dp <- fam[["learn"]][["dpars"]]
  ps <- lapply(dp, function(k) fam[["init_dpars"]][[k]](1, list()))
  stats::setNames(ps, dp)
}

# The objective at ONE parameter vector. Every fit below starts from the
# same constant inits, so stopping the optimizer before its first step
# compares the objective itself rather than wherever three optimizer
# paths happened to land. That is what makes the comparison bitwise.
ct_at_start <- function(form, fam, data) {
  suppressWarnings(
    as.numeric(logLik(frm(form, family = fam, data = data,
                          control = frmtmb_control(
                            optCtrl = list(iter.max = 0L,
                                           eval.max = 1L))))))
}

test_that("the unchosen column cannot reach the objective at all", {
  d <- ct_design()
  fam <- bandit2arm_delta(subject = id, trial = trial)
  two <- ct_at_start(bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
                     fam, d)
  one <- ct_at_start(bf(choice | reward(rec, rec) ~ 1, tau ~ 1), fam, d)
  noise <- ct_at_start(bf(choice | reward(p1, p2) ~ 1, tau ~ 1), fam, d)
  # bitwise, not to a tolerance: the unchosen entry is multiplied by a
  # zero indicator, and zero times a finite number is exactly zero
  expect_identical(one, two)
  expect_identical(noise, two)
})

test_that("the received outcome alone gives the same fit", {
  d <- ct_design()
  two <- frm(bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
             family = bandit2arm_delta(subject = id, trial = trial),
             data = d)
  one <- frm(bf(choice | reward(rec, rec) ~ 1, tau ~ 1),
             family = bandit2arm_delta(subject = id, trial = trial),
             data = d)
  # the CONVERGED values agree to testthat's relative default rather
  # than bitwise: the objective is identical but gradient accumulation
  # over a tape whose constants changed can move the last place
  expect_equal(as.numeric(logLik(one)), as.numeric(logLik(two)))
  expect_equal(unlist(fixef(one)), unlist(fixef(two)))
})

test_that("prl_fictitious() reads one column too, bitwise", {
  # This family was exempted from the guard on the ground that its
  # counterfactual update reads both columns. It does not: the update
  # forms c1 * reward1 + c2 * reward2 with the CHOSEN indicators and
  # flips the sign, so the unchosen column never enters. The exemption
  # was a sixth silent wrong answer, and this is what disproves it.
  d <- ct_design(seed = 2L)
  fam <- prl_fictitious(subject = id, trial = trial)
  two <- ct_at_start(bf(choice | reward(pay1, pay2) ~ 1, bias ~ 1,
                        tau ~ 1), fam, d)
  one <- ct_at_start(bf(choice | reward(rec, rec) ~ 1, bias ~ 1,
                        tau ~ 1), fam, d)
  noise <- ct_at_start(bf(choice | reward(p1, p2) ~ 1, bias ~ 1,
                          tau ~ 1), fam, d)
  expect_identical(one, two)
  expect_identical(noise, two)
})

test_that("simulate() refuses the duplicated column and names it", {
  d <- ct_design()
  one <- frm(bf(choice | reward(rec, rec) ~ 1, tau ~ 1),
             family = bandit2arm_delta(subject = id, trial = trial),
             data = d)
  expect_error(simulate(one, nsim = 1L, seed = 3L), "reward\\(\\)")
  expect_error(simulate(one, nsim = 1L, seed = 3L),
               "bandit2arm_delta", fixed = TRUE)
  # the sentence has to say which column is absent, not only that
  # something is, and it has to say that newdata is not a way round it
  expect_error(simulate(one, nsim = 1L, seed = 3L), "second")
  expect_error(simulate(one, nsim = 1L, seed = 3L), "newdata")
  # and newdata really is not a way round it: the formula names one
  # column twice, so a real schedule supplied there is read through it
  expect_error(simulate(one, newdata = d, nsim = 1L, seed = 3L),
               "reward\\(\\)")
  # the de novo route reaches the same slot and refuses there too
  expect_error(
    frm_simulate(bf(choice | reward(rec, rec) ~ 1, tau ~ 1), d,
                 family = bandit2arm_delta(subject = id, trial = trial),
                 newparams = list(alpha_Intercept = 0,
                                  tau_Intercept = 0), nsim = 1L,
                 seed = 3L),
    "reward\\(\\)")
})

test_that("a real two-column schedule still draws", {
  d <- ct_design()
  two <- frm(bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
             family = bandit2arm_delta(subject = id, trial = trial),
             data = d)
  s <- simulate(two, nsim = 1L, seed = 4L)
  expect_equal(length(s[[1L]]), nrow(d))
  expect_true(all(s[[1L]] %in% c(1, 2)))
})

test_that("every reward() family refuses a duplicated schedule", {
  # All four, including prl_fictitious() and rlddm(). rlddm() refuses
  # simulate() for its own reason and frm_task_simulate() is the route
  # its refusal message names, so this is the route that has to check.
  d <- ct_design()
  dup <- d
  dup$pay2 <- dup$pay1
  for (nm in c("bandit2arm_delta", "bandit2arm_dual",
               "prl_fictitious", "rlddm")) {
    fam <- do.call(nm, list(subject = quote(id), trial = quote(trial)))
    ps <- ct_pars(fam)
    expect_error(frm_task_simulate(fam, dup, pars = ps, seed = 5L),
                 "reward\\(\\)", info = nm)
    expect_no_error(frm_task_simulate(fam, d, pars = ps, seed = 5L))
  }
})

test_that("every payoff() family refuses a duplicated schedule", {
  # igt_pvl_delta(), igt_orl(), bandit4arm2_kalman_filter() and
  # ts_par7(). The last has no simulate() either, and drew 800 rows
  # from a degenerate schedule before this guard reached it.
  igt <- frm_task_design("igt", n_subject = 4L, n_trial = 40L, seed = 7L)
  rst <- frm_task_design("bandit4arm_restless", n_subject = 4L,
                         n_trial = 40L, seed = 7L)
  two <- frm_task_design("twostep", n_subject = 4L, n_trial = 40L,
                         seed = 7L)
  dupe <- function(z) {
    for (k in 2:4) z[[paste0("pay", k)]] <- z$pay1
    z
  }
  cases <- list(
    list(fam = igt_pvl_delta(subject = id, trial = trial), d = igt),
    list(fam = igt_orl(subject = id, trial = trial), d = igt),
    list(fam = bandit4arm2_kalman_filter(subject = id, trial = trial),
         d = rst),
    list(fam = ts_par7(subject = id, trial = trial), d = two))
  for (z in cases) {
    nm <- z$fam[["family"]]
    ps <- ct_pars(z$fam)
    expect_error(frm_task_simulate(z$fam, dupe(z$d), pars = ps,
                                   seed = 8L),
                 "payoff\\(\\)", info = nm)
    expect_no_error(frm_task_simulate(z$fam, z$d, pars = ps, seed = 8L))
  }
})

test_that("a partial duplicate is not the signature", {
  # three of four decks duplicated: the data still says what one deck
  # not played would have paid, so the draw is not refused
  d <- frm_task_design("igt", n_subject = 4L, n_trial = 40L, seed = 7L)
  part <- d
  part$pay2 <- part$pay1
  part$pay3 <- part$pay1
  expect_no_error(
    frm_task_simulate(igt_pvl_delta(subject = id, trial = trial), part,
                      pars = list(alpha = 0.3, shape = 0.5,
                                  lambda = 1, tau = 1), seed = 8L))
})

test_that("stage2() is not read as a payoff schedule", {
  # ts_par7() carries stage2(state, choice), two columns that are not
  # what each option would have paid. A design whose state and choice
  # agree on every trial must still draw.
  d <- frm_task_design("twostep", n_subject = 4L, n_trial = 40L,
                       seed = 9L)
  d$state2 <- d$choice2
  expect_no_error(
    frm_task_simulate(ts_par7(subject = id, trial = trial), d,
                      pars = list(w = 0.5, alpha1 = 0.3, tau1 = 2,
                                  alpha2 = 0.3, tau2 = 2, lambda = 0.5,
                                  pers = 0), seed = 9L))
})

## ---- the rule, not only its answers on the eight families ----------
## Everything above pins what the guard DOES. These pin the rule that
## produces it, which is the part a derived guard can lose silently: a
## schedule registered under a name the derivation does not know
## returns NULL and nothing says so. The only place that is visible is
## this package's own registration table.

test_that("every multi-column term this package registers is classified", {
  multi <- names(ln_aterms)[ln_aterms >= 2L]
  expect_gt(length(multi), 0)
  known <- c(ln_schedule_terms, names(ln_not_schedule_terms))
  # A term registered without being put in one of the two vectors is a
  # family that can draw from a degenerate schedule in silence. Item
  # 5.1's bandit_delta(n_option = K) is the one that will meet this,
  # because a K-armed schedule for K outside 2 and 4 must register a
  # name of its own.
  expect_equal(setdiff(multi, known), character(0))
  # and the two vectors must not disagree with each other
  expect_equal(intersect(ln_schedule_terms,
                         names(ln_not_schedule_terms)), character(0))
  # every exclusion carries its reason, because the reason is the only
  # thing a reviewer can disagree with
  expect_true(all(nzchar(ln_not_schedule_terms)))
  expect_true(all(nchar(ln_not_schedule_terms) > 40))
})

test_that("the derivation reads the schedule terms and only those", {
  # a schedule term of any arity, and a non-schedule term of the same
  # shape, put through the predicate directly
  expect_equal(ln_counterfactual_of(c("reward1", "reward2")),
               list(c("reward1", "reward2")))
  expect_equal(ln_counterfactual_of(paste0("payoff", 1:4)),
               list(paste0("payoff", 1:4)))
  expect_null(ln_counterfactual_of(c("stage21", "stage22")))
  # a term the derivation has never heard of is NOT guarded, which is
  # the honest limit of a derived rule and the reason for the test
  # above rather than a claim that it cannot happen
  expect_null(ln_counterfactual_of(c("outcome1", "outcome2")))
  expect_null(ln_counterfactual_of(c("vreal1", "vreal2")))
  # a single column of a schedule term is not a group
  expect_null(ln_counterfactual_of(c("reward1", "dec")))
})

test_that("all eight families derive the groups the docs claim", {
  want <- list(
    bandit2arm_delta = c("reward1", "reward2"),
    bandit2arm_dual = c("reward1", "reward2"),
    prl_fictitious = c("reward1", "reward2"),
    rlddm = c("reward1", "reward2"),
    bandit4arm2_kalman_filter = paste0("payoff", 1:4),
    igt_pvl_delta = paste0("payoff", 1:4),
    igt_orl = paste0("payoff", 1:4),
    ts_par7 = paste0("payoff", 1:4))
  for (nm in names(want)) {
    fam <- do.call(nm, list(subject = quote(id), trial = quote(trial)))
    got <- fam[["learn"]][["counterfactual"]]
    expect_equal(got, list(want[[nm]]), info = nm)
  }
  # ts_par7() names stage2() as well and must not have picked it up
  ts <- ts_par7(subject = id, trial = trial)
  expect_false(any(grepl("stage2",
                         unlist(ts[["learn"]][["counterfactual"]]))))
})

## ---- the opt-out argument ------------------------------------------

# A family built out of a shipped one's own recursion, so that the
# argument can be exercised without inventing a learning rule.
ct_probe <- function(...) {
  lrn <- bandit2arm_delta(subject = id, trial = trial)[["learn"]]
  ln_family("probe", lrn[["subject_expr"]], lrn[["trial_expr"]],
            dpars = c("alpha", "tau"),
            links = list(alpha = "logit", tau = "log"),
            primary = "alpha",
            inits = list(alpha = function(y, aterms) 0.3,
                         tau = function(y, aterms) 1),
            aterms = c("reward1", "reward2"), spec = lrn[["spec"]],
            data_map = c(reward1 = "pay1", reward2 = "pay2"), ...)
}

test_that("counterfactual = TRUE is refused rather than obeyed", {
  # TRUE is the spelling a maintainer reaches for to mean "yes, guard
  # this one", and before this it turned the guard OFF. So did NA, 0,
  # "no" and list(). Each one now stops the family being built, which
  # is where a mistake in it can still be seen.
  for (bad in list(TRUE, NA, 0, "no", list(), c("reward1"),
                   list(character(0)), list(1:2))) {
    expect_error(ct_probe(counterfactual = bad), "counterfactual",
                 info = paste(class(bad), length(bad)))
  }
  expect_error(ct_probe(counterfactual = TRUE), "does NOT take TRUE",
               fixed = TRUE)
})

test_that("the three legal spellings all mean what they say", {
  # absent: derived from the terms the family names
  expect_equal(ct_probe()[["learn"]][["counterfactual"]],
               list(c("reward1", "reward2")))
  # identical(FALSE): the opt-out, and no shipped family takes it
  expect_null(ct_probe(counterfactual = FALSE)[["learn"]][[
    "counterfactual"]])
  # written out, as a vector or as a list of them
  expect_equal(
    ct_probe(counterfactual = c("reward1", "reward2"))[["learn"]][[
      "counterfactual"]],
    list(c("reward1", "reward2")))
  expect_equal(
    ct_probe(counterfactual = list(c("reward1", "reward2")))[[
      "learn"]][["counterfactual"]],
    list(c("reward1", "reward2")))
})

test_that("the opt-out really does open the draw, and only it", {
  # the guard is what stands between a duplicated schedule and a draw,
  # so the opt-out has to be shown to remove it rather than assumed to
  d <- ct_design()
  dup <- d
  dup$pay2 <- dup$pay1
  ps <- list(alpha = 0.4, tau = 3)
  expect_error(frm_task_simulate(ct_probe(), dup, pars = ps, seed = 2L),
               "reward\\(\\)")
  out <- frm_task_simulate(ct_probe(counterfactual = FALSE), dup,
                           pars = ps, seed = 2L)
  expect_equal(nrow(out[[1L]]), nrow(d))
})
