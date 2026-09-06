## Every family against an independent Stan program of the same model,
## at frmtmb's own estimates. See helper-stan.R for why this is an
## identity rather than an agreement, and helper-stan-programs.R for the
## programs.
##
## The fixtures all have the same shape: real between-subject spread in
## the primary parameter, a condition effect on it, and one random
## intercept whose standard deviation the Stan program takes as data.
## Without real spread the fit collapses to sd = 0, the subject effects
## are all zero, and the gradient check asserts nothing.

ln_fix_2arm <- function(fam, task = "bandit2arm", ns = 12L, nt = 60L,
                        seed = 101L, pars = NULL, dpar = "alpha",
                        rest = list(tau ~ 1), sign_pay = FALSE,
                        graded = FALSE) {
  set.seed(seed)
  d <- frm_task_design(task, n_subject = ns, n_trial = nt, seed = seed)
  if (graded) {
    # GRADED payoffs, and for one family they are load-bearing. On
    # binary payoffs the value store stays inside the payoff range, so
    # the sign of the prediction error and the sign of the outcome agree
    # on every trial and bandit2arm_dual's two splits are provably the
    # same model. A binary fixture would therefore give the `outcome`
    # row of the identity table numbers identical to the `pe` row, and
    # the `outcome` branch of ln_stan_code_dual() would never be
    # exercised where it differs. test-families.R makes the same point
    # from the R side.
    d$pay1 <- stats::runif(nrow(d), -1, 1)
    d$pay2 <- stats::runif(nrow(d), -1, 1)
  } else if (sign_pay) {
    d$pay1 <- 2 * d$pay1 - 1
    d$pay2 <- 2 * d$pay2 - 1
  }
  d$cond <- factor(rep(c("a", "b"), length.out = ns)[as.integer(d$id)])
  cb <- as.numeric("b" == as.character(d$cond[!duplicated(d$id)]))
  u <- stats::rnorm(ns, 0, 0.6)
  pars[[dpar]] <- stats::plogis(stats::qlogis(0.35) + 0.9 * cb + u)
  d$choice <- frm_task_simulate(fam, d, pars = pars, seed = seed)[[1L]]$choice
  form <- do.call(frmtmb::bf,
                  c(list(choice | reward(pay1, pay2) ~ cond + (1 | id)),
                    rest))
  list(d = d, fit = frmtmb::frm(form, family = fam, data = d))
}

test_that("bandit2arm_delta reproduces a Stan program of the same model", {
  skip_unless_stan()
  o <- ln_fix_2arm(bandit2arm_delta(subject = id, trial = trial),
                   pars = list(tau = 3))
  r <- ln_lp_check(o$fit, ln_stan_code_delta(),
                   ln_stan_data(o$fit, o$d, ~ cond,
                                list(pay1 = as.numeric(o$d$pay1),
                                     pay2 = as.numeric(o$d$pay2))),
                   label = "bandit2arm_delta")
  # a deliberately non-stationary point, so the agreement cannot be an
  # artifact of both sides sitting at an optimum
  p2 <- o$fit$obj$env$last.par.best
  p2[names(p2) == "b"] <- p2[names(p2) == "b"] + 0.3
  r2 <- ln_lp_check(o$fit, ln_stan_code_delta(),
                    ln_stan_data(o$fit, o$d, ~ cond,
                                 list(pay1 = as.numeric(o$d$pay1),
                                      pay2 = as.numeric(o$d$pay2))),
                    par = p2, check_grad = FALSE,
                    label = "bandit2arm_delta (displaced)")
  expect_lt(abs(r2$const), 1e-6 * max(1, abs(r2$ours)))
  expect_false(isTRUE(all.equal(r$ours, r2$ours)))
})

test_that("bandit2arm_dual reproduces Stan under both splits", {
  skip_unless_stan()
  got <- list()
  for (sp in c("pe", "outcome")) {
    fam <- bandit2arm_dual(subject = id, trial = trial, split = sp)
    o <- ln_fix_2arm(fam, seed = 102L, dpar = "Arew",
                     pars = list(Apun = 0.15, tau = 3),
                     rest = list(Apun ~ 1, tau ~ 1), graded = TRUE)
    # The pe split selects its rate with sign(pe), so on GRADED payoffs
    # the joint log density has a kink in the random effects wherever a
    # prediction error crosses zero, and TMB's inner Newton solve stops
    # short of the conditional mode: measured 3.97e-02 here against
    # 1.44e-15 for the outcome split on the same data, whose selector is
    # data and therefore smooth. That is a property of the model, not of
    # the comparison. The VALUE identity is unaffected, because both
    # sides are evaluated at the same point, and it holds at 5.7e-14.
    tg <- if (sp == "pe") 0.05 else 1e-4
    r <- ln_lp_check(o$fit, ln_stan_code_dual(sp),
                     ln_stan_data(o$fit, o$d, ~ cond,
                                  list(pay1 = as.numeric(o$d$pay1),
                                       pay2 = as.numeric(o$d$pay2))),
                     tol_grad = tg,
                     label = paste0("bandit2arm_dual (", sp, ")"))
    got[[sp]] <- r$ours
  }
  # the two splits are DIFFERENT models on graded payoffs, so the two
  # rows of the identity table are two checks rather than one repeated
  expect_false(isTRUE(all.equal(got[["pe"]], got[["outcome"]])))
})

test_that("prl_fictitious reproduces a Stan program of the same model", {
  skip_unless_stan()
  o <- ln_fix_2arm(prl_fictitious(subject = id, trial = trial),
                   task = "reversal", seed = 103L,
                   pars = list(bias = 0.2, tau = 3),
                   rest = list(bias ~ 1, tau ~ 1), sign_pay = TRUE)
  ln_lp_check(o$fit, ln_stan_code_fict(),
              ln_stan_data(o$fit, o$d, ~ cond,
                           list(pay1 = as.numeric(o$d$pay1),
                                pay2 = as.numeric(o$d$pay2))),
              label = "prl_fictitious")
})

test_that("bandit4arm2_kalman_filter reproduces a Stan program", {
  skip_unless_stan()
  set.seed(104)
  ns <- 10L
  fam <- bandit4arm2_kalman_filter(subject = id, trial = trial)
  d <- frm_task_design("bandit4arm_restless", n_subject = ns, n_trial = 50,
                       seed = 104)
  d$cond <- factor(rep(c("a", "b"), length.out = ns)[as.integer(d$id)])
  cb <- as.numeric("b" == as.character(d$cond[!duplicated(d$id)]))
  tau <- exp(log(0.15) + 0.4 * cb + stats::rnorm(ns, 0, 0.3))
  d$choice <- frm_task_simulate(
    fam, d, pars = list(tau = tau, lambda = 0.98, center = 50, mu0 = 50,
                        sigma0 = 10, sigmaD = 3), seed = 104)[[1L]]$choice
  fit <- frmtmb::frm(
    frmtmb::bf(choice | payoff(pay1, pay2, pay3, pay4) ~ cond + (1 | id),
               lambda ~ 1, center ~ 1, mu0 ~ 1, sigma0 ~ 1, sigmaD ~ 1),
    family = fam, data = d)
  ln_lp_check(fit, ln_stan_code_kalman(),
              ln_stan_data(fit, d, ~ cond,
                           list(pay = as.matrix(d[, paste0("pay", 1:4)]),
                                sigma_o = 4)),
              label = "bandit4arm2_kalman_filter")
})

test_that("igt_pvl_delta reproduces a Stan program of the same model", {
  skip_unless_stan()
  set.seed(105)
  ns <- 12L
  fam <- igt_pvl_delta(subject = id, trial = trial)
  d <- frm_task_design("igt", n_subject = ns, n_trial = 60, seed = 105)
  d$cond <- factor(rep(c("a", "b"), length.out = ns)[as.integer(d$id)])
  cb <- as.numeric("b" == as.character(d$cond[!duplicated(d$id)]))
  al <- stats::plogis(stats::qlogis(0.3) + 0.8 * cb + stats::rnorm(ns, 0, 0.5))
  d$choice <- frm_task_simulate(
    fam, d, pars = list(alpha = al, shape = 0.4, lambda = 1.5, tau = 1),
    seed = 105)[[1L]]$choice
  fit <- frmtmb::frm(
    frmtmb::bf(choice | payoff(pay1, pay2, pay3, pay4) ~ cond + (1 | id),
               shape ~ 1, lambda ~ 1, tau ~ 1),
    family = fam, data = d)
  ln_lp_check(fit, ln_stan_code_pvl(),
              ln_stan_data(fit, d, ~ cond,
                           list(pay = as.matrix(d[, paste0("pay", 1:4)]))),
              label = "igt_pvl_delta")
})

test_that("ts_par7 reproduces a Stan program of the same model", {
  skip_unless_stan()
  set.seed(106)
  ns <- 12L
  fam <- ts_par7(subject = id, trial = trial)
  d <- frm_task_design("twostep", n_subject = ns, n_trial = 60, seed = 106)
  d$cond <- factor(rep(c("a", "b"), length.out = ns)[as.integer(d$id)])
  cb <- as.numeric("b" == as.character(d$cond[!duplicated(d$id)]))
  w <- stats::plogis(0.9 * cb + stats::rnorm(ns, 0, 0.6))
  d <- frm_task_simulate(
    fam, d, pars = list(w = w, alpha1 = 0.4, tau1 = 3, alpha2 = 0.4,
                        tau2 = 3, lambda = 0.6, pers = 0.2),
    seed = 106)[[1L]]
  fit <- frmtmb::frm(
    frmtmb::bf(choice | stage2(state2, choice2) +
                 payoff(pay1, pay2, pay3, pay4) ~ cond + (1 | id),
               alpha1 ~ 1, tau1 ~ 1, alpha2 ~ 1, tau2 ~ 1, lambda ~ 1,
               pers ~ 1),
    family = fam, data = d)
  ln_lp_check(fit, ln_stan_code_ts(),
              ln_stan_data(fit, d, ~ cond,
                           list(pay = as.matrix(d[, paste0("pay", 1:4)]),
                                s2 = as.integer(d$state2),
                                c2 = as.integer(d$choice2), pc = 0.7)),
              label = "ts_par7")
})
