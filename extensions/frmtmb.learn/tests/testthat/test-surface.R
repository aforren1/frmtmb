## The rest of the frmtmb surface on a learning fit. A family is not
## finished when it fits; it is finished when the methods a user reaches
## for next either work or refuse for a stated reason.

ln_surface_fit <- function(seed = 61L) {
  d <- frm_task_design("reversal", n_subject = 12L, n_trial = 50L,
                       seed = seed)
  set.seed(seed)
  a <- stats::plogis(stats::qlogis(0.35) + stats::rnorm(12L, 0, 0.6))
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = id, trial = trial), d,
    pars = list(alpha = a, tau = 3), seed = seed)[[1L]]$choice
  list(d = d,
       fit = frmtmb::frm(
         frmtmb::bf(choice | reward(pay1, pay2) ~ after_reversal + (1 | id),
                    tau ~ 1),
         family = bandit2arm_delta(subject = id, trial = trial), data = d))
}

test_that("summary, fixef, ranef and logLik work", {
  o <- ln_surface_fit()
  f <- o$fit
  expect_s3_class(f, "frmtmb_fit")
  expect_no_error(summary(f))
  fx <- frmtmb::fixef(f)
  expect_setequal(names(fx), c("alpha", "tau"))
  expect_true(is.finite(as.numeric(stats::logLik(f))))
  expect_true(is.finite(stats::AIC(f)))
  expect_equal(stats::nobs(f), nrow(o$d))
  expect_gte(length(frmtmb::ranef(f)), 1L)
  expect_true(is.finite(frmtmb::VarCorr(f)[[1L]][1L, 1L]))
})

test_that("predict gives the learning parameters, and a mean on the data", {
  o <- ln_surface_fit()
  f <- o$fit
  p <- stats::predict(f, type = "link")
  expect_length(p, nrow(o$d))
  expect_true(all(is.finite(p)))
  for (dp in c("alpha", "tau")) {
    v <- stats::predict(f, type = "link", dpar = dp)
    expect_true(all(is.finite(v)))
  }
  # the link scale reaches new data, because a linear predictor is
  # rowwise and belongs to the core
  nd <- o$d[1:5, ]
  expect_true(all(is.finite(
    stats::predict(f, newdata = nd, type = "link", dpar = "alpha"))))
  # The response scale is refused, and that is the design rather than a
  # gap: a nominal option code has no mean, so core has nothing to
  # return and says so in the family's own name.
  expect_error(stats::predict(f, type = "response"), "declares no mean")
  expect_error(stats::fitted(f), "bandit2arm_delta")
  # what replaces it. This fit is hierarchical, so sum(log(p)) is the
  # CONDITIONAL data log-likelihood and sits ABOVE the Laplace marginal
  # logLik() reports; the exact identity is asserted on a fit with no
  # random effects in test-engine.R.
  tr <- frm_value_trace(f)
  expect_length(tr$p, nrow(o$d))
  expect_true(all(tr$p > 0 & tr$p < 1))
  expect_gt(sum(log(tr$p)), as.numeric(stats::logLik(f)))
})

test_that("every residual type refuses, each for its own stated reason", {
  f <- ln_surface_fit()$fit
  # 'response' and 'pearson' refuse first, on the missing mean
  expect_error(stats::residuals(f, type = "response"), "declares no mean")
  expect_error(stats::residuals(f, type = "pearson"), "declares no mean")
  # 'deviance' and 'osa' refuse earlier still, in this package's own
  # words, for reasons that would hold even if a mean existed
  expect_error(stats::residuals(f, type = "deviance"),
               "loglik slot returns one total")
  expect_error(stats::residuals(f, type = "osa"),
               "no registered observation vector")
})

test_that("conditional_effects and importance refuse by name", {
  o <- ln_surface_fit()
  expect_error(frmtmb::conditional_effects(o$fit),
               "synthetic grid this function builds does not have")
  expect_error(
    frmtmb::frm(
      frmtmb::bf(choice | reward(pay1, pay2) ~ after_reversal + (1 | id),
                 tau ~ 1),
      family = bandit2arm_delta(subject = id, trial = trial), data = o$d,
      importance = 50),
    "importance")
})

test_that("the reshaping addition terms refuse in the family's words", {
  o <- ln_surface_fit()
  d <- o$d
  d$w <- 1
  expect_error(
    frmtmb::frm(
      frmtmb::bf(choice | reward(pay1, pay2) + weights(w) ~ 1, tau ~ 1),
      family = bandit2arm_delta(subject = id, trial = trial), data = d),
    "reshapes a per-row likelihood contribution")
})

test_that("simulate works, and refuses for the family that cannot", {
  o <- ln_surface_fit()
  s <- stats::simulate(o$fit, nsim = 2, seed = 3)
  expect_length(s, 2L)
  expect_true(all(unlist(s) %in% c(1L, 2L)))

  dt <- frm_task_design("twostep", n_subject = 5L, n_trial = 30L, seed = 62L)
  dt <- frm_task_simulate(
    ts_par7(subject = id, trial = trial), dt,
    pars = list(w = 0.5, alpha1 = 0.4, tau1 = 3, alpha2 = 0.4, tau2 = 3,
                lambda = 0.5, pers = 0), seed = 62)[[1L]]
  ft <- frmtmb::frm(
    frmtmb::bf(choice | stage2(state2, choice2) +
                 payoff(pay1, pay2, pay3, pay4) ~ 1,
               alpha1 ~ 1, tau1 ~ 1, alpha2 ~ 1, tau2 ~ 1, lambda ~ 1,
               pers ~ 1),
    family = ts_par7(subject = id, trial = trial), data = dt)
  expect_error(stats::simulate(ft, nsim = 1),
               "stage-two state the environment answers with")
})

test_that("par_template and set_prior reach every learning parameter", {
  o <- ln_surface_fit()
  tmpl <- frmtmb::par_template(o$fit)
  nms <- names(unlist(tmpl))
  expect_true(length(nms) > 0)
  gp <- frmtmb::get_prior(
    frmtmb::bf(choice | reward(pay1, pay2) ~ after_reversal, tau ~ 1),
    family = bandit2arm_delta(subject = id, trial = trial), data = o$d)
  expect_true(is.data.frame(gp))
  # the PRIMARY distributional parameter's rows carry an empty dpar,
  # which is core's spelling for "the main linear predictor"; every
  # other parameter is named
  expect_true("tau" %in% gp$dpar)
  expect_true(any(!nzchar(gp$dpar)))
  # a prior on the learning rate's slope shrinks it
  f2 <- frmtmb::frm(
    frmtmb::bf(choice | reward(pay1, pay2) ~ after_reversal + (1 | id),
               tau ~ 1),
    family = bandit2arm_delta(subject = id, trial = trial), data = o$d,
    prior = frmtmb::set_prior("normal(0, 0.05)", class = "b",
                              dpar = "alpha"))
  a1 <- abs(unlist(frmtmb::fixef(f2))[["alpha.after_reversalafter"]])
  a0 <- abs(unlist(frmtmb::fixef(o$fit))[["alpha.after_reversalafter"]])
  expect_lt(a1, a0)
})

test_that("refit rebuilds the block and reproduces the fit", {
  o <- ln_surface_fit()
  # refit() takes a new response; handing it the observed one must
  # reproduce the original fit, which is what says the block was
  # rebuilt correctly rather than carried over stale
  r <- frmtmb::refit(o$fit, newresp = o$d$choice)
  expect_equal(as.numeric(stats::logLik(r)),
               as.numeric(stats::logLik(o$fit)), tolerance = 1e-6)
  # and a refit on a DRAW is the recovery study's inner loop
  s <- stats::simulate(o$fit, nsim = 1, seed = 5)[[1L]]
  r2 <- frmtmb::refit(o$fit, newresp = s)
  expect_true(is.finite(as.numeric(stats::logLik(r2))))
})

test_that("frm_value_trace refuses a fit it cannot replay", {
  set.seed(63)
  d <- data.frame(y = stats::rnorm(40), x = stats::runif(40))
  f <- frmtmb::frm(frmtmb::bf(y ~ x), family = stats::gaussian(), data = d)
  expect_error(frm_value_trace(f), "no value recursion to replay")
})

test_that("the compatibility rows this package registered resolve", {
  fams <- c("bandit2arm_delta", "bandit2arm_dual", "prl_fictitious",
            "bandit4arm2_kalman_filter", "ts_par7", "igt_pvl_delta")
  feats <- frmtmb::frm_compat_features()
  expect_true(all(fams %in% feats$key))
  expect_true(all(c("frm_value_trace", "frm_task_simulate") %in% feats$key))
  # the three addition terms joined the vocabulary with their registration
  expect_true(all(c("reward", "payoff", "stage2") %in%
                    feats$key[feats$kind == "aterm"]))
  ok <- c("works", "conditional", "refused", "untested", "broken")
  for (nm in fams) {
    for (ft in c("importance", "cens()", "fitted", "predict", "residuals",
                 "s()", "us", "simulate", "prior")) {
      r <- frmtmb::frm_compat(nm, ft)
      st <- if (is.data.frame(r)) r$status[1L] else r$status
      expect_true(st %in% ok, info = paste(nm, ft))
    }
  }
  # the seam is named by name in the row it belongs to
  imp <- frmtmb::frm_compat("bandit2arm_delta", "importance")
  expect_identical(imp$status, "refused")
  expect_match(imp$note, "finest factorization")
  # simulate is the one capability that differs between the families
  expect_identical(frmtmb::frm_compat("bandit2arm_delta", "simulate")$status,
                   "works")
  expect_identical(frmtmb::frm_compat("ts_par7", "simulate")$status,
                   "refused")
  expect_identical(frmtmb::frm_compat("igt_pvl_delta", "fitted")$status,
                   "refused")
  expect_match(frmtmb::frm_compat("igt_pvl_delta", "fitted")$note,
               "NOMINAL")
})

test_that("the reference table matches the families that exist", {
  d <- frm_learn_families()
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 6L)
  for (nm in d$family) {
    expect_true(exists(nm, envir = asNamespace("frmtmb.learn")),
                info = nm)
  }
  # and the parameter lists in it are the families' own
  mk <- list(bandit2arm_delta = bandit2arm_delta(subject = id),
             bandit2arm_dual = bandit2arm_dual(subject = id),
             prl_fictitious = prl_fictitious(subject = id),
             bandit4arm2_kalman_filter =
               bandit4arm2_kalman_filter(subject = id),
             ts_par7 = ts_par7(subject = id),
             igt_pvl_delta = igt_pvl_delta(subject = id))
  for (i in seq_len(nrow(d))) {
    listed <- trimws(strsplit(d$pars[i], ",")[[1L]])
    expect_setequal(listed, mk[[d$family[i]]]$learn$dpars)
    # the hBayesDM column has one entry per parameter, in the same order
    expect_length(trimws(strsplit(d$hbayesdm_pars[i], ",")[[1L]]),
                  length(listed))
  }
})
