## Every family against the model written out longhand.
##
## The engine vectorizes across subjects, masks padded trials and
## selects the chosen option by multiplying by an indicator. Each
## reference below does the opposite of all three: one subject and one
## trial at a time, no mask, `if` and `[` for selection. It is the same
## spelling the Stan programs use, so a bug shared between this and the
## engine would have to be a bug in the model statement rather than in
## the vectorization.
##
## The comparison point is sum(log(frm_value_trace(fit)$p)), on fits
## with NO random effects, where that equals logLik(fit) exactly.

ln_sub_rows <- function(d) {
  lapply(levels(d$id), function(s) {
    r <- which(d$id == s)
    r[order(d$trial[r])]
  })
}

ln_soft <- function(u, k) u[k] - log(sum(exp(u)))

ln_ref_2arm <- function(d, step, init = c(0, 0)) {
  ll <- 0
  for (rows in ln_sub_rows(d)) {
    q <- init
    for (i in rows) {
      ll <- ll + step$ll(q, d, i)
      q <- step$up(q, d, i)
    }
  }
  ll
}

ln_ref_kalman <- function(d, p, sigma_o = 4) {
  ll <- 0
  for (rows in ln_sub_rows(d)) {
    mu <- rep(p$mu0, 4)
    v <- rep(p$sigma0^2, 4)
    for (i in rows) {
      ll <- ll + ln_soft(p$tau * mu, d$choice[i])
      k <- d$choice[i]
      r <- d[[paste0("pay", k)]][i]
      g <- v[k] / (v[k] + sigma_o^2)
      mu[k] <- mu[k] + g * (r - mu[k])
      v[k] <- v[k] * (1 - g)
      mu <- p$lambda * mu + (1 - p$lambda) * p$center
      v <- p$lambda^2 * v + p$sigmaD^2
    }
  }
  ll
}

ln_ref_pvl <- function(d, p) {
  ll <- 0
  for (rows in ln_sub_rows(d)) {
    q <- rep(0, 4)
    for (i in rows) {
      ll <- ll + ln_soft(p$tau * q, d$choice[i])
      k <- d$choice[i]
      x <- d[[paste0("pay", k)]][i]
      mag <- if (x == 0) 0 else abs(x)^p$shape
      ut <- if (x >= 0) mag else -p$lambda * mag
      q[k] <- q[k] + p$alpha * (ut - q[k])
    }
  }
  ll
}

ln_ref_ts <- function(d, p, pc = 0.7) {
  ll <- 0
  for (rows in ln_sub_rows(d)) {
    qmf <- c(0, 0)
    q2 <- rep(0, 4)
    rp <- c(0, 0)
    for (i in rows) {
      b1 <- max(q2[1], q2[2])
      b2 <- max(q2[3], q2[4])
      mb <- c(pc * b1 + (1 - pc) * b2, (1 - pc) * b1 + pc * b2)
      u1 <- p$tau1 * (p$w * mb + (1 - p$w) * qmf + p$pers * rp)
      ll <- ll + ln_soft(u1, d$choice[i])
      base <- 2 * (d$state2[i] - 1)
      u2 <- p$tau2 * q2[base + 1:2]
      ll <- ll + ln_soft(u2, d$choice2[i])
      m <- base + d$choice2[i]
      d1 <- q2[m] - qmf[d$choice[i]]
      d2 <- d[[paste0("pay", m)]][i] - q2[m]
      qmf[d$choice[i]] <- qmf[d$choice[i]] + p$alpha1 * (d1 + p$lambda * d2)
      q2[m] <- q2[m] + p$alpha2 * d2
      rp <- c(d$choice[i] == 1, d$choice[i] == 2) * 1
    }
  }
  ll
}

# the fitted parameters of a no-random-effect fit, on natural scales
ln_nat <- function(fit) {
  b <- unlist(frmtmb::fixef(fit))
  nm <- sub("[.][(]Intercept[)]$", "", names(b))
  fam <- frmtmb::single_response(fit, "a fit")$family
  lk <- fam$links
  out <- list()
  for (i in seq_along(b)) {
    # the family's own link object, applied here rather than a switch on
    # its name, so the reference cannot disagree with the fit about what
    # "logit" means
    out[[nm[i]]] <- lk[[nm[i]]]$linkinv(b[[i]])
  }
  out
}

ln_engine_ll <- function(fit) sum(log(frm_value_trace(fit)$p))

test_that("bandit2arm_delta matches the model written out longhand", {
  d <- frm_task_design("bandit2arm", n_subject = 6L, n_trial = 40L,
                       seed = 71L)
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = id, trial = trial), d,
    pars = list(alpha = 0.4, tau = 3), seed = 71)[[1L]]$choice
  fit <- frmtmb::frm(
    frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
    family = bandit2arm_delta(subject = id, trial = trial), data = d)
  p <- ln_nat(fit)
  ref <- ln_ref_2arm(d, list(
    ll = function(q, d, i) ln_soft(p$tau * q, d$choice[i]),
    up = function(q, d, i) {
      k <- d$choice[i]
      pay <- if (k == 1) d$pay1[i] else d$pay2[i]
      q[k] <- q[k] + p$alpha * (pay - q[k])
      q
    }))
  expect_equal(ref, ln_engine_ll(fit), tolerance = 1e-9)
  expect_equal(ref, as.numeric(stats::logLik(fit)), tolerance = 1e-9)
})

test_that("prl_fictitious matches the model written out longhand", {
  d <- frm_task_design("reversal", n_subject = 6L, n_trial = 40L,
                       seed = 72L)
  d$pay1 <- 2 * d$pay1 - 1
  d$pay2 <- 2 * d$pay2 - 1
  d$choice <- frm_task_simulate(
    prl_fictitious(subject = id, trial = trial), d,
    pars = list(alpha = 0.3, bias = 0.2, tau = 3), seed = 72)[[1L]]$choice
  fit <- frmtmb::frm(
    frmtmb::bf(choice | reward(pay1, pay2) ~ 1, bias ~ 1, tau ~ 1),
    family = prl_fictitious(subject = id, trial = trial), data = d)
  p <- ln_nat(fit)
  ref <- ln_ref_2arm(d, list(
    ll = function(q, d, i) ln_soft(c(p$tau * q[1] + p$bias, p$tau * q[2]),
                                   d$choice[i]),
    up = function(q, d, i) {
      oc <- if (d$choice[i] == 1) d$pay1[i] else d$pay2[i]
      t1 <- if (d$choice[i] == 1) oc else -oc
      c(q[1] + p$alpha * (t1 - q[1]), q[2] + p$alpha * (-t1 - q[2]))
    }))
  expect_equal(ref, ln_engine_ll(fit), tolerance = 1e-9)
})

test_that("bandit4arm2_kalman_filter matches the model written out longhand", {
  d <- frm_task_design("bandit4arm_restless", n_subject = 6L, n_trial = 40L,
                       seed = 73L)
  fam <- bandit4arm2_kalman_filter(subject = id, trial = trial)
  d$choice <- frm_task_simulate(
    fam, d, pars = list(tau = 0.15, lambda = 0.98, center = 50, mu0 = 50,
                        sigma0 = 10, sigmaD = 3), seed = 73)[[1L]]$choice
  fit <- frmtmb::frm(
    frmtmb::bf(choice | payoff(pay1, pay2, pay3, pay4) ~ 1, lambda ~ 1,
               center ~ 1, mu0 ~ 1, sigma0 ~ 1, sigmaD ~ 1),
    family = fam, data = d)
  ref <- ln_ref_kalman(d, ln_nat(fit))
  expect_equal(ref, ln_engine_ll(fit), tolerance = 1e-8)
})

test_that("igt_pvl_delta matches the model written out longhand", {
  d <- frm_task_design("igt", n_subject = 6L, n_trial = 50L, seed = 74L)
  fam <- igt_pvl_delta(subject = id, trial = trial)
  d$choice <- frm_task_simulate(
    fam, d, pars = list(alpha = 0.3, shape = 0.4, lambda = 1.5, tau = 1),
    seed = 74)[[1L]]$choice
  fit <- frmtmb::frm(
    frmtmb::bf(choice | payoff(pay1, pay2, pay3, pay4) ~ 1, shape ~ 1,
               lambda ~ 1, tau ~ 1),
    family = fam, data = d)
  ref <- ln_ref_pvl(d, ln_nat(fit))
  expect_equal(ref, ln_engine_ll(fit), tolerance = 1e-8)
})

test_that("ts_par7 matches the model written out longhand", {
  d <- frm_task_design("twostep", n_subject = 12L, n_trial = 80L, seed = 75L)
  fam <- ts_par7(subject = id, trial = trial)
  d <- frm_task_simulate(
    fam, d, pars = list(w = 0.5, alpha1 = 0.4, tau1 = 3, alpha2 = 0.4,
                        tau2 = 3, lambda = 0.6, pers = 0.2),
    seed = 75)[[1L]]
  fit <- frmtmb::frm(
    frmtmb::bf(choice | stage2(state2, choice2) +
                 payoff(pay1, pay2, pay3, pay4) ~ 1,
               alpha1 ~ 1, tau1 ~ 1, alpha2 ~ 1, tau2 ~ 1, lambda ~ 1,
               pers ~ 1),
    family = fam, data = d)
  ref <- ln_ref_ts(d, ln_nat(fit))
  expect_equal(ref, ln_engine_ll(fit), tolerance = 1e-8)
})
