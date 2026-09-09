## Phase 0 of dev/extension-gaps-plan.md: the learn row.
##
## 100 subjects x 200 trials, correlated random effects on every
## parameter: bandit2arm_delta with `(1 | p | id)` on both of its
## parameters, then rlddm the same way. 20,000 rows.
##
## What the row decides: whether the per-trial recursion scales to the
## population designs the literature runs.
##
## Two choices this file made where the plan left them open.
##
## rlddm has five parameters and the plan says "the same way". The
## correlated block here covers alpha, drift, bs and ndt, and `bias` is
## held at 0.5. That is what the reinforcement-learning literature does
## with a start point between two ARMS rather than a correct and an
## error response, and ?rlddm says so; estimating it would add a fifth
## row and column to the correlation matrix for a parameter the design
## has no information about.
##
## The truth for the correlated block is a diagonal one, that is, a
## correlation of zero between the subject deviations. A non-zero truth
## would make the row a recovery study of the correlation, which is
## Phase 2's item 2.2, not Phase 0's cost measurement.
##
## See dev/scale-findings.md for the numbers this produced.

# Each standard deviation is stated ON THE LINK the family estimates
# that parameter on, so that the fitted component has a truth to be
# compared with: `alpha` is logit, `tau` and `bs` are log, `drift` is
# IDENTITY. The one exception is `ndt`, whose link is a scaled logit
# whose bound is the data's own fastest response. Its deviation is
# drawn on the LOG scale and no link-scale truth is stated for it,
# because there is none to state: the link's bound is a property of the
# drawn data. That exception is not a convenience; it is the thing the
# rlddm row is measuring.
learn_truth <- list(alpha = 0.35, tau = 3, sd_alpha = 0.5, sd_tau = 0.3,
                    drift = 2.5, bs = 1.5, ndt = 0.25, sd_drift = 1.0,
                    sd_bs = 0.2, sd_ndt = 0.15)

# One Wald interval, by the name confint() gives it.
learn_ci <- function(ci, key) {
  j <- grep(key, rownames(ci), fixed = TRUE)
  if (!length(j)) return(c(NA_real_, NA_real_))
  as.numeric(ci[j[1L], 1:2])
}

# The bandit design: 100 subjects x 200 trials of a stationary
# two-armed bandit, with a subject deviation on the learning rate and
# on the inverse temperature.
learn_bandit_data <- function(seed = 20260908L,
                              ns = if (scale_small()) 5L else 100L,
                              nt = if (scale_small()) 30L else 200L) {
  tr <- learn_truth
  d <- frm_task_design("bandit2arm", n_subject = ns, n_trial = nt,
                       seed = seed)
  set.seed(seed + 1L)
  ua <- stats::rnorm(ns, 0, tr$sd_alpha)
  ut <- stats::rnorm(ns, 0, tr$sd_tau)
  i <- as.integer(d$id)
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = id, trial = trial), d,
    pars = list(alpha = stats::plogis(stats::qlogis(tr$alpha) + ua[i]),
                tau = exp(log(tr$tau) + ut[i])),
    seed = seed)[[1L]]$choice
  d
}

# The rlddm design on the same task frame: the learning rule feeds a
# diffusion, so a trial carries a response TIME and the boundary it
# reached.
learn_rlddm_data <- function(seed = 20260908L,
                             ns = if (scale_small()) 5L else 100L,
                             nt = if (scale_small()) 30L else 200L) {
  tr <- learn_truth
  d <- frm_task_design("bandit2arm", n_subject = ns, n_trial = nt,
                       seed = seed)
  set.seed(seed + 2L)
  i <- as.integer(d$id)
  ua <- stats::rnorm(ns, 0, tr$sd_alpha)
  ud <- stats::rnorm(ns, 0, tr$sd_drift)
  ub <- stats::rnorm(ns, 0, tr$sd_bs)
  un <- stats::rnorm(ns, 0, tr$sd_ndt)
  s <- frm_task_simulate(
    rlddm(subject = id, trial = trial), d,
    pars = list(alpha = stats::plogis(stats::qlogis(tr$alpha) + ua[i]),
                drift = tr$drift + ud[i],
                bs = tr$bs * exp(ub[i]),
                ndt = tr$ndt * exp(un[i]),
                bias = 0.5),
    seed = seed)[[1L]]
  s
}

learn_scale_run <- function(row, form, fam, d, truth_key, truth_value,
                            ci_key, sd_truth) {
  scale_mem_reset()
  bd <- scale_build(form, family = fam, data = d)
  g0 <- scale_grad(bd$dry$obj, bd$dry$obj$par)
  ctl <- scale_control(bd$dry$obj, bd$dry$obj$par, g0$calls)
  bd$dry <- NULL

  fit <- NULL
  t_fit <- scale_elapsed(fit <- frm(form, family = fam, data = d,
                                    se = TRUE))
  g1 <- scale_grad(fit$obj, fit$opt$par)
  mem <- scale_mem_peak_mb()

  ci <- suppressWarnings(stats::confint(fit))
  b <- unlist(fixef(fit))
  vc <- VarCorr(fit)
  sds <- sqrt(diag(vc[[1L]]))
  # the whole correlated block, by name, because the question this row
  # asks is whether a correlated block on EVERY parameter survives the
  # design, and a component that collapsed to zero is the answer
  sd_all <- paste(rownames(vc[[1L]]),
                  formatC(sds, digits = 4, format = "g"),
                  sep = "=", collapse = ";")
  cr <- stats::cov2cor(vc[[1L]])
  cor_max <- max(abs(cr[lower.tri(cr)]))
  i_k <- learn_ci(ci, ci_key)
  scale_record(
    row, rows = nrow(d), n_par = length(fit$opt$par),
    build_s = bd$build_s, frame_s = bd$frame_s,
    grad_start_s = g0$seconds, grad_start_calls = g0$calls,
    grad_opt_s = g1$seconds, control_ratio = ctl,
    fit_s = t_fit, mem_mb = mem,
    logLik = as.numeric(stats::logLik(fit)),
    key = truth_key, est = unname(b[truth_key]), truth = truth_value,
    est_lo = i_k[1L], est_hi = i_k[2L],
    sd_all = sd_all, sd_true = sd_truth, n_sd = length(sds),
    cor_max_abs = cor_max, cor_true = 0,
    diag = scale_diag(fit))
  list(fit = fit, ci = ci, interval = i_k)
}

test_that("the learn bandit scale row fits and reports its cost", {
  skip_unless_scale()
  scale_row_on("learn")
  d <- learn_bandit_data()
  form <- bf(choice | reward(pay1, pay2) ~ 1 + (1 | p | id),
             tau ~ 1 + (1 | p | id))
  r <- learn_scale_run("learn", form,
                       bandit2arm_delta(subject = id, trial = trial), d,
                       "alpha.(Intercept)",
                       stats::qlogis(learn_truth$alpha), "(Intercept)",
                       "alpha=0.5;tau=0.3")
  expect_true(is.finite(as.numeric(stats::logLik(r$fit))))
  tr <- stats::qlogis(learn_truth$alpha)
  expect_true(r$interval[1L] <= tr && r$interval[2L] >= tr)
})

test_that("the learn rlddm scale row fits and reports its cost", {
  skip_unless_scale()
  scale_row_on("learn-rlddm")
  d <- learn_rlddm_data()
  form <- bf(rt | dec(choice) + reward(pay1, pay2) ~ 1 + (1 | p | id),
             drift ~ 1 + (1 | p | id), bs ~ 1 + (1 | p | id),
             ndt ~ 1 + (1 | p | id), bias = 0.5)
  r <- learn_scale_run("learn-rlddm", form,
                       rlddm(subject = id, trial = trial), d,
                       "alpha.(Intercept)",
                       stats::qlogis(learn_truth$alpha), "(Intercept)",
                       "alpha=0.5;drift=1;bs=0.2;ndt=no link truth")
  expect_true(is.finite(as.numeric(stats::logLik(r$fit))))
})
