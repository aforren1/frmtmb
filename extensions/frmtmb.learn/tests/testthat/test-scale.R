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
## ITEM 1.0b. The rlddm row now carries `ndt_group(id)`. Through 0.3.0
## `ndt`'s bound was the fastest response in the WHOLE data set, which
## with 100 learners is the fastest response of all of them, so a
## subject deviation on `ndt` was a deviation on a fraction of somebody
## else's floor. This row asserted the consequence: conv=1, a maximum
## gradient of 6.36e+09, a Hessian that was not positive definite and
## four NaN standard errors at a log-likelihood of -8333.1. With the
## grouping each row is bounded by its own learner's fastest response
## and the assertions below are the recovery ones.
##
## See dev/scale-findings.md for the numbers this produced, and
## dev/rlddm-findings.md for the before and after.

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
  # The per-subject non-decision times the draw used, so that the fit's
  # own spread can be compared with them on the natural scale. The
  # deviation is drawn on the LOG scale just above; it is the recorded
  # TRUTH that is natural, and that is the part that matters, because
  # `ndt` is estimated as a fraction of a bound the drawn data sets.
  attr(s, "ndt_subject") <- tr$ndt * exp(un)
  # THE FOUR DEVIATIONS THE MODEL ACTUALLY ESTIMATES, which are not the
  # four the design drew. Three of them are: `alpha` on the logit,
  # `drift` on the identity, `bs` on the log. The fourth is not. Under a
  # per-group bound the model estimates `qlogis(ndt_i / floor_i)`, and
  # `floor_i` is that learner's own fastest response, which is a draw
  # from its WHOLE parameter vector: a learner with a higher boundary
  # separation is slower, so its floor is later, so it needs a LOWER
  # fraction to express the same non-decision time. The correlation
  # between the `bs` and `ndt` deviations is therefore negative BEFORE
  # any fit, and the design's diagonal truth is not the truth of the
  # block this row estimates. dev/rlddm-findings.md and the review's
  # `dev/rev-rlddm-cormax.R` carry the measurement.
  #
  # Recorded as a matrix rather than as a number so that a reader can
  # see which pair carries it, and computed from the truths and the
  # floors alone, with no fit.
  fl <- as.numeric(tapply(s$rt, s$id, min))[match(levels(s$id),
                                                  levels(s$id))]
  attr(s, "dev_drawn") <- cbind(alpha = ua, drift = ud, bs = ub,
                                ndt = un)
  attr(s, "dev_fitted_par") <- cbind(
    alpha = ua, drift = ud, bs = ub,
    ndt = stats::qlogis(pmin(pmax(attr(s, "ndt_subject") / fl, 1e-8),
                             1 - 1e-8)))
  attr(s, "own_floor") <- fl
  s
}

# The maximum absolute off-diagonal correlation of a set of deviations,
# which is the quantity `cor_max_abs` reports for the fitted block.
learn_cor_max <- function(m) {
  cr <- stats::cor(m)
  max(abs(cr[lower.tri(cr)]))
}

learn_scale_run <- function(row, form, fam, d, truth_key, truth_value,
                            ci_key, sd_truth, extra = NULL,
                            cor_true = 0) {
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
  # A row with something of its own to say says it here, so that both
  # rows still go through one runner and one record.
  ex <- if (is.null(extra)) list() else extra(fit, d)
  do.call(scale_record, c(list(
    row, rows = nrow(d), n_par = length(fit$opt$par),
    build_s = bd$build_s, frame_s = bd$frame_s,
    grad_start_s = g0$seconds, grad_start_calls = g0$calls,
    grad_opt_s = g1$seconds, control_ratio = ctl,
    fit_s = t_fit, mem_mb = mem,
    logLik = as.numeric(stats::logLik(fit)),
    key = truth_key, est = unname(b[truth_key]), truth = truth_value,
    est_lo = i_k[1L], est_hi = i_k[2L],
    sd_all = sd_all, sd_true = sd_truth, n_sd = length(sds),
    cor_max_abs = cor_max, cor_true = cor_true),
    ex, list(diag = scale_diag(fit))))
  list(fit = fit, ci = ci, interval = i_k, extra = ex)
}

# The rlddm row's own quantities: what the per-subject bound is worth,
# measured against the truths the draw used and against each learner's
# own fastest response.
learn_rlddm_extra <- function(fit, d) {
  one <- d[match(levels(d$id), as.character(d$id)), , drop = FALSE]
  hat <- as.numeric(suppressWarnings(
    frmtmb.eam::ndt_time(fit, newdata = one)))
  own <- as.numeric(tapply(d$rt, d$id, min))[match(as.character(one$id),
                                                   levels(d$id))]
  tru <- as.numeric(attr(d, "ndt_subject"))[match(as.character(one$id),
                                                  levels(d$id))]
  bd <- frmtmb::single_response(fit)[["family"]][["ndt_bound"]]
  fl <- mean(bd[["floors"]])
  # which pair carries cor_max_abs, by name, so that a reader does not
  # have to reconstruct it from sd_all
  cr <- stats::cov2cor(VarCorr(fit)[[1L]])
  lo <- which(abs(cr) == max(abs(cr[lower.tri(cr)])), arr.ind = TRUE)
  cor_pair <- paste(rownames(cr)[lo[1L, 1L]], "vs",
                    colnames(cr)[lo[1L, 2L]])
  nd <- suppressWarnings(stats::predict(
    fit, newdata = one[1L, , drop = FALSE], dpar = "ndt",
    type = "response", re.form = NA, se.fit = TRUE))
  pop <- as.numeric(nd$fit[1L]) * fl
  pop_se <- as.numeric(nd$se.fit[1L]) * fl
  list(ndt_pop = pop, ndt_pop_se = pop_se,
       ndt_pop_true = learn_truth$ndt,
       # RECORDED, not asserted, and dev/rlddm-findings.md says why:
       # `plogis(b0) * mean(floors)` is not the population non-decision
       # time, and its standard error covers only the fraction. The gap
       # to the stated 0.25 decomposes into a lognormal mean-against-
       # median offset, the realized draw's own deviation, the fit's
       # error, and this conversion, and only the third is the fit.
       ndt_pop_z = scale_z(pop, pop_se, learn_truth$ndt),
       ndt_floor_mean = fl,
       ndt_frac = as.numeric(nd$fit[1L]),
       ndt_sub_mean = mean(hat), ndt_sub_mean_true = mean(tru),
       # the two statistics that DO separate the model from its floors,
       # each as a ratio to a spread this run measured
       ndt_mean_err_ratio = abs(mean(hat) - mean(tru)) / stats::sd(tru),
       ndt_rmse_ratio = sqrt(mean((hat - tru)^2)) / stats::sd(tru),
       ndt_rmse_ms = 1000 * sqrt(mean((hat - tru)^2)),
       ndt_cor_true = stats::cor(hat, tru),
       sd_ndt_natural = stats::sd(hat),
       sd_ndt_natural_true = stats::sd(tru),
       # A COMPARISON BESIDE IT, because this field without one is what
       # produced 1.0a's B5 and produced it again here. `ndt_i` is a
       # fraction times a floor, so the floors alone give a spread even
       # with the random effect exactly zero: this is the spread a
       # COMMON fraction of each learner's own floor would give, which
       # needs no second fit. Measured separately by the review's
       # `dev/rev-rlddm-floors.R`, a fit with the random effect on `ndt`
       # switched OFF returns 0.0380049 against a truth of 0.0393041
       # where the full model returns 0.0333713, so on this statistic
       # the random effect is WORSE than leaving it out. Read the
       # per-learner error and the correlation instead.
       sd_ndt_floor_only = stats::sd(as.numeric(nd$fit[1L]) *
                                       attr(d, "own_floor")),
       cor_max_pair = cor_pair,
       ndt_margin_min_ms = 1000 * min(own - hat),
       ndt_below_own_floor = sum(own - hat > 0),
       n_subjects = length(hat),
       n_bad_se = sum(!is.finite(fit$sdr$sd)))
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
  # ITEM 1.0b: ndt_group(id) is what makes `ndt`'s bound each LEARNER's
  # own fastest response rather than the fastest of all one hundred.
  form <- bf(rt | dec(choice) + reward(pay1, pay2) +
               ndt_group(id) ~ 1 + (1 | p | id),
             drift ~ 1 + (1 | p | id), bs ~ 1 + (1 | p | id),
             ndt ~ 1 + (1 | p | id), bias = 0.5)
  r <- learn_scale_run("learn-rlddm", form,
                       rlddm(subject = id, trial = trial), d,
                       "alpha.(Intercept)",
                       stats::qlogis(learn_truth$alpha), "(Intercept)",
                       "alpha=0.5;drift=1;bs=0.2;ndt=no link truth",
                       extra = learn_rlddm_extra,
                       # NOT ZERO, and the row said zero until punch
                       # round 2. The design draws four independent
                       # deviations, but the block this row fits is not
                       # the block it drew: see learn_rlddm_data(). The
                       # truth is computed per replicate from the truths
                       # and the floors, with no fit, so it moves with
                       # the seed instead of being pinned to one.
                       cor_true = learn_cor_max(
                         attr(d, "dev_fitted_par")))
  expect_true(is.finite(as.numeric(stats::logLik(r$fit))))
  if (scale_small()) {
    # at five learners by thirty trials there is not enough data for
    # any of this, and the row is a code path rather than a measurement
    expect_true(r$extra$ndt_below_own_floor == r$extra$n_subjects)
  } else {
    # THE PLAN'S ACCEPTANCE CRITERION for item 1.0b, and it is what the
    # row asserted the negation of through 0.3.0: conv=1, a maximum
    # gradient of 6.36e+09, pdHess FALSE and four NaN standard errors.
    expect_identical(r$fit$opt$convergence, 0L)
    expect_true(isTRUE(frmtmb::diagnose(r$fit, quiet = TRUE)$pdHess))
    expect_true(all(is.finite(r$fit$sdr$sd)))
    expect_identical(r$extra$n_bad_se, 0L)
    # A converged fit that recovers nothing is not the point, so what
    # follows is recovery. The population non-decision time is RECORDED
    # rather than asserted: `plogis(b0) * mean(floors)` is not that
    # quantity and its standard error covers only the fraction, so a z
    # against the design's stated 0.25 tests the instrument as much as
    # the fit. dev/rlddm-findings.md decomposes it. What IS asserted is
    # the per-learner recovery the plan's item 2.1 asks for, each stated
    # as a ratio to a spread this run measured.
    #
    # No learner past its own fastest response. This one CANNOT FAIL
    # and is kept for what it does check. Under a per-group bound
    # ndt_i = plogis(eta_i) * floor_i and plogis < 1 for every finite
    # eta, so ndt_i < floor_i identically; what it tests is that the
    # floor LOOKUP paired each learner with its own floor and not with
    # somebody else's, which is the defect item 1.0a's review found by
    # droplevels(). It is not a check on the fit and the plan's
    # acceptance criterion should not read it as one. The margin below
    # is the informative number.
    expect_identical(r$extra$ndt_below_own_floor, r$extra$n_subjects)
    # and one that CAN fail: the fit does not sit at the wall. A bound
    # that was binding would put the tightest learner within a hair of
    # its own floor, which is where the global-bound arm sat. Stated
    # against the spread of the drawn truths rather than as a constant.
    expect_gt(r$extra$ndt_margin_min_ms /
                (1000 * r$extra$sd_ndt_natural_true), 0.1)
    # the per-learner error is smaller than the spread of the truths it
    # is estimating, that is, the fit beats reporting one number for
    # everybody
    expect_lt(r$extra$ndt_rmse_ratio, 1)
    # and the population mean is inside that same spread
    expect_lt(r$extra$ndt_mean_err_ratio, 1)
    # the fitted per-learner non-decision times track the drawn ones.
    # Stated as a correlation rather than as a tolerance on any
    # difference, because a correlation has a scale of its own.
    expect_gt(r$extra$ndt_cor_true, 0.5)
    # the learning rate's interval covers its truth
    tr <- stats::qlogis(learn_truth$alpha)
    expect_true(all(is.finite(r$interval)))
    expect_true(r$interval[1L] <= tr && r$interval[2L] >= tr)
  }
})
