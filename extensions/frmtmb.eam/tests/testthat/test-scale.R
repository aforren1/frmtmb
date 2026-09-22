## Phase 0 of dev/extension-gaps-plan.md: the eam row.
##
## 30 subjects x 400 trials, two conditions, random effects on the drift
## rate, the boundary separation and the non-decision time, and a second
## run with wiener(variability = "sv"). 12,000 rows.
##
## What the row decides: whether a hierarchical DDM is minutes or hours,
## and whether the bounded non-decision-time link tolerates a random
## effect at all. The second question is the one item 2.1 of the plan
## expects trouble on, because ndt's link is bounded by the GLOBAL
## fastest response, so a subject deviation on it is a deviation on a
## fraction of one subject's floor. The tier RECORDS that recovery and
## does not assert it: Phase 0 measures cost, Phase 2 validates
## recovery, and a tier that fails on the question it was built to ask
## would report nothing.
##
## See dev/scale-findings.md for the numbers this produced.

# The truth, in one place, because both designs draw from it.
eam_truth <- list(mu0 = 0.4, mu_cond = 0.9, bs = 1.4, ndt = 0.25,
                  sd_mu = 0.35, sd_lbs = 0.20, sd_lndt = 0.12,
                  sv = 0.4)

# 30 subjects x 400 trials. The subject deviations are drawn on the
# scale the family estimates each parameter on, except the non-decision
# time: its link is bounded by a quantity the data have not produced
# yet, so its deviation is drawn on the natural scale and the truth is
# stated there too.
eam_scale_data <- function(seed = 20260908L, sv = 0,
                           ns = if (scale_small()) 4L else 30L,
                           nt = if (scale_small()) 40L else 400L) {
  tr <- eam_truth
  set.seed(seed)
  u_mu <- stats::rnorm(ns, 0, tr$sd_mu)
  u_bs <- stats::rnorm(ns, 0, tr$sd_lbs)
  u_nd <- stats::rnorm(ns, 0, tr$sd_lndt)
  s <- rep(seq_len(ns), each = nt)
  cond <- rep(rep(0:1, each = nt / 2L), times = ns)
  d <- ddm_simulate(ns * nt,
                    mu = tr$mu0 + tr$mu_cond * cond + u_mu[s],
                    bs = tr$bs * exp(u_bs[s]),
                    ndt = tr$ndt * exp(u_nd[s]),
                    bias = 0.5, sv = sv)
  d$s <- factor(s)
  d$cond <- factor(cond, labels = c("a", "b"))
  # the per-subject non-decision times the draw used, so that the fit's
  # own spread can be compared with them on the NATURAL scale. The
  # deviation itself is drawn on the LOG scale, just above; it is the
  # recorded TRUTH that is natural, and that is the part that matters,
  # because the link ndt is estimated on is a scaled logit whose bound
  # is the global fastest response and a standard deviation on that
  # scale has no truth stated in advance.
  attr(d, "ndt_subject") <- tr$ndt * exp(u_nd)
  d
}

# ITEM 1.0a. `ndt_group(s)` is what makes the non-decision time's bound
# each SUBJECT's own fastest response rather than the whole data set's.
# The bound-lifted arm cannot carry it, because max_ndt and ndt_group()
# set the same bound to different things and the pair is refused, so
# that arm asks for the form without it.
eam_scale_form <- function(group = TRUE) {
  if (group) {
    bf(rt | dec(upper) + ndt_group(s) ~ cond + (1 | s), bs ~ 1 + (1 | s),
       ndt ~ 1 + (1 | s), bias = 0.5)
  } else {
    bf(rt | dec(upper) ~ cond + (1 | s), bs ~ 1 + (1 | s),
       ndt ~ 1 + (1 | s), bias = 0.5)
  }
}

# One row of the row's own data, for a population-level prediction.
d1_of <- function(r) {
  eam_scale_data(sv = 0)[1L, , drop = FALSE]
}

# The Wald interval a row reports, by the name confint() gives it.
eam_ci <- function(ci, key) {
  j <- grep(key, rownames(ci), fixed = TRUE)
  if (!length(j)) return(c(NA_real_, NA_real_))
  as.numeric(ci[j[1L], 1:2])
}

# One design, measured. Everything the plan's Phase 0 asks for: the
# tape build, one gradient, the whole frm() call including sdreport(),
# the peak R heap and the estimates against the truth.
eam_scale_run <- function(row, fam, sv, group = TRUE) {
  d <- eam_scale_data(sv = sv)
  form <- eam_scale_form(group)
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
  b <- unlist(fixef_by_dpar(fit))
  # ONLY the grouped link reports a fraction. Under `ndt_group()` the
  # link is a plain logit on a fraction of the row's own bound, so the
  # population non-decision time is that fraction on the mean of the
  # bounds the fit used, and the delta-method standard error rides with
  # it. Under the scalar bound the bound stays INSIDE the link, a scaled
  # logit onto (0, ub), so `frm_linpred()` already reports a time; this line
  # used to multiply by `ub` in that case too and the `eam-unbounded`
  # row recorded `time * ub`, 0.1315 s where `ndt_time()` on the same
  # fit said 0.2921. The expectation below is what now catches it.
  nd <- suppressWarnings(
    frm_linpred(fit, newdata = d[1L, , drop = FALSE], dpar = "ndt",
                   type = "response", re_formula = NA, se.fit = TRUE))
  ndt_frac <- as.numeric(nd$fit[1L])
  ndt_frac_se <- as.numeric(nd$se.fit[1L])
  rsp <- frmtmb::single_response(fit)
  bnd <- rsp[["family"]][["ndt_bound"]]
  to_time <- if (is.null(bnd[["floors"]])) 1 else mean(bnd[["floors"]])
  ndt_hat <- ndt_frac * to_time
  ndt_se <- ndt_frac_se * to_time
  vc <- varcorr_matrices(fit)
  sds <- vapply(vc, function(m) sqrt(m[1L, 1L]), numeric(1))
  # by name as well as by position, so that a reader can check the
  # position the three named fields below assume
  sd_all <- paste(names(vc), formatC(sds, digits = 4, format = "g"),
                  sep = "=", collapse = ";")
  # the fitted per-subject non-decision time, one row per subject, on
  # the natural scale where the simulator's spread is stated
  one <- d[match(levels(d$s), as.character(d$s)), , drop = FALSE]
  ndt_sub <- suppressWarnings(as.numeric(ndt_time(fit, newdata = one)))
  # how much room each subject's fitted non-decision time has left
  # against its OWN fastest response, which is the quantity the
  # per-subject bound exists to respect
  own <- as.numeric(tapply(d$rt, d$s, min))[match(as.character(one$s),
                                                  levels(d$s))]
  margin_ms <- 1000 * (own - ndt_sub)
  # THE SAME QUANTITY BY TWO ROUTES, which checks this file's own
  # arithmetic rather than the model: the population prediction put back
  # on the response's scale, and `ndt_time()` on the subjects. The two
  # differ by the random effects alone, so the tolerance is built from
  # the spread this run measured and the standard error it reported, and
  # never from a constant. The ungrouped row used to fail it by a factor
  # of its own bound: `ndt` was recorded as 0.1315 s where `ndt_time()`
  # on the same fit said 0.2921 (dev/eamhier-findings.md).
  testthat::expect_lt(abs(ndt_hat - mean(ndt_sub)),
                      3 * (stats::sd(ndt_sub) + ndt_se))

  tr <- eam_truth
  i_cond <- eam_ci(ci, "condb")
  scale_record(
    row, rows = nrow(d), n_par = length(fit$opt$par),
    build_s = bd$build_s, frame_s = bd$frame_s,
    grad_start_s = g0$seconds, grad_start_calls = g0$calls,
    grad_opt_s = g1$seconds, control_ratio = ctl,
    fit_s = t_fit, mem_mb = mem,
    logLik = as.numeric(stats::logLik(fit)),
    mu0 = unname(b["mu.(Intercept)"]), mu0_true = tr$mu0,
    mu_cond = unname(b["mu.condb"]), mu_cond_true = tr$mu_cond,
    mu_cond_lo = i_cond[1L], mu_cond_hi = i_cond[2L],
    bs = exp(unname(b["bs.(Intercept)"])), bs_true = tr$bs,
    ndt = ndt_hat, ndt_se = ndt_se, ndt_true = tr$ndt,
    # BOTH z scores, because the row is read for two different
    # questions and the assertions above use the second. `ndt_z` scores
    # against the population constant, which is the MEDIAN of the
    # per-subject draw and is what a reader comparing rows across seeds
    # wants; `ndt_z_drawn` scores against the mean of the 30 subjects
    # this seed drew, which is the estimand the fit is actually
    # estimating. They differ by more than the fit's own standard
    # error: the first exceeds 4 on 8 of 60 replicates and the second
    # on 0 of 60 (dev/eamhier-findings.md).
    ndt_z = scale_z(ndt_hat, ndt_se, tr$ndt),
    ndt_z_drawn = scale_z(ndt_hat, ndt_se,
                          mean(attr(d, "ndt_subject"))),
    ndt_frac = ndt_frac, ndt_frac_se = ndt_frac_se,
    # 1 when `frm_linpred()` already reported a time, the mean floor when it
    # reported a fraction, so the row says which scale it was read on
    ndt_to_time = to_time,
    ndt_sub_mean = mean(ndt_sub),
    ndt_sub_mean_true = mean(attr(d, "ndt_subject")),
    ndt_margin_min_ms = min(margin_ms),
    ndt_below_own_floor = sum(margin_ms > 0),
    n_subjects = length(ndt_sub),
    sd_mu = unname(sds[1L]), sd_mu_true = tr$sd_mu,
    sd_bs = unname(sds[2L]), sd_bs_true = tr$sd_lbs,
    sd_ndt_link = unname(sds[3L]),
    sd_ndt_natural = stats::sd(ndt_sub),
    sd_ndt_natural_true = stats::sd(attr(d, "ndt_subject")),
    sd_all = sd_all,
    diag = scale_diag(fit))

  list(fit = fit, ci = ci, b = b, t_fit = t_fit)
}

test_that("the eam scale row fits and reports its cost", {
  skip_unless_scale()
  scale_row_on("eam")
  r <- eam_scale_run("eam", wiener(), sv = 0)
  expect_true(is.finite(as.numeric(stats::logLik(r$fit))))
  i <- eam_ci(r$ci, "condb")
  if (scale_small()) {
    # at four subjects the ndt component collapses instead of running
    # away, so the fit converges and the ordinary assertion applies
    expect_true(i[1L] <= eam_truth$mu_cond &&
                  i[2L] >= eam_truth$mu_cond)
  } else {
    # ITEM 1.0a of dev/extension-gaps-plan.md, LANDED. Through 0.6.0 this
    # fit could not converge: `ndt`'s scaled logit was bounded by the
    # GLOBAL fastest response and 20 of the 30 subjects have a true
    # `ndt` above it, so every standard error was NaN and this row
    # asserted the defect. `ndt_group(s)` now bounds each subject by its
    # own fastest response, and the assertions are the recovery ones.
    expect_identical(r$fit$opt$convergence, 0L)
    expect_true(all(is.finite(i)))
    expect_true(i[1L] <= eam_truth$mu_cond &&
                  i[2L] >= eam_truth$mu_cond)
    # The plan's target for this row is a log-likelihood of at least
    # -7027.4, the value the bound-lifted arm reached, against -7148.8
    # with the global bound; the run RECORDS it in the `logLik` field
    # above and `dev/ndt-findings.md` carries the comparison. It is not
    # asserted here, and deliberately: helper-scale.R's own rule is
    # that nothing absolute goes in a tier assertion, and a
    # log-likelihood carried over from another arm at another seed is
    # the most brittle constant of all. What IS asserted is what this
    # run measures about itself.
    expect_true(all(is.finite(r$fit$sdr$sd)))
    # The population non-decision time, as a z against the fit's OWN
    # standard error, which is what scale_z() is here for, and against
    # the DRAW's own mean rather than against `eam_truth$ndt`.
    #
    # WHY the target moved. `ndt` is drawn as 0.25 * exp(u) with u
    # normal of standard deviation 0.12, so 0.25 is the MEDIAN of the
    # per-subject non-decision times and their mean is 0.2518. With 30
    # subjects the sample mean of that log-normal carries a standard
    # error of about 5.5 ms, more than twice the 2.5 ms this fit
    # reports, so scoring against the constant tests the draw and not
    # the fit. Measured over 60 replicates of this design
    # (dev/eamhier-findings.md): against 0.25 the assertion passes on
    # 52 of 60 and reaches z = 10.10 at seed 20260936, whose 30
    # subjects have a mean ndt of 0.27013 and whose fit returns
    # 0.27545; against the draw's own mean, 60 of 60 with a maximum of
    # 2.60. The eam-sv row below already scores against the draw, and
    # is 60 of 60 there at a maximum of 2.98.
    d <- eam_scale_data(sv = 0)
    one <- d[match(levels(d$s), as.character(d$s)), , drop = FALSE]
    nd <- suppressWarnings(frm_linpred(
      r$fit, newdata = d1_of(r), dpar = "ndt", type = "response",
      re_formula = NA, se.fit = TRUE))
    bd <- frmtmb::single_response(r$fit)[["family"]][["ndt_bound"]]
    fl <- mean(bd[["floors"]])
    expect_lt(scale_z(as.numeric(nd$fit[1L]) * fl,
                      as.numeric(nd$se.fit[1L]) * fl,
                      mean(attr(d, "ndt_subject"))), 4)
    # ITEM 1.0a's note, and item 2.1's first carry-in: what separates
    # this model from its floors is the PER-SUBJECT error, so that is
    # what is asserted, against the between-subject spread this same
    # run measured rather than against a constant. Over the same 60
    # replicates the ratio is under 1 on 60 of 60 and 0.468 at worst.
    hat <- as.numeric(ndt_time(r$fit, newdata = one))
    t0 <- attr(d, "ndt_subject")
    expect_lt(sqrt(mean((hat - t0)^2)), stats::sd(t0))
    # and the constraint the per-subject bound exists to respect: every
    # subject's fitted non-decision time below its OWN fastest response
    own <- as.numeric(tapply(d$rt, d$s, min))
    expect_true(all(hat < own))
  }
})

test_that("the eam scale row fits with the ndt bound lifted", {
  skip_unless_scale()
  scale_row_on("eam-unbounded")
  # The SAME data as the `eam` row, fitted with the non-decision-time
  # bound raised above every subject's true value.
  #
  # This arm exists because "the fit did not converge" is not a cause.
  # wiener() bounds ndt by the GLOBAL fastest response, so a subject
  # whose true ndt is above that bound cannot be represented at any
  # value of the random effect. Raising `max_ndt` past the largest true
  # ndt puts the truth back inside the parameter space, at the price of
  # rows the density cannot reach, which is what `allow_unreachable`
  # opts into. If this arm converges and the `eam` arm does not, the
  # bound is the cause and not the optimizer.
  fam <- wiener(max_ndt = 0.45, allow_unreachable = TRUE)
  r <- eam_scale_run("eam-unbounded", fam, sv = 0, group = FALSE)
  expect_true(is.finite(as.numeric(stats::logLik(r$fit))))
})

test_that("the eam scale row fits with across-trial drift variability", {
  skip_unless_scale()
  scale_row_on("eam-sv")
  r <- eam_scale_run("eam-sv", wiener(variability = "sv"),
                     sv = eam_truth$sv)
  expect_true(is.finite(as.numeric(stats::logLik(r$fit))))
  i <- eam_ci(r$ci, "condb")
  if (scale_small()) {
    expect_true(i[1L] <= eam_truth$mu_cond && i[2L] >= eam_truth$mu_cond)
  } else {
    # WHAT THIS ROW NO LONGER ASSERTS, and the measurement that says
    # why. It used to assert that the condition effect's Wald interval
    # covers 0.9. Under the per-subject bound it does not, at this one
    # seed, and the bound is NOT what broke it. The same 12,000 rows,
    # both parameterizations, `dev/ndt-scripts/ndt-sv-arm.R`:
    #
    #   per-group bound: mu.condb 0.8223 (0.7621, 0.8825), sv 0.3105,
    #     population ndt within 3.8 ms of the truth, logLik -6926.65
    #   global bound:    mu.condb 0.8565 (0.7936, 0.9195), sv 0.5589,
    #     population ndt 22.7 ms BELOW the truth, logLik -7129.88
    #
    # BOTH underestimate the condition effect on this draw, by 0.078
    # and 0.044 against a truth of 0.9, and the global arm covers only
    # because its interval is 8 percent wider. Neither recovers `sv`,
    # which is 0.31 and 0.56 against 0.4. What moved is the
    # drift-against-variability trade-off, not the bound: with `ndt`
    # pinned 22.7 ms low the drift spread is absorbed by `sv`, and with
    # `ndt` right it is not. That is a question for Phase 2's 60
    # replicates, which is where recovery is decided, so this row
    # RECORDS the estimate and asserts only what one draw can carry.
    #
    # Item 2.1 has since measured the trade-off over replicates, in a
    # design with NO hierarchy and NO per-group bound: 60 single-subject
    # replicates of the same 12,000 trials recover `sv` (57 of 59
    # covered once one collapsed fit is set aside), the estimated `sv`
    # and the estimated condition effect correlate at 0.59 (Spearman),
    # and every replicate that misses on the condition effect has a low
    # `sv`. So a single seed's 0.31 was inside this estimator's
    # ordinary spread. And THIS row's arm, at 60 replicates, covers the
    # condition effect 57 times and `sv` 57 times, both nominal.
    # dev/eamhier-findings.md.
    expect_identical(r$fit$opt$convergence, 0L)
    expect_true(all(is.finite(i)))
    # the non-decision time itself does recover under the bound, which
    # is what this lane changed
    d <- eam_scale_data(sv = eam_truth$sv)
    one <- d[match(levels(d$s), as.character(d$s)), , drop = FALSE]
    own <- as.numeric(tapply(d$rt, d$s, min))
    hat <- as.numeric(ndt_time(r$fit, newdata = one))
    expect_true(all(hat < own))
    # the per-subject error against the between-subject spread the same
    # run measured, as in the `eam` row above: under 1 on 60 of 60
    # replicates of this arm, 0.483 at worst
    t0 <- attr(d, "ndt_subject")
    expect_lt(sqrt(mean((hat - t0)^2)), stats::sd(t0))
    # a z against the fit's own standard error rather than an absolute
    # 20 ms, which is the rule helper-scale.R states above scale_z()
    nd <- suppressWarnings(frm_linpred(
      r$fit, newdata = one[1L, , drop = FALSE], dpar = "ndt",
      type = "response", re_formula = NA, se.fit = TRUE))
    bd <- frmtmb::single_response(r$fit)[["family"]][["ndt_bound"]]
    fl <- mean(bd[["floors"]])
    expect_lt(scale_z(as.numeric(nd$fit[1L]) * fl,
                      as.numeric(nd$se.fit[1L]) * fl,
                      mean(attr(d, "ndt_subject"))), 4)
  }
})
