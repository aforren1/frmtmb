## Phase 0 of dev/extension-gaps-plan.md: the spline curves row.
##
## `s(t) + s(t, subject, bs = "fs")` at 40 subjects, 200 points each,
## then frm_curve() and frm_curve_feature() on it. 8,000 rows.
##
## What the row decides: the memoized joint-precision solve was 6.9 s at
## 8006 random coefficients in the docs, so what is unmeasured is the
## FEATURE SEARCH, which calls predict() once per Newton step. Both the
## curve and the feature are timed separately from the fit here.
##
## The plan's spline row names only the curves design, so the Royston
## and Parmar survival design in the same package's "Realistic scale"
## line is Phase 2's item 2.5 and is not run here.
##
## The subject curves are a peak whose HEIGHT and LOCATION vary by
## subject, which is the movement-science design the package's own
## fs test uses, at 40 subjects and 200 points instead of 20 and 30.
##
## See dev/scale-findings.md for the numbers this produced.

spline_truth <- list(n_sub = 40L, n_t = 200L, height = 1, loc = 0.45,
                     width = 0.16, sd_h = 0.12, sd_loc = 0.04,
                     sigma = 0.06, k_pop = 20L, k_fs = 6L)

spline_peak <- function(t, h, s, w) h * exp(-0.5 * ((t - s) / w)^2)

spline_scale_data <- function(seed = 20260908L) {
  tr <- spline_truth
  set.seed(seed)
  n_sub <- if (scale_small()) 6L else tr$n_sub
  n_t <- if (scale_small()) 30L else tr$n_t
  sub <- rep(seq_len(n_sub), each = n_t)
  d <- data.frame(subject = factor(sub),
                  t = rep(seq(0, 1, length.out = n_t), times = n_sub))
  h <- stats::rnorm(n_sub, tr$height, tr$sd_h)
  s <- stats::rnorm(n_sub, tr$loc, tr$sd_loc)
  d$v <- spline_peak(d$t, h[sub], s[sub], tr$width) +
    stats::rnorm(nrow(d), 0, tr$sigma)
  d
}

test_that("the spline curves scale row fits and reports its cost", {
  skip_unless_scale()
  scale_row_on("spline")
  tr <- spline_truth
  d <- spline_scale_data()
  form <- frmtmb::bf(v ~ s(t, k = tr$k_pop) +
                       s(t, subject, bs = "fs", k = tr$k_fs))
  scale_mem_reset()

  bd <- scale_build(form, family = stats::gaussian(), data = d)
  g0 <- scale_grad(bd$dry$obj, bd$dry$obj$par)
  ctl <- scale_control(bd$dry$obj, bd$dry$obj$par, g0$calls)
  bd$dry <- NULL

  fit <- NULL
  # the fs fit reaches a maximum absolute gradient of about 1.5e-3 on
  # some platforms, which the package's own fs test already records;
  # scale_diag() reports whatever this run reached
  t_fit <- scale_elapsed(
    fit <- suppressWarnings(frmtmb::frm(form, family = stats::gaussian(),
                                        data = d, se = TRUE)))
  g1 <- scale_grad(fit$obj, fit$opt$par)

  g <- data.frame(t = seq(0.05, 0.95, length.out = 80))
  # THE COLD CALL FIRST, and on purpose.
  #
  # frm_curve() memoizes the joint-precision solve on the fit object, so
  # the solve is paid ONCE per fit by whichever curve call happens
  # first, and every later curve call on that fit is cheap. That means
  # `scale_interleave()`'s per-arm "first" is only cold for arm 1: by
  # the time arm 3 runs, arms 1 and 2 have warmed the cache. Measuring
  # the feature search here, before anything else touches the curve
  # machinery, is what makes its cold number a real one.
  t_solve <- scale_elapsed(
    frm_curve_feature(fit, var = "t", type = "maximum", newdata = g))
  # interleaved and replicated from here on, all arms warm
  pf <- scale_interleave(list(
    curve_sim = function() {
      frm_curve(fit, newdata = g, re.form = NA, simultaneous = TRUE,
                nsim = 10000L, seed = 1)
    },
    curve_pw = function() {
      frm_curve(fit, newdata = g, re.form = NA, simultaneous = FALSE)
    },
    feature = function() {
      frm_curve_feature(fit, var = "t", type = "maximum", newdata = g)
    }))
  cv <- frm_curve(fit, newdata = g, re.form = NA, simultaneous = TRUE,
                  nsim = 10000L, seed = 1)
  ft <- frm_curve_feature(fit, var = "t", type = "maximum", newdata = g)
  mem <- scale_mem_peak_mb()

  scale_record(
    "spline", rows = nrow(d), n_coef = length(fit$estimates$b),
    n_par = length(fit$opt$par),
    build_s = bd$build_s, frame_s = bd$frame_s,
    grad_start_s = g0$seconds, grad_start_calls = g0$calls,
    grad_opt_s = g1$seconds, control_ratio = ctl,
    fit_s = t_fit, mem_mb = mem,
    rounds = pf$rounds,
    # the memoized solve, measured once on a cold fit through the
    # cheapest call that pays it
    solve_cold_s = t_solve,
    # everything below is WARM: the solve is already cached
    curve_sim_s = pf$seconds[["curve_sim"]],
    curve_sim_spread = pf$spread[["curve_sim"]],
    curve_pw_s = pf$seconds[["curve_pw"]],
    feature_s = pf$seconds[["feature"]],
    feature_spread = pf$spread[["feature"]],
    curve_n_predict = attr(cv, "check")$n_predict,
    curve_cov_rel_error = attr(cv, "check")$cov_rel_error,
    logLik = as.numeric(stats::logLik(fit)),
    peak_t = ft$.estimate, peak_t_se = ft$.se, peak_t_true = tr$loc,
    peak_t_lo = ft$.lower_ci, peak_t_hi = ft$.upper_ci,
    peak_v = ft$.value, peak_v_true = tr$height,
    sigma = stats::sigma(fit), sigma_true = tr$sigma,
    diag = scale_diag(fit))

  expect_true(is.finite(as.numeric(stats::logLik(fit))))
  # the location of the peak is what the feature search exists to find,
  # and its own interval is the tolerance
  expect_true(ft$.lower_ci < tr$loc && ft$.upper_ci > tr$loc)
})
