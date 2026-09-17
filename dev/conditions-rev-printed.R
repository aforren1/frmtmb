# Reviewer, lane wt-conditions (round 2): printed output of real
# converted sites, base build against lane build, on the same triggering
# input. For each arm and mode this writes ONE script that runs every
# trigger as its own top-level expression in a fresh Rscript, so
# deferred warnings print where R prints them, and captures stdout and
# stderr together (2>&1) and stderr alone.
#
# Modes:
#   warn0     try(EXPR) at top level, options(warn = 0): error text from
#             try(), deferred "Warning message:" blocks, messages
#   warn1     the same with options(warn = 1): immediate warnings
#   suppress  suppressWarnings(suppressMessages(try(EXPR)))
#   catch     print(tryCatch(EXPR, error =, warning =, message =)) and a
#             withCallingHandlers() that muffles and prints what it saw
#   record    in-process: class, message and call of the first condition
#
#   Rscript dev/conditions-rev-printed.R
root <- "C:/Users/adf44/source/r/frmtmb-wt-conditions"
out <- file.path(root, "dev/conditions-rev-log/printed")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
rscript <- "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
libs <- c(base = "C:/Users/adf44/source/r/rellib-r3",
          lane = "C:/Users/adf44/source/r/conditions-lib")

setup <- '
suppressPackageStartupMessages({
  library(frmtmb); library(frmtmb.eam); library(frmtmb.ode)
  library(frmtmb.spline); library(frmtmb.latent); library(frmtmb.learn)
  library(frmtmb.coupling); library(frmtmb.sample)
})
options(frmtmb.notices = TRUE, width = 80)
set.seed(20260917)
d <- data.frame(y = rnorm(60), x = rnorm(60), z = rnorm(60),
                g = factor(rep(1:10, each = 6)),
                one = factor(rep("a", 60)), obs = factor(1:60))
d$n <- rep(10L, 60); d$k <- rbinom(60, 10, 0.4)
d$yb <- rbinom(60, 1, 0.5)
d$day <- as.Date("2020-01-01") + 0:59
d$w <- ifelse(seq_len(60) %% 7 == 0, NA, d$x)
fg <- frm(y ~ x, data = d, family = gaussian())
fb <- frm(k | trials(n) ~ x, data = d, family = binomial())
fm <- frm(y ~ x + (1 | g), data = d, family = gaussian())
fsig <- frm(bf(y ~ x, sigma ~ z), data = d, family = gaussian())
dat_rt <- data.frame(rt = 0.5 + rexp(40), upper = rbinom(40, 1, 0.5))
p64x <- rnorm(64); p64y <- rnorm(64)
'

triggers <- c(
  # core errors, frm() pipeline, mostly paste-built and multi-argument
  frm_family_unknown = 'frm(y ~ x, data = d, family = "not_a_family")',
  frm_missing_var = 'frm(y ~ nope, data = d, family = gaussian())',
  frm_start_unknown = 'frm(y ~ x, data = d, start = list(bogus = 1))',
  frm_start_len = 'frm(y ~ x, data = d, start = list(beta = 1:5))',
  frm_trials_missing = 'frm(k ~ x, data = d, family = binomial())',
  frm_empty_data = 'frm(y ~ x, data = d[0, ], family = gaussian())',
  frm_reml_nl = 'frm(bf(y ~ a * x, a ~ 1, nl = TRUE), data = d, REML = TRUE)',
  frm_optimizer_fn = 'frm(y ~ x, data = d, control = frmtmb_control(optimizer = function(...) stop("inner optimizer broke")))',
  frm_optimizer_name = 'frm(y ~ x, data = d, control = frmtmb_control(optimizer = "nosuch"))',
  bf_bad_dpar = 'bf(y ~ x, nosuchdpar ~ 1)',
  set_prior_class = 'frm(y ~ x, data = d, prior = set_prior("normal(0, 1)", class = "nosuch"))',
  set_prior_bad_dist = 'set_prior("nosuchdist(1, 2)")',
  predict_newdata_missing = 'predict(fg, newdata = data.frame(z = 1))',
  predict_bad_arg = 'predict(fg, re_frmula = NA)',
  fitted_typo = 'fitted(fg, re_frmula = NA)',
  coef_typo = 'coef(fg, robst = TRUE)',
  nobs_extra = 'nobs(fg, 7)',
  print_digits = 'print(fg, digits = 3)',
  anova_mix = 'anova(fg, frm(y ~ x, data = d, REML = TRUE))',
  anova_space = 'anova(fg, frm(y ~ z, data = d))',
  confint_parm = 'confint(fg, parm = "nosuchparm")',
  hypothesis_bad = 'hypothesis(fg, "nosuchpar = 0")',
  ranef_fixed = 'ranef(fg)',
  vcov_cluster_len = 'vcov_cluster(fg, cluster = 1:3)',
  simulate_nsim = 'simulate(fg, nsim = -1)',
  cond_eff_bad = 'conditional_effects(fg, effects = "nosuch")',
  bf_bad_dpar2 = 'frm(bf(y ~ x, nosuch ~ 1), data = d)',
  confint_level = 'confint(fg, level = 2)',
  frm_data_chr = 'frm(y ~ x, data = "notdf")',
  prior_coef = 'frm(y ~ x, data = d, prior = set_prior("normal(0, 1)", class = "b", coef = "nosuch"))',
  predict_type_matcharg = 'predict(fg, type = "bogus")',
  control_matcharg = 'frmtmb_control(check_nlev_1 = "bogus")',
  set_prior_stopifnot = 'set_prior("normal(0)")',
  spline_deriv_newdata = 'frm_curve_deriv(fg, var = "x")',
  coupling_sfreq = 'frm_cross_spectrum(p64x, p64y, sfreq = -1)',
  ode_lincmt_both = 'frm_lincmt(parms = list(ka = 1, CL = 2, ke = 0.2, V = 10), times = 1:5, ncmt = 1, depot = TRUE)',
  coupling_frange = 'frm_cross_spectrum(p64x, p64y, frange = c(1, 0))',
  spline_curve_empty = 'frm_curve(fg, newdata = data.frame())',
  getME_mv_Zt = 'lme4::getME(frm(mvbf(bf(y ~ x + (1 | g)), bf(z ~ x + (1 | g))), data = d), "Zt")',
  frm_periodogram_short = 'frm_periodogram(1)',
  # extensions
  eam_wiener_linkinv = 'wiener()$links[["ndt"]]$linkinv(0)',
  eam_wiener_lpdf = 'wiener_lpdf(0.5, 1, 1.5, 0.5, 2)',
  eam_wiener_nodec = 'frm(bf(rt ~ 1, bias = 0.5), family = wiener(), data = dat_rt)',
  eam_wiener_badcol = 'frm(bf(rt | dec(one) ~ 1, bias = 0.5), family = wiener(), data = dat_rt)',
  ode_lincmt_parms = 'frm_lincmt(parms = list(1, 0.2, 10), times = 1:5, ncmt = 1, depot = TRUE)',
  ode_lincmt_missing = 'frm_lincmt(parms = list(ka = 1, ke = 0.2), times = 1:5, ncmt = 1, depot = TRUE)',
  spline_curve_nonfit = 'frm_curve(d, newdata = d)',
  spline_curve_level = 'frm_curve(fg, newdata = d, level = 1.5)',
  latent_hmm_reml = 'frm(bf(y ~ 1), family = hmm(K = 2, gaussian(), time = x, group = g), data = d, REML = TRUE)',
  latent_hmm_weights = 'frm(bf(y | weights(z) ~ 1), family = hmm(K = 2, gaussian(), time = x, group = g), data = d)',
  coupling_len = 'frm_cross_spectrum(p64x, p64y[1:10])',
  coupling_segments = 'frm_cross_spectrum(p64x, p64y, segments = 0)',
  coupling_whole = 'frm_cross_spectrum(p64x, p64y, segments = 2.5)',
  coupling_coherence = 'frm_coherence(fg)',
  learn_ctprobe = 'frm_task_simulate(1, d)',
  sample_loglik = 'log_lik(fg, pointwise = TRUE)',
  sample_nonfit = 'frm_sample(1)',
  # warnings
  w_single_level = 'invisible(frm(y ~ x + (1 | one), data = d, family = gaussian()))',
  w_obs_level = 'invisible(frm(y ~ x + (1 | obs), data = d, family = gaussian()))',
  w_sigma_varies = 'sigma(fsig)',
  w_nonfinite_cov = 'frmtmb:::solve_joint_precision(matrix(0, 2, 2))',
  w_structure = 'invisible(frmtmb_structure(loglik = function(y, dpars, aterms, block, extra) 1))',
  w_ps_span = 'invisible(predict(frm(y ~ ps(x, k = 6), data = d), newdata = data.frame(x = c(-50, 0, 50))))',
  w_two_in_call = '{ invisible(frm(y ~ x + (1 | one), data = d)); sigma(fsig) }',
  w_lapply = 'invisible(lapply(1:2, function(i) sigma(fsig)))',
  w_ode_states = 'invisible(frm_ode(function(t, y, p) lapply(y, function(v) -p[[1]] * v), init = as.list(rep(1, 8)), times = c(0.5, 1), parms = list(0.3), output = 1L))',
  w_sample_probs = 'NULL',
  w_ar_levels = 'invisible(frm(bf(y ~ x + ar(x, g, cov = TRUE)) + gaussian(), data = d, dry_run = "frame"))',
  w_coupling_df = 'NULL',
  w_under_warn2 = 'local({ op <- options(warn = 2); on.exit(options(op)); sigma(fsig) })',
  # messages
  m_na_rows = 'invisible(frm(y ~ w, data = d, family = gaussian()))',
  m_dates = 'invisible(frm(bf(y ~ day) + gaussian(), data = d, dry_run = "frame"))',
  m_rank_def = 'invisible(frm(y ~ x + I(2 * x), data = d))',
  m_bernoulli = 'invisible(frm(k | trials(n) ~ x, d, family = mixture(binomial, binomial), dry_run = "frame"))',
  m_categorical = 'invisible(frm(bf(ych ~ x), family = categorical(), data = transform(d, ych = sample(c("a", "b", "c"), 60, TRUE))))',
  m_anova_refit = 'invisible(anova(frm(y ~ x + (1 | g), data = d, REML = TRUE), frm(y ~ 1 + (1 | g), data = d, REML = TRUE), refit = TRUE))',
  m_confint_internal = 'invisible(confint(fm, "sd_g__Intercept"))',
  m_in_trycatch = 'tryCatch(invisible(frm(y ~ w, data = d)), error = function(e) "no")',
  m_nolf_notice = 'invisible(frm(y ~ x + (1 | g), data = d, family = gaussian()))'
)

wrap <- list(
  warn0 = function(e) sprintf("try(%s)", e),
  warn1 = function(e) sprintf("try(%s)", e),
  suppress = function(e) sprintf("suppressWarnings(suppressMessages(try(%s)))", e),
  catch = function(e) paste0(
    sprintf("print(tryCatch({%s; \"none\"}, error = function(c) paste(\"E\", class(c)[1L]), warning = function(c) paste(\"W\", conditionMessage(c)), message = function(c) paste(\"M\", conditionMessage(c))))\n", e),
    sprintf("try(withCallingHandlers(%s, warning = function(w) { cat(\"WCH-W:\", conditionMessage(w), \"\\n\"); invokeRestart(\"muffleWarning\") }, message = function(m) { cat(\"WCH-M:\", conditionMessage(m)); invokeRestart(\"muffleMessage\") }))", e)),
  record = function(e) NULL)

write_script <- function(arm, mode) {
  lines <- c(sprintf(".libPaths(c('%s', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))", libs[[arm]]),
             setup,
             if (mode == "warn1") "options(warn = 1)")
  if (mode == "record") {
    lines <- c(lines, "rec <- list()")
    for (nm in names(triggers)) {
      lines <- c(lines, sprintf(
        "rec[[%s]] <- tryCatch({%s; NULL}, condition = function(c) c(class = paste(class(c), collapse = '/'), message = conditionMessage(c), call = paste(deparse(conditionCall(c)), collapse = ' ')))",
        deparse(nm), triggers[[nm]]))
    }
    lines <- c(lines, sprintf("saveRDS(rec, '%s/record-%s.rds')", out, arm))
  } else {
    for (nm in names(triggers)) {
      lines <- c(lines,
                 sprintf("cat('\\n#### %s\\n', file = stderr()); cat('\\n#### %s\\n')", nm, nm),
                 wrap[[mode]](triggers[[nm]]))
    }
  }
  f <- file.path(out, sprintf("%s-%s.R", mode, arm))
  writeLines(lines, f)
  f
}

jobs <- expand.grid(arm = names(libs), mode = names(wrap),
                    stringsAsFactors = FALSE)
# REV_JOBS=lane:record runs one job, for checking the triggers first
sel <- Sys.getenv("REV_JOBS")
if (nzchar(sel)) {
  jobs <- jobs[paste(jobs$arm, jobs$mode, sep = ":") %in%
                 strsplit(sel, ",")[[1L]], ]
}
for (i in seq_len(nrow(jobs))) {
  f <- write_script(jobs$arm[i], jobs$mode[i])
  both <- suppressWarnings(system2(rscript, shQuote(f), stdout = TRUE,
                                   stderr = TRUE))
  writeLines(both, sub("[.]R$", ".both.txt", f))
  if (jobs$mode[i] != "record") {
    err <- suppressWarnings(system2(rscript, shQuote(f), stdout = FALSE,
                                    stderr = TRUE))
    writeLines(err, sub("[.]R$", ".stderr.txt", f))
  }
  cat("ran", basename(f), "\n")
}
