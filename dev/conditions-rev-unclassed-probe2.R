# Reviewer, lane wt-conditions (round 2): the user-reachable match.arg()
# and stopifnot() sites that dev/conditions-rev-unclassed-probe.R did not
# reach, plus other unclassed refusal shapes, on the lane build. A probe
# counts only when its message is the base R text of that site.
#   Rscript dev/conditions-rev-unclassed-probe2.R      (seed 20260917)
.libPaths(c("C:/Users/adf44/source/r/conditions-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb); library(frmtmb.spline); library(frmtmb.eam)
  library(frmtmb.latent); library(frmtmb.ode); library(frmtmb.coupling)
  library(frmtmb.learn); library(frmtmb.sample)
})
set.seed(20260917)
d <- data.frame(y = rnorm(40), x = rnorm(40), g = gl(8, 5))
fit <- frm(y ~ x, data = d, family = gaussian())
fake_draws <- structure(list(), class = "frmtmb_draws")
probes <- list(
  check_custom_family = quote(check_custom_family(1)),
  cluster_scores_nonfit = quote(cluster_scores(lm(y ~ x, d), d$g)),
  register_prior_defaults = quote(frmtmb_register_prior_defaults(1)),
  periodogram_detrend = quote(frm_periodogram(rnorm(64), detrend = "bogus")),
  anova_multiple_method = quote(frmtmb:::anova.frmtmb_multiple(
    structure(list(), class = "frmtmb_multiple"), method = "bogus")),
  gddm_control_tridiagonal = quote(gddm_control(tridiagonal = "bogus")),
  gddm_lapse = quote(gddm(lapse = "bogus")),
  gddm_simulate_lapse = quote(gddm_simulate(10, lapse = "bogus")),
  bandit2arm_dual_split = quote(bandit2arm_dual(subject = id,
                                                split = "bogus")),
  task_design_task = quote(frm_task_design(task = "bogus")),
  lincmt_output = quote(frm_lincmt(parms = list(ke = 0.2, V = 10),
                                   times = 1:5, output = "bogus")),
  ode_on_error = quote(frm_ode(function(t, y, p) list(-y[[1]]),
                               init = list(1), times = 1:3,
                               on_error = "bogus")),
  loo_compare_criterion = quote(loo_compare(fake_draws, fake_draws,
                                            criterion = "bogus")),
  hypothesis_draws_scope = quote(hypothesis(fake_draws, "x = 0",
                                            scope = "bogus")),
  pp_check_prefix = quote(pp_check(fake_draws, prefix = "bogus")),
  predictive_error_method = quote(predictive_error(fake_draws,
                                                   method = "bogus")),
  curve_feature_type = quote(frm_curve_feature(fit, var = "x",
                                               type = "bogus")),
  royston_parmar_scale = quote(royston_parmar(scale = "bogus")),
  rp_floored_action = quote(rp_floored(fit, action = "bogus")),
  # other unclassed shapes a user meets
  hmm_starts_nonfit = quote(hmm_starts(1)),
  frm_no_formula = quote(frm(data = d)),
  getME_mv = quote(lme4::getME(frm(mvbf(bf(y ~ x), bf(x ~ 1)), data = d),
                               "Zt"))
)
res <- lapply(names(probes), function(nm) {
  e <- tryCatch({eval(probes[[nm]]); NULL}, error = identity)
  data.frame(probe = nm,
             class = if (is.null(e)) "no error" else
               paste(class(e), collapse = "/"),
             message = if (is.null(e)) "" else
               substr(gsub("\\s+", " ", conditionMessage(e)), 1, 90),
             stringsAsFactors = FALSE)
})
out <- do.call(rbind, res)
options(width = 220)
print(out, right = FALSE)
cat("\nprobes", nrow(out), " not frmtmb_error",
    sum(!grepl("frmtmb_error", out$class)), "\n")
