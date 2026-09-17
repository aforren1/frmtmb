# Lane wt-conditions, punch round 1: every user-reachable refusal the
# review listed as unclassed (dev/reviews/20260917-conditions.md), one
# probe per site, including every distribution of parse_prior_dist()'s
# switch. Records the class vector, whether tryCatch(frmtmb_error = )
# catches it, and whether the extension subclass is the expected one.
#   Rscript dev/conditions-punch-probe.R <lib> <outfile>   (seed 20260917)
av <- commandArgs(trailingOnly = TRUE)
.libPaths(c(av[1L], "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb); library(frmtmb.spline); library(frmtmb.eam)
  library(frmtmb.latent); library(frmtmb.ode); library(frmtmb.coupling)
  library(frmtmb.learn); library(frmtmb.sample)
})
set.seed(20260917)
d <- data.frame(y = rnorm(40), x = rnorm(40), g = gl(8, 5),
                f = factor(rep(c("a", "b"), 20)))
fit <- frm(y ~ x, data = d, family = gaussian())
fitf <- frm(y ~ x + f, data = d, family = gaussian())
fake_draws <- structure(list(), class = "frmtmb_draws")
mv <- frm(mvbf(bf(y ~ x), bf(x ~ 1)), data = d)
short <- rnorm(3)
dl <- d
dl$lc <- replicate(40, list(1:2), simplify = FALSE)

# name = list(expression, package whose subclass the refusal should carry)
P <- function(e, pkg = "frmtmb") list(e = substitute(e), pkg = pkg)
probes <- list(
  # match.arg(), core
  predict_type = P(predict(fit, type = "bogus")),
  fitted_scale = P(fitted(fit, scale = "bogus")),
  residuals_type = P(residuals(fit, type = "bogus")),
  confint_method = P(confint(fit, method = "bogus")),
  drop1_test = P(drop1(fit, test = "bogus")),
  hypothesis_method = P(hypothesis(fit, "x = 0", method = "bogus")),
  hypothesis_scope = P(hypothesis(fit, "x = 0", scope = "bogus")),
  ce_band = P(conditional_effects(fit, band = "bogus")),
  ce_method = P(conditional_effects(fit, method = "bogus")),
  control_nlev = P(frmtmb_control(check_nlev_1 = "bogus")),
  control_olre = P(frmtmb_control(check_olre = "bogus")),
  default_prior_route = P(default_prior(y ~ x, data = d, route = "bogus")),
  get_prior_route = P(get_prior(y ~ x, data = d, route = "bogus")),
  validate_prior_route = P(validate_prior(set_prior("normal(0, 1)"),
                                          y ~ x, data = d,
                                          route = "bogus")),
  vcov_cluster_type = P(vcov_cluster(fit, ~ g, type = "bogus")),
  periodogram_taper = P(frm_periodogram(rnorm(64), taper = "bogus")),
  periodogram_detrend = P(frm_periodogram(rnorm(64), detrend = "bogus")),
  anova_multiple_method = P(frmtmb:::anova.frmtmb_multiple(
    structure(list(), class = "frmtmb_multiple"), method = "bogus")),
  anova_multiple_use = P(frmtmb:::anova.frmtmb_multiple(
    structure(list(), class = "frmtmb_multiple"), use = "bogus")),
  # match.arg(), extensions
  cross_spectrum_window = P(frm_cross_spectrum(rnorm(64), rnorm(64),
                                               window = "bogus"),
                            "frmtmb.coupling"),
  gddm_control_tridiagonal = P(gddm_control(tridiagonal = "bogus"),
                               "frmtmb.eam"),
  gddm_lapse = P(gddm(lapse = "bogus"), "frmtmb.eam"),
  gddm_simulate_lapse = P(gddm_simulate(10, lapse = "bogus"), "frmtmb.eam"),
  hmm_init = P(hmm(2, init = "bogus"), "frmtmb.latent"),
  bandit2arm_dual_split = P(bandit2arm_dual(subject = id, split = "bogus"),
                            "frmtmb.learn"),
  task_design_task = P(frm_task_design(task = "bogus"), "frmtmb.learn"),
  lincmt_output = P(frm_lincmt(parms = list(ke = 0.2, V = 10), times = 1:5,
                               output = "bogus"), "frmtmb.ode"),
  ode_on_error = P(frm_ode(function(t, y, p) list(-y[[1]]), init = list(1),
                           times = 1:3, on_error = "bogus"), "frmtmb.ode"),
  loo_compare_criterion = P(loo_compare(fake_draws, fake_draws,
                                        criterion = "bogus"),
                            "frmtmb.sample"),
  hypothesis_draws_scope = P(hypothesis(fake_draws, "x = 0",
                                        scope = "bogus"), "frmtmb.sample"),
  pp_check_prefix = P(pp_check(fake_draws, prefix = "bogus"),
                      "frmtmb.sample"),
  predictive_error_method = P(predictive_error(fake_draws, method = "bogus"),
                              "frmtmb.sample"),
  curve_feature_type = P(frm_curve_feature(fit, var = "x", type = "bogus"),
                         "frmtmb.spline"),
  royston_parmar_scale = P(royston_parmar(scale = "bogus"), "frmtmb.spline"),
  rp_floored_action = P(rp_floored(fit, action = "bogus"), "frmtmb.spline"),
  # stopifnot() on user input
  set_prior_numeric = P(set_prior(1)),
  set_prior_normal = P(set_prior("normal(0)")),
  set_prior_student_t = P(set_prior("student_t(3, 0)")),
  set_prior_cauchy = P(set_prior("cauchy(0)")),
  set_prior_exponential_arity = P(set_prior("exponential(1, 2)")),
  set_prior_exponential_rate = P(set_prior("exponential(-1)")),
  set_prior_lkj = P(set_prior("lkj()")),
  set_prior_logistic = P(set_prior("logistic(0)")),
  set_prior_gamma = P(set_prior("gamma(1)")),
  set_prior_inv_gamma = P(set_prior("inv_gamma(1)")),
  set_prior_beta = P(set_prior("beta(1)")),
  prior_plus_string = P(set_prior("normal(0, 1)") + "x"),
  frm_prior_list = P(frm(y ~ x, data = d, prior = list(1))),
  diagnose_nonfit = P(diagnose(lm(y ~ x, d))),
  frmtmb_family_bad = P(frmtmb_family(family = 1, dpars = "mu",
                                      links = list(mu = "identity"),
                                      lpdf = function(y, dpars) y)),
  check_custom_family = P(check_custom_family(1)),
  cluster_scores_nonfit = P(cluster_scores(lm(y ~ x, d), ~ g)),
  vcov_cluster_nonfit = P(vcov_cluster(lm(y ~ x, d), ~ g)),
  register_prior_defaults = P(frmtmb:::frmtmb_register_prior_defaults(1)),
  # other unclassed errors
  getME_mv = P(lme4::getME(mv, "Zt")),
  hmm_starts_nonfit = P(hmm_starts(1), "frmtmb.latent"),
  task_simulate_nonfamily = P(frm_task_simulate(1, d), "frmtmb.learn"),
  predict_missing_column = P(predict(fit, newdata = d[, c("y", "g")])),
  predict_new_level = P(predict(fitf, newdata = transform(
    d[1:4, ], f = factor(c("a", "zz", "a", "b"))))),
  frm_unknown_variable = P(frm(y ~ nope, data = d)),
  frm_data_not_df = P(frm(y ~ x, data = "notdf")),
  frm_nl_unknown = P(frm(bf(y ~ b0 * no_such_thing, b0 ~ 1, nl = TRUE),
                         data = d, dry_run = "frame",
                         start = list(beta = 1))),
  frm_nl_short = P(frm(bf(y ~ b0 * short, b0 ~ 1, nl = TRUE), data = d,
                       dry_run = "frame", start = list(beta = 1))),
  frm_list_column = P(frm(y ~ lc, data = dl)),
  # the control: a refusal converted before this round
  control_converted = P(frm(y ~ x, data = d, family = "not_a_family"))
)
res <- lapply(names(probes), function(nm) {
  e <- tryCatch({eval(probes[[nm]]$e); NULL}, error = identity)
  pkg <- probes[[nm]]$pkg
  own <- if (pkg == "frmtmb") "frmtmb_error" else
    paste0(gsub(".", "_", pkg, fixed = TRUE), "_error")
  cls <- if (is.null(e)) "no error" else paste(class(e), collapse = "/")
  data.frame(probe = nm,
             frmtmb_error = inherits(e, "frmtmb_error"),
             own_class = inherits(e, own) &&
               (pkg != "frmtmb" || !grepl("^frmtmb_[a-z]+_error",
                                          class(e)[1L])),
             class = cls,
             message = if (is.null(e)) "" else
               substr(gsub("\\s+", " ", conditionMessage(e)), 1, 110),
             stringsAsFactors = FALSE)
})
out <- do.call(rbind, res)
sink(av[2L])
options(width = 300)
print(out[, c("probe", "frmtmb_error", "own_class", "class")], right = FALSE)
cat("\n")
for (i in seq_len(nrow(out))) cat(sprintf("%-28s %s\n", out$probe[i],
                                          out$message[i]))
cat("\nprobes", nrow(out), " frmtmb_error", sum(out$frmtmb_error),
    " not frmtmb_error", sum(!out$frmtmb_error),
    " expected subclass", sum(out$own_class), "\n")
sink()
