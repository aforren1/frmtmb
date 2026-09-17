# Reviewer, lane wt-conditions: user calls that reach a match.arg() or
# stopifnot() refusal inside a frmtmb function, and whether
# tryCatch(frmtmb_error = ) catches them on the lane build.
#   Rscript dev/conditions-rev-unclassed-probe.R      (seed 20260917)
.libPaths(c("C:/Users/adf44/source/r/conditions-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb); library(frmtmb.spline); library(frmtmb.eam)
  library(frmtmb.latent); library(frmtmb.ode); library(frmtmb.coupling)
  library(frmtmb.learn)
})
set.seed(20260917)
d <- data.frame(y = rnorm(40), x = rnorm(40), g = gl(8, 5))
fit <- frm(y ~ x, data = d, family = gaussian())
probes <- list(
  # match.arg on a user-supplied argument
  predict_type = quote(predict(fit, type = "bogus")),
  fitted_scale = quote(fitted(fit, scale = "bogus")),
  residuals_type = quote(residuals(fit, type = "bogus")),
  confint_method = quote(confint(fit, method = "bogus")),
  drop1_test = quote(drop1(fit, test = "bogus")),
  ce_band = quote(conditional_effects(fit, band = "bogus")),
  ce_method = quote(conditional_effects(fit, method = "bogus")),
  hypothesis_method = quote(hypothesis(fit, "x = 0", method = "bogus")),
  hypothesis_scope = quote(hypothesis(fit, "x = 0", scope = "bogus")),
  control_nlev = quote(frmtmb_control(check_nlev_1 = "bogus")),
  default_prior_route = quote(default_prior(y ~ x, data = d,
                                            route = "bogus")),
  vcov_cluster_type = quote(vcov_cluster(fit, d$g, type = "bogus")),
  periodogram_taper = quote(frm_periodogram(rnorm(64), taper = "bogus")),
  cross_spectrum_window = quote(frm_cross_spectrum(rnorm(64), rnorm(64),
                                                   window = "bogus")),
  hmm_init = quote(hmm(2, init = "bogus")),
  # stopifnot on user input
  set_prior_normal0 = quote(set_prior("normal(0)")),
  set_prior_student3 = quote(set_prior("student_t(3, 0)")),
  set_prior_exponential = quote(set_prior("exponential(-1)")),
  set_prior_gamma1 = quote(set_prior("gamma(1)")),
  set_prior_numeric = quote(set_prior(1)),
  prior_plus_string = quote(set_prior("normal(0, 1)") + "x"),
  diagnose_nonfit = quote(diagnose(lm(y ~ x, d))),
  vcov_cluster_nonfit = quote(vcov_cluster(lm(y ~ x, d), d$g)),
  frm_prior_list = quote(frm(y ~ x, data = d, prior = list(1))),
  frmtmb_family_bad = quote(frmtmb_family(family = 1, dpars = "mu",
                                          lpdf = function(y, mu) mu)),
  # the control: a converted refusal
  control_converted = quote(frm(y ~ x, data = d, family = "not_a_family"))
)
res <- lapply(names(probes), function(nm) {
  e <- tryCatch(eval(probes[[nm]]), error = identity,
                warning = identity)
  caught <- tryCatch({eval(probes[[nm]]); "no error"},
                     frmtmb_error = function(e) "caught",
                     error = function(e) "NOT caught")
  data.frame(probe = nm,
             class = if (inherits(e, "condition")) class(e)[1] else "none",
             frmtmb_error = caught,
             message = if (inherits(e, "condition"))
               substr(gsub("\\s+", " ", conditionMessage(e)), 1, 70) else "",
             stringsAsFactors = FALSE)
})
out <- do.call(rbind, res)
options(width = 200)
print(out, right = FALSE)
cat("\nprobes", nrow(out), " NOT caught by frmtmb_error",
    sum(out$frmtmb_error == "NOT caught"), "\n")
