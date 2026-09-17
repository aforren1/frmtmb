# Reviewer recheck round 1: does a duplicate-slot prior reach a fit
# through an entry point other than frm()? Lane build, seed 20260916.
#   Rscript dev/priorform-rev2-entry.R
.libPaths(c("C:/Users/adf44/source/r/priorform-lib", "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(20260916)
n <- 120
d <- data.frame(x = rnorm(n), g = factor(rep(1:12, 10)))
d$y <- 1 + 0.5 * d$x + rnorm(12)[d$g] + rnorm(n)
dup <- set_prior("normal(0, 1)", coef = "x") + set_prior("", coef = "x", lb = 0.9)
fit <- frm(y ~ x + (1 | g), data = d)
try1 <- function(label, expr) {
  r <- tryCatch({force(expr); "ACCEPTED"},
                error = function(e) paste("refused:", substr(conditionMessage(e), 1, 60)))
  cat(sprintf("%-26s %s\n", label, r))
}
try1("frm", frm(y ~ x + (1 | g), data = d, prior = dup))
try1("update(prior =)", update(fit, prior = dup))
try1("validate_prior", validate_prior(dup, y ~ x + (1 | g), data = d))
try1("frm_simulate", frm_simulate(y ~ x + (1 | g), d, nsim = 2, seed = 1, prior = set_prior("normal(0, 1)", class = "b") + set_prior("normal(0, 2)", class = "b") + set_prior("exponential(1)", class = "sd") + set_prior("exponential(1)", class = "sigma") + set_prior("normal(0, 1)", class = "Intercept")))
try1("par_template", par_template(y ~ x + (1 | g), data = d, prior = dup))
try1("hypothesis/prior_summary", NULL)
fitp <- frm(y ~ x + (1 | g), data = d, prior = set_prior("normal(0, 1)", coef = "x"))
try1("update(fitp, prior=dup)", update(fitp, prior = dup))
try1("frm_bootstrap(update fn)", frm_bootstrap(fitp, nboot = 2, seed = 1))
try1("frm_sample(fit, prior=dup) pre", tryCatch(frm_sample(fit, prior = dup, chains = 1, iter = 10, refresh = 0), error = function(e) stop(e)))
# a MAP fit's own prior and the call's prior on the same slot: brms's
# update() keeps the new row and drops the old one
fitc <- frm(y ~ x + (1 | g), data = d, prior = set_prior("normal(0, 1)", coef = "x"))
r <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fitc, set_prior("", coef = "x", lb = 0.9), base = fitc$prior))
cat("MAP normal(0,1) coef x + call bounds-only coef x ->", paste(vapply(r$ri$entries, function(e) paste0(e$comp, e$idx, "=", e$dist$kind), ""), collapse = " "),
    "| lower:", paste(names(r$ri$lower), r$ri$lower), "\n")
r2 <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fitc, set_prior("normal(0, 5)", class = "b"), base = fitc$prior))
cat("MAP coef x normal(0,1) + call class b normal(0,5) ->", paste(vapply(r2$ri$entries, function(e) paste0(e$comp, e$idx, "=", paste(unlist(e$dist), collapse="/")), ""), collapse = " "), "\n")
