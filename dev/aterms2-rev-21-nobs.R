# Reviewer, re-check after punch round 1: every consumer of nobs() and
# of the likelihood's row count on a UNIVARIATE subset() model `f`,
# against the same model fitted on d[s, ] (`g`). nobs(f) is now the
# data's rows (brms); everything that needs the FITTED rows must still
# get them. Seeds: data 2101, simulate 5, sampler 3.
# Log: dev/aterms2-rev-log-21-nobs.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
options(mc.cores = 1)
q <- function(e) suppressWarnings(suppressMessages(e))
tr <- function(e) tryCatch(q(e), error = function(e) paste("ERROR:",
                                                          conditionMessage(e)))
same <- function(label, a, b) {
  a <- tr(a)
  b <- tr(b)
  num <- is.numeric(unlist(a)) && is.numeric(unlist(b)) &&
    length(unlist(a)) == length(unlist(b))
  rd <- if (num) {
    x <- as.numeric(unlist(a)); y <- as.numeric(unlist(b))
    max(abs(x - y)) / max(1e-300, max(abs(y)))
  } else NA
  cat(sprintf("%-34s f: %-28s g: %-28s rel %s\n", label,
              substr(paste(format(unlist(a)[seq_len(min(3,
                length(unlist(a))))]), collapse = " "), 1, 28),
              substr(paste(format(unlist(b)[seq_len(min(3,
                length(unlist(b))))]), collapse = " "), 1, 28),
              if (is.na(rd)) "-" else format(rd, digits = 3)))
  invisible(list(a, b))
}
set.seed(2101)
n <- 80
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n),
                g = factor(rep(1:8, length.out = n)),
                s = rep(c(TRUE, TRUE, FALSE, FALSE), length.out = n))
d$y <- 1 + d$x + 0.5 * d$z + rnorm(8, 0, 0.5)[d$g] + rnorm(n)
# NA outside the subset: harmless, the row stays in the data's count;
# NA inside it: the row is dropped from both counts
d$x[which(!d$s)[1:3]] <- NA
d$x[which(d$s)[5]] <- NA
ds_ <- d[d$s, ]
f <- q(frm(y | subset(s) ~ x + z + w + (1 | g), data = d))
g <- q(frm(y ~ x + z + w + (1 | g), data = ds_))
cat("rows: data", n, " data after NA rule", n - 1, " fitted",
    length(f$frame$y[[1]]), "\n\n")
same("nobs()", nobs(f), nobs(g))
same("nobs(resp = 'y')", nobs(f, resp = "y"), nobs(g, resp = "y"))
same("logLik nobs attr", attr(logLik(f), "nobs"), attr(logLik(g), "nobs"))
same("logLik", as.numeric(logLik(f)), as.numeric(logLik(g)))
same("AIC", AIC(f), AIC(g))
same("BIC", BIC(f), BIC(g))
same("df.residual", df.residual(f), df.residual(g))
same("extractAIC", extractAIC(f), extractAIC(g))
same("extractAIC k = log(n)", extractAIC(f, k = log(nobs(f))),
     extractAIC(g, k = log(nobs(g))))
same("AIC(f, f2) data.frame", AIC(f, q(update(f, . ~ . - w))),
     AIC(g, q(update(g, . ~ . - w))))
same("BIC(f, f2) data.frame", BIC(f, q(update(f, . ~ . - w))),
     BIC(g, q(update(g, . ~ . - w))))
same("anova(f2, f) LRT", anova(q(update(f, . ~ . - w)), f),
     anova(q(update(g, . ~ . - w)), g))
same("drop1", drop1(f), drop1(g))
same("drop1 test = Chisq", drop1(f, test = "Chisq"),
     drop1(g, test = "Chisq"))
st_f <- tr(stats::step(f, trace = 0))
st_g <- tr(stats::step(g, trace = 0))
cat(sprintf("%-34s f: %s\n%-34s g: %s\n", "step() formula",
            if (is.character(st_f)) st_f else deparse1(formula(st_f)),
            "", if (is.character(st_g)) st_g else deparse1(formula(st_g))))
st_f <- tr(stats::step(f, trace = 0, k = log(nobs(f))))
st_g <- tr(stats::step(g, trace = 0, k = log(nobs(g))))
cat(sprintf("%-34s f: %s\n%-34s g: %s\n", "step(k = log(nobs)) formula",
            if (is.character(st_f)) st_f else deparse1(formula(st_f)),
            "", if (is.character(st_g)) st_g else deparse1(formula(st_g))))
sa_f <- tr(MASS::stepAIC(f, trace = 0))
sa_g <- tr(MASS::stepAIC(g, trace = 0))
cat(sprintf("%-34s f: %s\n%-34s g: %s\n", "MASS::stepAIC formula",
            if (is.character(sa_f)) sa_f else deparse1(formula(sa_f)),
            "", if (is.character(sa_g)) sa_g else deparse1(formula(sa_g))))
same("simulate(seed = 5)", simulate(f, nsim = 2, seed = 5),
     simulate(g, nsim = 2, seed = 5))
same("dharma_residuals scaled", dharma_residuals(f, seed = 5)$scaledResiduals,
     dharma_residuals(g, seed = 5)$scaledResiduals)
same("summary()$nobs", summary(f)$nobs, summary(g)$nobs)
same("vcov_cluster(~ g)", vcov_cluster(f, cluster = ~ g),
     vcov_cluster(g, cluster = ~ g))
same("vcov_cluster(fitted-row vector)",
     vcov_cluster(f, cluster = ds_$g[!is.na(ds_$x)]),
     vcov_cluster(g, cluster = ds_$g[!is.na(ds_$x)]))
same("vcov_cluster(data-row vector)", vcov_cluster(f, cluster = d$g),
     vcov_cluster(g, cluster = ds_$g[!is.na(ds_$x)]))
same("influence cooks", cooks.distance(influence(f)),
     cooks.distance(influence(g)))
same("frm_bootstrap(R = 4, seed 1)",
     frm_bootstrap(f, R = 4, seed = 1)$t, frm_bootstrap(g, R = 4, seed = 1)$t)
same("insight::n_obs", insight::n_obs(f), insight::n_obs(g))
same("insight::get_df residual", insight::get_df(f, type = "residual"),
     insight::get_df(g, type = "residual"))
same("insight::get_loglikelihood nobs",
     attr(insight::get_loglikelihood(f), "nobs"),
     attr(insight::get_loglikelihood(g), "nobs"))
same("emmeans df", summary(emmeans::emmeans(f, ~ 1))$df,
     summary(emmeans::emmeans(g, ~ 1))$df)
same("sandwich::sandwich", sandwich::sandwich(f), sandwich::sandwich(g))

cat("\n== frmtmb.sample: same model, same seed\n")
fs <- q(frm(y | subset(s) ~ x + z + w, data = d))
gs <- q(frm(y ~ x + z + w, data = ds_))
fs2 <- q(frm(y | subset(s) ~ x + z, data = d))
gs2 <- q(frm(y ~ x + z, data = ds_))
smp <- function(m) q(frm_sample(m, chains = 1, iter = 400, refresh = 0,
                                seed = 3))
dfs <- smp(fs); dgs <- smp(gs); dfs2 <- smp(fs2); dgs2 <- smp(gs2)
same("draws identical (max rel)", dfs$draws, dgs$draws)
same("nobs(draws)", nobs(dfs), nobs(dgs))
lf <- q(loo(dfs)); lg <- q(loo(dgs)); lf2 <- q(loo(dfs2)); lg2 <- q(loo(dgs2))
same("loo estimates", lf$estimates, lg$estimates)
same("loo pointwise dim", dim(lf$pointwise), dim(lg$pointwise))
same("waic estimates", q(waic(dfs))$estimates, q(waic(dgs))$estimates)
same("loo_compare", loo::loo_compare(lf, lf2), loo::loo_compare(lg, lg2))
same("loo_model_weights stacking",
     loo::loo_model_weights(list(lf, lf2), method = "stacking"),
     loo::loo_model_weights(list(lg, lg2), method = "stacking"))
same("loo_model_weights pseudobma",
     loo::loo_model_weights(list(lf, lf2), method = "pseudobma", BB = FALSE),
     loo::loo_model_weights(list(lg, lg2), method = "pseudobma", BB = FALSE))
same("bayes_R2", bayes_R2(dfs), bayes_R2(dgs))
same("log_lik dim", dim(log_lik(dfs)), dim(log_lik(dgs)))

cat("\n== multivariate subset model\n")
d$s2 <- rep(c(TRUE, FALSE, TRUE, TRUE, TRUE), length.out = n)
d$y2 <- d$z + rnorm(n)
m <- q(frm(bf(y | subset(s) ~ z) + bf(y2 | subset(s2) ~ w), data = d,
           family = gaussian()))
cat("nobs()", nobs(m), " nobs(resp = y)", nobs(m, resp = "y"),
    " nobs(resp = y2)", nobs(m, resp = "y2"), " logLik nobs",
    attr(logLik(m), "nobs"), " df.residual", df.residual(m),
    " rows y", length(m$frame$y$y), " y2", length(m$frame$y$y2), "\n")
cat("nobs(m, 'y') positional:", tr(nobs(m, "y")), "\n")
cat("nobs(m, resp = 'nope'):", tr(nobs(m, resp = "nope")), "\n")
cat("nobs(m, resp = c('y','y2')):", tr(nobs(m, resp = c("y", "y2"))), "\n")
cat("nobs(m, 7):", tr(nobs(m, 7)), "\n")
dm <- smp(m)
cat("nobs(draws, resp = 'y2'):", tr(nobs(dm, resp = "y2")),
    " nobs(draws):", tr(nobs(dm)), "\n")
