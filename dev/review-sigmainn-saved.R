# Reviewer check: fits and draws saved by frmtmb 0.76.0 (rellib-r14),
# read by the lane build. "save" writes them; "read" runs every reader
# and says which return a number and which refuse.
# Usage: Rscript dev/review-sigmainn-saved.R save|read|readr14
arm <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sigmainn"
out <- file.path(wt, "dev", "review-sigmainn-out")
.libPaths(c(if (identical(arm, "read")) "C:/Users/adf44/source/r/wt-sigmainn-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(wt, "dev", "stan-cache"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "; frmtmb.sample",
    format(packageVersion("frmtmb.sample")), "\n")
set.seed(20261028)
ng <- 20; k <- 8
d <- data.frame(t = rep(seq_len(k), ng), g = factor(rep(seq_len(ng), each = k)))
d$x <- rnorm(nrow(d)); d$x2 <- rnorm(nrow(d))
d$y <- 1 + 0.5 * d$x + unlist(lapply(seq_len(ng), function(i)
  as.numeric(arima.sim(list(ar = 0.8), k, sd = 0.5))))
fo <- list(
  covT = bf(y ~ x + ar(time = t, gr = g, cov = TRUE)),
  covTsx = bf(y ~ x + ar(time = t, gr = g, cov = TRUE), sigma ~ x2),
  covF = bf(y ~ x + ar(time = t, gr = g)),
  cosy = bf(y ~ x + cosy(time = t, gr = g)),
  plain = bf(y ~ x))
if (arm == "save") {
  fits <- lapply(fo, function(f) frm(f, data = d))
  saveRDS(fits, file.path(out, "saved-r14-fits.rds"))
  dr <- frm_sample(fo$covT, data = d, chains = 2, iter = 600, seed = 7,
                   refresh = 0)
  saveRDS(dr, file.path(out, "saved-r14-draws.rds"))
  cat("saved\n")
  quit(save = "no")
}
fits <- readRDS(file.path(out, "saved-r14-fits.rds"))
dr <- readRDS(file.path(out, "saved-r14-draws.rds"))
show <- function(v) {
  if (is.numeric(v) && length(v)) return(sprintf("value %s", paste(signif(head(as.numeric(v), 3), 6), collapse = " ")))
  if (is.data.frame(v) || is.list(v)) {
    num <- unlist(Filter(is.numeric, if (is.data.frame(v)) as.list(v) else v))
    return(sprintf("value <%s> %s", class(v)[1], paste(signif(head(num, 3), 6), collapse = " ")))
  }
  sprintf("value <%s>", class(v)[1])
}
try1 <- function(lab, expr) {
  r <- tryCatch(suppressWarnings(show(force(expr))), error = function(e)
    paste("REFUSED:", substr(gsub("\n", " ", conditionMessage(e)), 1, 110)))
  cat(sprintf("  %-34s %s\n", lab, r))
}
nd <- d[1:4, ]
for (m in names(fits)) {
  f <- fits[[m]]
  cat("\n==", m, "\n")
  try1("print", capture.output(print(f)))
  try1("summary", summary(f))
  try1("sigma", sigma(f))
  try1("VarCorr", VarCorr(f))
  try1("fixef", fixef(f))
  try1("coef", coef(f))
  try1("vcov", vcov(f, full = TRUE))
  try1("confint", confint(f))
  try1("logLik", logLik(f))
  try1("AIC", AIC(f))
  try1("fitted", fitted(f))
  try1("fitted(dpar=sigma)", fitted(f, dpar = "sigma"))
  try1("frm_linpred(dpar=sigma)", frm_linpred(f, dpar = "sigma"))
  try1("predict()", predict(f))
  try1("predict(newdata)", predict(f, newdata = nd))
  try1("predict(type=response,se)", predict(f, type = "response", se.fit = TRUE)$se.fit)
  try1("residuals(pearson)", residuals(f, type = "pearson"))
  try1("residuals(response)", residuals(f, type = "response"))
  try1("simulate", simulate(f, nsim = 2, seed = 1)[[1]])
  try1("cond_effects(dpar=sigma)", conditional_effects(f, dpar = "sigma"))
  try1("cond_effects(pp)", conditional_effects(f, method = "posterior_predict"))
  try1("emmeans(dpar=sigma)", summary(emmeans::emmeans(f, ~ 1, dpar = "sigma")))
  try1("emmeans(mu)", summary(emmeans::emmeans(f, ~ 1)))
  try1("marginaleffects::predictions", marginaleffects::predictions(f))
  try1("insight::get_sigma", insight::get_sigma(f))
  try1("insight::get_parameters", insight::get_parameters(f))
  try1("hypothesis", hypothesis(f, "x > 0"))
  try1("prior_summary", prior_summary(f))
  try1("posterior_predict", posterior_predict(f, ndraws = 2))
  try1("log_lik", log_lik(f))
  try1("loo", loo(f))
  try1("cooks.distance", cooks.distance(f))
  try1("update", logLik(update(f)))
}
cat("\n== draws (frm_sample 0.23.0, covT)\n")
try1("print", capture.output(print(dr)))
try1("summary", summary(dr))
try1("as.matrix sigma", as.matrix(dr)[, "sigma"])
try1("sigma(draws)", sigma(dr))
try1("posterior_predict", posterior_predict(dr, ndraws = 2))
try1("log_lik", log_lik(dr))
try1("loo", loo(dr))
try1("posterior_epred", posterior_epred(dr, ndraws = 2))
try1("bayes_R2", bayes_R2(dr))
try1("VarCorr", VarCorr(dr))
