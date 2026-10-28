# Reviewer check: readers of sigma the lane's table does not list, on a
# fresh AR(1) cov = TRUE fit (ar about 0.8), lane against rellib-r14.
# Usage: Rscript dev/review-sigmainn-readers.R after|before
arm <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sigmainn"
.libPaths(c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-sigmainn-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(wt, "dev", "stan-cache"))
suppressPackageStartupMessages({
  library(frmtmb)
})
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(20261028)
ng <- 20; k <- 8
d <- data.frame(t = rep(seq_len(k), ng), g = factor(rep(seq_len(ng), each = k)))
d$x <- rnorm(nrow(d)); d$x2 <- rnorm(nrow(d))
d$y <- 1 + 0.5 * d$x + unlist(lapply(seq_len(ng), function(i)
  as.numeric(arima.sim(list(ar = 0.8), k, sd = 0.5))))
fit <- frm(bf(y ~ x + ar(time = t, gr = g, cov = TRUE)), data = d)
ac <- fit$frame$autocor[[1]]
arv <- frmtmb:::autocor_natural(fit$estimates$thetaac[ac$theta_idx], ac)[[1]]
sg <- exp(fit$estimates$betad[[1]])
cat(sprintf("ar %.6f  sigma parameter %.6f  sigma/sqrt(1-ar^2) %.6f  sigma*sqrt(1-ar^2) %.6f\n",
            arv, sg, sg / sqrt(1 - arv^2), sg * sqrt(1 - arv^2)))
r <- function(lab, expr) {
  v <- tryCatch(force(expr), error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-40s %s\n", lab, if (is.character(v)) substr(v[1], 1, 140) else
    paste(signif(v, 6), collapse = " ")))
}
set.seed(1)
r("sigma()", sigma(fit))
r("insight::get_sigma", as.numeric(insight::get_sigma(fit)))
r("performance::r2", unlist(performance::r2(fit)))
r("performance::rmse", performance::rmse(fit))
r("predict() Est.Error mean (marg sd + se)", {
  p <- predict(fit, ndraws = 4000); mean(p[, "Est.Error"]) })
r("posterior_predict? sd of (draw - fitted)", {
  pp <- predict(fit, summary = FALSE, ndraws = 4000)
  sd(as.numeric(sweep(pp, 2, fitted(fit)))) })
r("conditional_effects(pp) band /3.92", {
  ce <- conditional_effects(fit, effects = "x", method = "posterior_predict",
                            ndraws = 4000)[[1]]
  mean((ce$upper__ - ce$lower__) / (2 * qnorm(0.975))) })
r("conditional_effects(dpar=sigma) est", {
  ce <- conditional_effects(fit, effects = "x", dpar = "sigma")[[1]]
  range(ce$estimate__) })
r("residuals(pearson) sd", sd(residuals(fit, type = "pearson")))
r("residuals(deviance) sd", sd(residuals(fit, type = "deviance")))
r("cooks.distance sum", sum(cooks.distance(fit)))
r("influence hat sum", sum(influence(fit)$hat))
r("dfbetas max abs", max(abs(dfbetas(fit))))
r("rstandard sd", sd(rstandard(fit)))
r("simulate sd(resid)", { s <- simulate(fit, nsim = 400, seed = 2)
  sd(as.matrix(s) - fitted(fit)) })
r("emmeans dpar sigma (log)", summary(emmeans::emmeans(fit, ~ 1, dpar = "sigma"))$emmean)
r("marginaleffects predictions(dpar? sigma)", {
  marginaleffects::avg_predictions(fit, dpar = "sigma")$estimate })
r("R2 / bayes_R2 (ML)", tryCatch(unlist(bayes_R2(fit)), error = function(e) stop(e)))
r("frm_simulate sd", {
  s <- frm_simulate(fit, nsim = 200, seed = 3); sd(as.matrix(s) - fitted(fit)) })
r("dharma? (skip)", NA)
