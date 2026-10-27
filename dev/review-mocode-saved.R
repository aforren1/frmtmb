## Reviewer probe: every method I can reach on the 0.75.0 fit of
## y ~ mo(k), k in 1..4, that dev/mocode-saved.R saved
## (dev/mocode-out/old-k.rds), read with the lane library.
##   Rscript dev/review-mocode-saved.R
.libPaths(c("C:/Users/adf44/source/r/wt-mocode-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "lib:",
    dirname(dirname(system.file(package = "frmtmb"))), "\n")
fit <- readRDS(
  "C:/Users/adf44/source/r/frmtmb-wt-mocode/dev/mocode-out/old-k.rds")
run <- function(lab, expr) {
  w <- character(0)
  r <- withCallingHandlers(
    tryCatch(suppressMessages(expr), error = function(e) e),
    warning = function(c) {
      w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
    })
  out <- if (inherits(r, "error")) {
    paste("REFUSED/ERROR:", substr(gsub("\n", " ", conditionMessage(r)),
                                    1, 110))
  } else {
    v <- tryCatch(as.numeric(unlist(r)), error = function(e) NA,
                  warning = function(e) NA)
    v <- v[is.finite(v) | is.nan(v)]
    paste("RETURNED:", paste(signif(head(v, 5), 6), collapse = " "))
  }
  cat(sprintf("%-30s %s\n", lab, out))
  if (length(w)) cat(sprintf("%-30s   warning: %s\n", "",
                             substr(unique(w)[1L], 1, 110)))
}
d <- fit$data
run("print()", utils::capture.output(print(fit)))
run("summary()", summary(fit)$fixed)
run("fixef()", fixef(fit))
run("coef()", coef(fit))
run("vcov()", vcov(fit))
run("logLik()", logLik(fit))
run("AIC()", AIC(fit))
run("BIC()", BIC(fit))
run("nobs()", nobs(fit))
run("confint()", confint(fit))
run("log_lik()", log_lik(fit))
run("loo()", loo(fit))
run("waic()", waic(fit))
run("fitted()", fitted(fit))
run("predict()", predict(fit))
run("posterior_epred()", posterior_epred(fit))
run("posterior_predict()", posterior_predict(fit, ndraws = 5))
run("posterior_linpred()", posterior_linpred(fit, ndraws = 5))
run("pp_check()", pp_check(fit, ndraws = 5))
run("hypothesis()", hypothesis(fit, "mok > 0"))
run("prior_summary()", prior_summary(fit))
run("variables()", variables(fit))
run("diagnose()", diagnose(fit))
run("frm_bootstrap()", frm_bootstrap(fit, fixef, nsim = 2, seed = 1))
run("drop1()", drop1(fit))
run("anova(old, refit)", anova(fit, update(fit)))
run("emmeans(~ k)", if (requireNamespace("emmeans", quietly = TRUE))
  summary(emmeans::emmeans(fit, ~ k)))
run("marginaleffects predictions",
    if (requireNamespace("marginaleffects", quietly = TRUE))
      marginaleffects::predictions(fit))
run("bayes_R2()", bayes_R2(fit))
run("loo_R2()", loo_R2(fit))
run("profile confint",
    confint(fit, method = "profile", parm = "mok"))
run("update(newdata = d)", fixef(update(fit, newdata = d)))
run("refit(fit, y)", fixef(refit(fit, d$y)))
