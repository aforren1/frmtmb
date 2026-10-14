# Second pass on the lane build: summary rows, CE, RE fits, the sampler
# routes. Output: dev/famsimplex-try2.log
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/stan-cache")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-sim.R")
show <- function(lab, expr) {
  cat("\n---", lab, "\n")
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r) && length(r) == 1) cat(r, "\n") else {
    cat("class:", class(r), " dim:", dim(r), "\n")
    if (!is.null(dimnames(r))) str(dimnames(r))
  }
  invisible(r)
}
dl <- fs_data_logistic_normal(seed = 13)
fit2 <- frm(bf(Y ~ x, sigmab ~ z), family = logistic_normal(), data = dl)
print(fit2)
show("ln ce", conditional_effects(fit2, "x"))
show("ln default_prior sample", print(default_prior(fit2, route = "sample")))
show("ln hypothesis", print(hypothesis(fit2, "lncor__b__c > 0")))
dd <- fs_data_dirichlet(seed = 11)
fit <- frm(bf(Y ~ x, phi ~ z), family = dirichlet(), data = dd)
ce <- show("dir ce", conditional_effects(fit, "x"))
print(head(ce[[1]]))
show("dir ce boot", conditional_effects(fit, "x", band = "boot", boot = 20))
dre <- fs_data_dirichlet(seed = 12, n = 400, ng = 20, phi_formula = FALSE)
fre <- frm(bf(Y ~ x + (1 | ID | g)), family = dirichlet(), data = dre)
print(fre)
print(VarCorr(fre))
show("dirre fitted", fitted(fre))
show("dirre ce", conditional_effects(fre, "x"))
lre <- fs_data_logistic_normal(seed = 14, n = 400, ng = 20,
                               sigma_formula = FALSE)
flre <- frm(bf(Y ~ x + (1 | ID | g)), family = logistic_normal(), data = lre)
print(flre)
show("default_prior lnre sample", print(default_prior(flre, route = "sample")))
cat("\n=== sampler\n")
ds <- frm_sample(fit2, chains = 2, iter = 600, refresh = 0, seed = 1)
print(head(variables(ds), 30))
show("posterior_epred ln", posterior_epred(ds, ndraws = 5))
show("posterior_predict ln", pp <- posterior_predict(ds, ndraws = 5))
show("log_lik ln", log_lik(ds, ndraws = 5))
show("pp_check ln", pp_check(ds))
show("predictive_error ln", predictive_error(ds, ndraws = 5))
ds2 <- frm_sample(fit, chains = 2, iter = 600, refresh = 0, seed = 1)
show("posterior_epred dir", posterior_epred(ds2, ndraws = 5))
show("posterior_predict dir", pp2 <- posterior_predict(ds2, ndraws = 5))
if (is.array(pp2)) cat("row sums:", range(apply(pp2, 1:2, sum)), "\n")
show("posterior_linpred dir", posterior_linpred(ds2, ndraws = 5))
show("log_lik dir", log_lik(ds2, ndraws = 5))
show("ce dir draws", conditional_effects(ds2, "x"))
show("loo dir", loo(ds2))
show("as_tmbstan", { a <- as_tmbstan(fit, chains = 1, iter = 300,
                                      refresh = 0); class(a) })
