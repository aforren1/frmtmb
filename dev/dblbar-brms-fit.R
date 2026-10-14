# brms 2.23.0 fits of three `||` models on lane dblbar's data (seed
# 20261013), 4 chains of 2000 iterations, sampler seed 1. Saves the
# draws' medians and sds, and the layout of VarCorr(), ranef(), coef()
# and summary(), to dev/dblbar-log/brms-fit-<case>.rds.
# Usage: Rscript dev/dblbar-brms-fit.R <case>
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
source(file.path(wt, "dev/dblbar-cases.R"))
suppressPackageStartupMessages(library(brms))
case <- commandArgs(TRUE)[1]
dd <- dblbar_data()
cs <- dblbar_cases[[case]]
form <- eval(str2lang(cs[1]))
fit <- brm(form, data = dd, family = eval(str2lang(cs[2])), chains = 4,
           iter = 2000, seed = 1, refresh = 0, cores = 4,
           backend = "rstan")
dr <- as_draws_df(fit)
v <- setdiff(variables(fit), c("lprior", "lp__"))
v <- v[!grepl("^(r_|z_|L_)", v)]
m <- as.matrix(dr[, v])
nd <- data.frame(x = c(0.3, -0.5), z = c(0.2, 0.1), f = c("b", "c"),
                 g = c("NEW1", "NEW2"), g1 = "p1", hg = "u")
pe_unc <- posterior_epred(fit, newdata = nd, allow_new_levels = TRUE,
                          sample_new_levels = "uncertainty", seed = 3)
pe_gau <- posterior_epred(fit, newdata = nd, allow_new_levels = TRUE,
                          sample_new_levels = "gaussian", seed = 3)
pe_pop <- posterior_epred(fit, newdata = nd, re_formula = NA)
out <- list(
  case = cs, variables = variables(fit),
  median = apply(m, 2, stats::median), sd = apply(m, 2, stats::sd),
  rhat = posterior::summarise_draws(dr[, v])$rhat,
  vc = VarCorr(fit), ranef_dimnames = lapply(ranef(fit), dimnames),
  coef_dimnames = lapply(coef(fit), dimnames),
  summary_random = lapply(summary(fit)$random, rownames),
  newlev = list(unc_mean = colMeans(pe_unc), unc_sd = apply(pe_unc, 2, sd),
                gau_mean = colMeans(pe_gau), gau_sd = apply(pe_gau, 2, sd),
                pop_mean = colMeans(pe_pop), pop_sd = apply(pe_pop, 2, sd)),
  re_formula_int = colMeans(posterior_epred(
    fit, newdata = dd[1:6, ], re_formula = ~ (1 | g)))
)
saveRDS(out, file.path(wt, sprintf("dev/dblbar-log/brms-fit-%s.rds", case)))
print(out[c("median", "sd", "rhat", "summary_random", "newlev")])
str(out$vc)
