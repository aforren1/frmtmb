# brms 2.23.0 fits of the simplex-response families on the data the
# agreement check reuses, and the shapes of every post-fit method.
# Output: dev/famsimplex-brms-fits2-<model>.log, dev/famsimplex-brms-<model>.rds
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
source("C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-sim.R")
out_dir <- "C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev"
args <- commandArgs(TRUE)
which <- if (length(args)) args else c("dir", "dirre", "ln", "lnre", "dir2")
shape <- function(lab, expr) {
  cat("\n---", lab, "\n")
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r) && length(r) == 1) cat(r, "\n") else {
    cat("class:", class(r), " dim:", dim(r), "\n")
    if (!is.null(dimnames(r))) str(dimnames(r))
  }
  invisible(r)
}
run <- function(name, formula, family, data, seed = 1) {
  cat("\n==========", name, "==========\n")
  t0 <- proc.time()
  fit <- brm(formula, data = data, family = family, chains = 4,
             cores = 4, iter = 2000, seed = seed, refresh = 0,
             backend = "rstan")
  cat("time:", (proc.time() - t0)[["elapsed"]], "\n")
  saveRDS(fit, file.path(out_dir, paste0("famsimplex-brmsfit-", name,
                                         ".rds")))
  print(summary(fit))
  dr <- as_draws_df(fit)
  vars <- grep("^(b_|sd_|cor_|phi|sigma|lncor|Intercept)",
               posterior::variables(dr), value = TRUE)
  summ <- posterior::summarise_draws(
    posterior::subset_draws(dr, variable = vars),
    "median", "sd", "mean", "rhat", "ess_bulk")
  print(as.data.frame(summ), digits = 6)
  saveRDS(list(summ = as.data.frame(summ), data = data,
               formula = formula, family = family$family,
               variables = posterior::variables(dr)),
          file.path(out_dir, paste0("famsimplex-brms-", name, ".rds")))
  cat("variables:", paste(posterior::variables(dr), collapse = " "), "\n")
  shape("fitted", fitted(fit))
  shape("fitted summary=FALSE", fitted(fit, summary = FALSE))
  shape("posterior_epred", posterior_epred(fit, ndraws = 5))
  shape("predict", predict(fit, ndraws = 50))
  pp <- shape("posterior_predict", posterior_predict(fit, ndraws = 5))
  if (is.array(pp)) cat("pp row sums:", range(apply(pp, 1:2, sum)), "\n")
  shape("posterior_linpred", posterior_linpred(fit, ndraws = 5))
  shape("posterior_linpred dpar", posterior_linpred(fit, ndraws = 5,
                                                    dpar = "mub"))
  ll <- shape("log_lik", log_lik(fit, ndraws = 5))
  shape("residuals", residuals(fit))
  shape("predictive_error", predictive_error(fit, ndraws = 5))
  ce <- shape("conditional_effects", conditional_effects(fit, "x"))
  if (is.list(ce)) print(head(ce[[1]]))
  ce2 <- shape("conditional_effects categorical",
        conditional_effects(fit, "x", categorical = TRUE))
  if (is.list(ce2)) print(head(ce2[[1]]))
  shape("pp_check", pp_check(fit, ndraws = 5))
  shape("loo", loo(fit))
  shape("bayes_R2", bayes_R2(fit))
  shape("VarCorr", VarCorr(fit))
  ps <- shape("prior_summary", prior_summary(fit))
  print(ps)
  invisible(fit)
}
if ("dir" %in% which) {
  dd <- fs_data_dirichlet(seed = 11)
  run("dir", bf(Y ~ x, phi ~ z), dirichlet(), dd)
}
if ("dirre" %in% which) {
  dd <- fs_data_dirichlet(seed = 12, n = 400, ng = 20, phi_formula = FALSE)
  run("dirre", bf(Y ~ x + (1 | ID | g)), dirichlet(), dd)
}
if ("ln" %in% which) {
  dd <- fs_data_logistic_normal(seed = 13)
  run("ln", bf(Y ~ x, sigmab ~ z), logistic_normal(), dd)
}
if ("lnre" %in% which) {
  dd <- fs_data_logistic_normal(seed = 14, n = 400, ng = 20,
                                sigma_formula = FALSE)
  run("lnre", bf(Y ~ x + (1 | ID | g)), logistic_normal(), dd)
}
if ("dir2" %in% which) {
  dd <- fs_data_dirichlet2(seed = 15)
  run("dir2", bf(Y ~ x), brmsfamily("dirichlet2"), dd)
}
