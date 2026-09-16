# Reviewer's output harness. The lane hashed VALUE and WARNINGS on 22
# calls and muffled messages. This one also hashes MESSAGES and the
# PRINTED representation, and widens the call list to things a brms
# user actually runs on a fit: named and positional argument forms,
# the summary and print surface, the plotting methods built rather
# than counted, and the core names frmtmb.sample re-exports, which the
# lane's list left out entirely.
#
#   Rscript dev/sgrev-brmsout.R <LIB> <arm> <outfile>
av <- commandArgs(trailingOnly = TRUE)
.libPaths(c(av[1], "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
arm <- av[2]
q <- function(e) suppressWarnings(suppressMessages(e))
if (arm == "none") q(library(brms))
if (arm == "S") { q(library(brms)); q(library(frmtmb.sample)) }
if (arm == "T") { q(library(frmtmb.sample)); q(library(brms)) }
if (arm == "U") { q(library(frmtmb.sample)); q(loadNamespace("brms")) }
if (arm == "SC") { q(library(brms)); q(library(frmtmb))
                   q(library(frmtmb.sample)) }
q(library(loo)); q(library(posterior)); q(library(bayesplot))
q(library(coda)); q(library(bridgesampling)); q(library(rstantools))

bfit <- readRDS("dev/stan-cache/samplegen-brmsfit.rds")
nd <- bfit$data[1:5, , drop = FALSE]

calls <- alist(
  # --- the lane's 22, for comparability -----------------------------
  as.mcmc = as.mcmc(bfit),
  log_lik = log_lik(bfit, ndraws = 40),
  log_posterior = log_posterior(bfit),
  loo_subsample = loo_subsample(bfit, observations = 50),
  mcmc_plot_n = nrow(mcmc_plot(bfit)$data),
  neff_ratio = neff_ratio(bfit),
  nsamples = nsamples(bfit),
  nuts_params = nuts_params(bfit),
  parnames = parnames(bfit),
  posterior_epred = posterior_epred(bfit, ndraws = 40),
  posterior_interval = posterior_interval(bfit),
  posterior_linpred = posterior_linpred(bfit, ndraws = 40),
  posterior_predict = posterior_predict(bfit, ndraws = 40),
  posterior_samples = posterior_samples(bfit),
  pp_mixture = pp_mixture(bfit),
  predictive_error = predictive_error(bfit, ndraws = 40),
  predictive_interval = predictive_interval(bfit, ndraws = 40),
  psis = psis(bfit),
  restructure = list(class(restructure(bfit)),
                     dim(as.matrix(restructure(bfit)))),
  rhat = rhat(bfit),
  stancode = stancode(bfit),
  standata = standata(bfit),
  # --- forms the lane did not run -----------------------------------
  # POSITIONAL second arguments, which is where a method whose
  # formals differ from brms's answers a different question
  as_mcmc_pos = as.mcmc(bfit, "b_Intercept"),
  post_interval_pos = posterior_interval(bfit, "b_Intercept"),
  post_samples_pos = dim(posterior_samples(bfit, "b_Intercept")),
  pp_mixture_pos = tryCatch(pp_mixture(bfit, nd), error = conditionMessage),
  epred_reform_na = posterior_epred(bfit, nd, NA, ndraws = 20),
  predict_reform_na = posterior_predict(bfit, nd, NA, ndraws = 20),
  linpred_pos = posterior_linpred(bfit, TRUE, nd, ndraws = 20),
  loglik_pos = dim(log_lik(bfit, nd)),
  # named-argument forms
  post_interval_prob = posterior_interval(bfit, prob = 0.5),
  pred_interval_prob = predictive_interval(bfit, prob = 0.5, ndraws = 40),
  epred_newdata = posterior_epred(bfit, newdata = nd, ndraws = 20),
  predict_newdata = posterior_predict(bfit, newdata = nd, ndraws = 20),
  loglik_newdata = log_lik(bfit, newdata = nd),
  stancode_v = stancode(bfit, version = FALSE),
  standata_nd = names(standata(bfit, newdata = nd)),
  mcmc_plot_hist = {
    p <- mcmc_plot(bfit, type = "hist")
    lapply(ggplot2::ggplot_build(p)$data, function(d) dim(d))
  },
  rhat_pars = rhat(bfit, pars = "b_Intercept"),
  neff_pars = neff_ratio(bfit, pars = "b_Intercept"),
  # --- the print and summary surface --------------------------------
  print_fit = paste(capture.output(print(bfit)), collapse = "\n"),
  summary_fit = paste(capture.output(print(summary(bfit))),
                      collapse = "\n"),
  # --- names frmtmb.sample RE-EXPORTS from core, which the lane's
  #     list omitted although they are on the search path in every arm
  fixef = fixef(bfit),
  ranef_dim = lapply(ranef(bfit), dim),
  VarCorr = lapply(VarCorr(bfit), function(z) lapply(z, dim)),
  loo_fit = {
    l <- loo(bfit)
    list(l$estimates, names(l))
  },
  waic_fit = waic(bfit)$estimates,
  bayes_R2 = bayes_R2(bfit),
  hypothesis = hypothesis(bfit, "Intercept = 0")$hypothesis,
  prior_summary = as.data.frame(prior_summary(bfit)),
  posterior_summary = posterior_summary(bfit),
  ndraws = ndraws(bfit),
  nchains = nchains(bfit),
  niterations = niterations(bfit),
  nvariables = nvariables(bfit),
  variables = variables(bfit),
  as_draws_df_dim = dim(as_draws_df(bfit)),
  as_draws_rvars_n = length(as_draws_rvars(bfit)),
  conditional_effects = {
    ce <- conditional_effects(bfit)
    lapply(ce, function(d) list(dim(d), names(d)))
  },
  pp_check_n = {
    p <- pp_check(bfit, ndraws = 10)
    lapply(ggplot2::ggplot_build(p)$data, function(d) dim(d))
  }
)

hash <- function(x) {
  f <- tempfile(); on.exit(unlink(f))
  saveRDS(x, f, compress = FALSE, version = 3)
  unname(tools::md5sum(f))
}
out <- list()
for (nm in names(calls)) {
  set.seed(20260915)
  warns <- character(); msgs <- character()
  val <- withCallingHandlers(
    tryCatch(eval(calls[[nm]]),
             error = function(e) structure(conditionMessage(e),
                                           class = "sgrev_err")),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  kind <- if (inherits(val, "sgrev_err")) "ERROR" else "value"
  pr <- tryCatch(paste(capture.output(print(val)), collapse = "\n"),
                 error = function(e) "<unprintable>")
  out[[nm]] <- list(kind = kind, val = hash(val), warn = hash(warns),
                    msg = hash(msgs), printed = hash(pr),
                    warn_txt = warns, msg_txt = msgs,
                    err_txt = if (kind == "ERROR") unclass(val) else "",
                    cls = paste(class(val), collapse = ","))
}
saveRDS(out, av[3])
cat("OK", arm, length(out), "\n")
