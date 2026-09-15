# The side the lane did not measure. Its output harness held brms to
# its own output on a brmsfit. The fix also changes WHICH GENERIC RUNS
# for a frmtmb_draws, because the owner's generic is now the one
# reached, and a generic can warn, deprecate or coerce before it
# dispatches. So: the same calls on the same draws object, BASE
# against FIX, in each load order, on value, warnings, messages and
# the printed form.
#
#   Rscript dev/sgrev-drawsout.R <LIB> <arm> <outfile>
av <- commandArgs(trailingOnly = TRUE)
.libPaths(c(av[1], "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
arm <- av[2]
q <- function(e) suppressWarnings(suppressMessages(e))
if (arm == "N") q(library(frmtmb.sample))
if (arm == "S") { q(library(brms)); q(library(frmtmb.sample)) }
if (arm == "T") { q(library(frmtmb.sample)); q(library(brms)) }
if (arm == "P") {
  for (p in c("bayesplot", "bridgesampling", "coda", "gratia", "loo",
              "posterior", "rstantools")) q(library(p, character.only = TRUE))
  q(library(frmtmb.sample))
}
if (arm == "Q") {
  q(library(frmtmb.sample))
  for (p in c("bayesplot", "bridgesampling", "coda", "gratia", "loo",
              "posterior", "rstantools")) q(library(p, character.only = TRUE))
}
ds <- readRDS("dev/stan-cache/sgrev-draws.rds")
nd <- ds$fit$frame$data[1:5, , drop = FALSE]

calls <- alist(
  as.mcmc = dim(as.matrix(as.mcmc(ds))),
  bayes_factor = bayes_factor(ds, ds),
  bridge_sampler = bridge_sampler(ds),
  post_prob = post_prob(ds),
  kfold = kfold(ds),
  loo_moment_match = loo_moment_match(ds),
  loo_subsample = loo_subsample(ds),
  psis = { p <- psis(ds); list(class(p), dim(p$log_weights)) },
  log_lik = log_lik(ds, ndraws = 20),
  nsamples = nsamples(ds),
  posterior_epred = posterior_epred(ds, ndraws = 20),
  posterior_interval = posterior_interval(ds),
  posterior_linpred = posterior_linpred(ds, ndraws = 20),
  posterior_predict = posterior_predict(ds, ndraws = 20),
  predictive_error = predictive_error(ds, ndraws = 20),
  predictive_interval = predictive_interval(ds, ndraws = 20),
  log_posterior = log_posterior(ds),
  neff_ratio = neff_ratio(ds),
  nuts_params = nuts_params(ds),
  rhat = rhat(ds),
  mcmc_plot = lapply(ggplot2::ggplot_build(mcmc_plot(ds))$data, dim),
  parnames = parnames(ds),
  posterior_samples = posterior_samples(ds),
  pp_mixture = pp_mixture(ds),
  reloo = reloo(ds),
  restructure = restructure(ds),
  stancode = stancode(ds),
  standata = standata(ds),
  # core names frmtmb.sample re-exports
  loo = { l <- loo(ds); list(l$estimates, names(l)) },
  waic = waic(ds)$estimates,
  fixef = fixef(ds),
  ranef = lapply(ranef(ds), dim),
  VarCorr = VarCorr(ds),
  hypothesis = hypothesis(ds, "x = 0")$hypothesis,
  bayes_R2 = bayes_R2(ds),
  pp_check = lapply(ggplot2::ggplot_build(pp_check(ds, ndraws = 10))$data,
                    dim),
  conditional_effects = lapply(conditional_effects(ds),
                               function(d) list(dim(d), names(d))),
  posterior_summary = posterior_summary(ds),
  ndraws = ndraws(ds), nchains = nchains(ds),
  niterations = niterations(ds), nvariables = nvariables(ds),
  variables = variables(ds),
  as_draws_df = dim(as_draws_df(ds)),
  summary_ds = paste(capture.output(print(summary(ds))), collapse = "\n"),
  print_ds = paste(capture.output(print(ds)), collapse = "\n")
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
    warning = function(w) { warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning") },
    message = function(m) { msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage") })
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
