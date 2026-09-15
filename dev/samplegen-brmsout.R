# brms held to OUTPUT, not only to reachability: the same calls on the
# same cached brmsfit, in a process with brms alone and in processes
# with frmtmb.sample loaded in each order. Each call's value, warnings
# and error are hashed; an arm agrees with the control when every hash
# matches.
#
# The fit is the 400-row lognormal brms fit dev/generics-scale-brms.R
# sampled with seed 2026 (2 chains x 2000), copied into this worktree's
# gitignored dev/stan-cache. Calls that would refit or recompile Stan
# (kfold, reloo, loo_moment_match, bridge_sampler, bayes_factor,
# post_prob) are left to the method-table probe, dev/samplegen-check.R.
#
#   Rscript dev/samplegen-brmsout.R <LIB> <arm>
#     none  library(brms)                          the control
#     S     library(brms); library(frmtmb.sample)
#     T     library(frmtmb.sample); library(brms)
#     U     library(frmtmb.sample); brms only loaded; calls through the
#           names frmtmb.sample exports
av <- commandArgs(trailingOnly = TRUE)
.libPaths(c(av[1], "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
arm <- av[2]
q <- function(e) suppressWarnings(suppressMessages(e))
if (arm == "none") q(library(brms))
if (arm == "S") { q(library(brms)); q(library(frmtmb.sample)) }
if (arm == "T") { q(library(frmtmb.sample)); q(library(brms)) }
if (arm == "U") { q(library(frmtmb.sample)); q(loadNamespace("brms")) }
bfit <- readRDS("dev/stan-cache/samplegen-brmsfit.rds")

calls <- alist(
  as.mcmc = as.mcmc(bfit),
  log_lik = log_lik(bfit, ndraws = 40),
  log_posterior = log_posterior(bfit),
  loo_subsample = loo_subsample(bfit, observations = 50),
  mcmc_plot = nrow(mcmc_plot(bfit)$data),
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
  restructure = restructure(bfit),
  rhat = rhat(bfit),
  stancode = stancode(bfit),
  standata = standata(bfit)
)

hash <- function(x) {
  f <- tempfile()
  on.exit(unlink(f))
  saveRDS(x, f, compress = FALSE, version = 3)
  unname(tools::md5sum(f))
}
for (nm in names(calls)) {
  set.seed(20260915)
  warns <- character()
  val <- withCallingHandlers(
    tryCatch(eval(calls[[nm]]),
             error = function(e) structure(conditionMessage(e),
                                           class = "err")),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  kind <- if (inherits(val, "err")) "ERROR" else "value"
  # restructure() returns the fit, which carries environments that no
  # two processes serialize alike; its class and draws are what it is
  if (nm == "restructure" && kind == "value") {
    val <- list(class(val), dim(as.matrix(val)))
  }
  cat(sprintf("%-20s %-5s %s %s  %s\n", nm, kind, hash(val),
              hash(warns),
              if (kind == "ERROR") substr(unclass(val), 1, 60) else ""))
}
cat("DONE\n")
