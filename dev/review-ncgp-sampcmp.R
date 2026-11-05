# Item 9: frm_sample() against frm() at one point, and the draws of this build against
# gpgap's route and brms's program at true sigma 1e-2. n 20, seed 6620. Args: arm(lane|gpgap|brms) chain
args <- commandArgs(TRUE)
arm <- args[1]; chain <- as.integer(args[2])
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
if (arm == "gpgap") {
  .libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib", REV_BASE, REV_USER))
  Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win", FRMTMB_STAN_CACHE = file.path(REV_WT, "dev/review-ncgp-stan-cache"))
  suppressPackageStartupMessages(library(frmtmb))
  cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
} else {
  rev_libs("lane")
  Sys.setenv(FRMTMB_STAN_CACHE = file.path(REV_WT, "dev/review-ncgp-stan-cache"))
}
set.seed(6620)
n <- 20L
x <- sort(runif(n, 0, 6)); y <- sin(2 * x) + 1e-2 * rnorm(n)
d <- data.frame(y = y, x = x)
if (arm == "brms") {
  suppressPackageStartupMessages(library(brms))
  b <- brm(y ~ gp(x), d, chains = 1, iter = 2000, seed = 100 + chain, refresh = 0, silent = 2, backend = "rstan", chain_id = chain,
           control = list(adapt_delta = 0.99, max_treedepth = 12))
  m <- as.matrix(b, variable = c("sigma", "sdgp_gpx", "lscale_gpx", "b_Intercept"))
  saveRDS(m, file.path(REV_OUT, sprintf("sampcmp-%s-%d.rds", arm, chain)))
  cat("RSM done", arm, chain, "\n")
} else {
  suppressPackageStartupMessages(library(frmtmb.sample))
  fit <- suppressWarnings(frm(y ~ gp(x), d))
  e <- fit$obj$env
  rp <- frmtmb.sample:::sample_resolve_priors(fit, "flat", defaults = FALSE)
  pa <- frmtmb.sample:::prior_augmented_obj(fit, rp$entries %||% list())
  e2 <- pa$env
  p1 <- e$last.par.best
  v1 <- e$f(p1)
  v2 <- tryCatch(e2$f(if (identical(names(e2$par), names(p1))) p1 else stop("order")), error = function(err) { cat("order differs", conditionMessage(err), "\n"); NA })
  cat(sprintf("RSM %s joint at frm's max: frm tape %.12f sampler tape %.12f diff %.2e; tape order %s\n", arm, v1, v2, v1 - v2, paste(unique(names(e2$par)), collapse = ",")))
  t0 <- Sys.time()
  dr <- suppressWarnings(suppressMessages(frm_sample(fit, chains = 1, iter = 2000, seed = 100 + chain, refresh = 0, control = list(adapt_delta = 0.99, max_treedepth = 12))))
  cat(sprintf("RSM %s chain %d wall %.0f s\n", arm, chain, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  saveRDS(dr$draws[, intersect(c("sigma", "sdgp_gpx", "lscale_gpx", "b_Intercept"), colnames(dr$draws))], file.path(REV_OUT, sprintf("sampcmp-%s-%d.rds", arm, chain)))
}
cat("DONE\n")
