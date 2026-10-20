# Lane sampleprior, item 4: does frm_sample()'s init jitter need the
# fit's covariance factor? The jitter today is isotropic, sd init_jitter
# = 0.25 on the sampler's coordinates (mode_inits()); it is NOT scaled
# by the standard errors at all. Arms, on illse's ill-conditioned
# designs (dev/illse-designs.R, no random effects):
#   iso     the default: mode + N(0, 0.25^2 I)
#   factor  mode + F z, z ~ N(0, I), F F' = V the fit's covariance held
#           as a factor (frm_joint_cov(): V = Vrest + root root'), so each
#           chain starts one standard error out along a random direction
#   factor2 mode + 2 F z
# Chain 1 starts at the mode in every arm. 4 chains, iter 1000, sampler
# seed = data seed. For each chain: the start's distance from the mode
# in standard errors (sqrt(d' V^-1 d), through F) and its log-density
# drop; for the run: divergences, max-treedepth hits, the largest R-hat
# and smallest bulk ESS over the outer parameters, and the elapsed time.
#   Rscript dev/sampleprior-jitter.R <design> <seed>
# design: poly7 | scales | collin6 | poly3 | poly4 | collin3 | offset5
a <- commandArgs(TRUE)
des <- a[1]
seed <- as.integer(a[2])
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sampleprior"
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
source(file.path(wt, "dev/illse-designs.R"))
mk <- switch(des, poly7 = des_poly(7), scales = des_scales(),
             poly3 = des_poly(3), poly4 = des_poly(4),
             collin3 = des_collin(1e-3), offset5 = des_offset(1e5),
             collin6 = des_collin(1e-6))
cs <- mk(seed)
fit <- suppressWarnings(frm(cs$formula, family = cs$family, data = cs$data))
jc <- suppressWarnings(frm_joint_cov(fit))
V <- jc$V
F <- if (!is.null(jc$root)) {
  Vr <- jc$Vrest
  if (is.null(Vr)) Vr <- matrix(0, nrow(V), nrow(V))
  e <- eigen((Vr + t(Vr)) / 2, symmetric = TRUE)
  keep <- e$values > 0
  cbind(e$vectors[, keep, drop = FALSE] %*%
          diag(sqrt(e$values[keep]), sum(keep)), jc$root)
} else {
  e <- eigen((V + t(V)) / 2, symmetric = TRUE)
  e$vectors %*% diag(sqrt(pmax(e$values, 0)))
}
cat(sprintf("%s seed %d: %d outer parameters, factor %d columns, root %s, ",
            cs$lab, seed, nrow(V), ncol(F), !is.null(jc$root)),
    sprintf("se range %.3g to %.3g, kappa(V) %.3g\n",
            min(sqrt(diag(V))), max(sqrt(diag(V))), kappa(V, exact = TRUE)))
mode <- unlist(fit$estimates[c("beta", "betad")])
stopifnot(length(mode) == nrow(V))
# distance of a start from the mode in standard errors, through F:
# d = F w has the least-norm w = F^+ d, and |w| is the Mahalanobis length
Fp <- MASS::ginv(F)
se_dist <- function(d) sqrt(sum((Fp %*% d)^2))
inits_of <- function(arm) {
  set.seed(seed)
  lapply(1:4, function(i) {
    if (i == 1L) return(mode)
    switch(arm,
      iso = mode + stats::rnorm(length(mode), 0, 0.25),
      factor = mode + as.vector(F %*% stats::rnorm(ncol(F))),
      factor2 = mode + 2 * as.vector(F %*% stats::rnorm(ncol(F))))
  })
}
for (arm in c("iso", "factor", "factor2")) {
  ini <- inits_of(arm)
  dist <- vapply(ini[-1], function(v) se_dist(v - mode), 1)
  t0 <- proc.time()[["elapsed"]]
  ds <- suppressWarnings(suppressMessages(frm_sample(
    fit, chains = 4, iter = 1000, seed = seed, refresh = 0,
    init = lapply(ini, unname))))
  el <- proc.time()[["elapsed"]] - t0
  sf <- ds$stanfit
  # the start's log density relative to the mode, on the sampled tape
  obj <- frmtmb.sample:::prior_augmented_obj(ds$fit,
    frmtmb.sample:::sample_resolve_priors(fit, NULL)$ri$entries)
  drop <- vapply(ini[-1], function(v) obj$fn(v) - obj$fn(mode), 1)
  arr <- posterior::as_draws_array(frmtmb.sample:::draws_raw_array(ds))
  outer <- setdiff(posterior::variables(arr), "lp__")
  sm <- posterior::summarise_draws(posterior::subset_draws(arr, outer),
                                   "rhat", "ess_bulk")
  cat(sprintf(paste0("  %-7s start se-dist %s | -logdens rise %s | div %d ",
                     "treedepth %d | max rhat %.3f min ess %.0f | %.1f s\n"),
              arm, paste(sprintf("%.3g", dist), collapse = " "),
              paste(sprintf("%.3g", drop), collapse = " "),
              rstan::get_num_divergent(sf), rstan::get_num_max_treedepth(sf),
              max(sm$rhat), min(sm$ess_bulk), el))
}
