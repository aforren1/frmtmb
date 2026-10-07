# Lane surface, item 9, second step: the same objective sampled from
# the mode init and from a random init, and the first warmup
# iterations' step sizes, to tell an init defect from a density defect.
#
#   Rscript dev/surface-gp-diag2.R lane|base > dev/surface-out/gp-diag2-<arm>.txt
arm <- commandArgs(TRUE)[1]
source("dev/surface-env.R")
surface_env(arm)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "| frmtmb", as.character(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb.sample")
set.seed(5)
n <- 60
d9 <- data.frame(x = round(stats::runif(n, 0, 6), 1))
d9$y <- 0.5 + sin(d9$x) + stats::rnorm(n, 0, 0.3)
fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = d9)
rp <- ns$sample_resolve_priors(fit, NULL, base = fit$prior, defaults = TRUE)
run <- function(tag, obj, init, ...) {
  sf <- suppressWarnings(tmbstan::tmbstan(obj, chains = 1, iter = 300,
                                          refresh = 0, seed = 4,
                                          init = init, ...))
  sp <- rstan::get_sampler_params(sf, inc_warmup = TRUE)[[1]]
  cat(sprintf("%-34s first stepsizes %s | post-warmup accept %.3f\n", tag,
              paste(format(sp[1:4, "stepsize__"], digits = 3),
                    collapse = " "),
              mean(sp[151:300, "accept_stat__"])))
}
obj <- ns$prior_augmented_obj(fit, rp$ri$entries)
lpb <- obj$env$last.par.best
init_mode <- ns$stan_init_arrays(ns$mode_inits(lpb, 1, 0.25))
run("mode init, prior obj", obj, init_mode)
obj2 <- ns$prior_augmented_obj(fit, rp$ri$entries)
run("random init, prior obj", obj2, "random")
obj3 <- ns$prior_augmented_obj(fit, rp$ri$entries)
# the same point, moved off the exact mode by a small jitter
set.seed(1)
init_j <- list(array(lpb + rnorm(length(lpb), 0, 1e-3), dim = length(lpb)))
run("mode + N(0, 1e-3) init", obj3, init_j)
obj4 <- ns$prior_augmented_obj(fit, rp$ri$entries)
# frm_sample's own call: does the last.par.best that frm_sample reads
# differ from the template's estimates?
cat("lpb vs unlist(fit$estimates) max abs diff:",
    format(max(abs(lpb - unlist(fit$estimates[names(
      fit$frame$par_template)])))), "\n")
cat("names of lpb:", paste(unique(names(lpb)), collapse = " "), "\n")
f0 <- obj4$env$f(lpb, order = 0)
cat("nll at lpb via obj4 before any fn call:", f0, "\n")
