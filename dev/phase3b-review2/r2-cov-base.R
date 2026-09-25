# Reviewer 2, item 6: Wald coverage of the hierarchical drift intercept
# on the BASE build (rellib-r3, frmtmb 0.62.0 / frmtmb.eam 0.10.0), a
# plain wiener() with no censoring and no contaminant, at the worker's
# design (dev/phase3b-eam-recovery7.R, simulate_one() with no arm
# edits; fit_one()'s "cens" formula without cens()).
# Usage: Rscript r2-cov-base.R <lib: base|lane> <seed> [<seed> ...]
# Writes dev/phase3b-review2/cov-<lib>/seed-<seed>.rds
source("dev/phase3b-review2/r2-prelude.R")
a <- commandArgs(trailingOnly = TRUE)
which_lib <- a[[1]]
r2_lib(which_lib)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
seeds <- as.integer(a[-1])
out_dir <- file.path("dev/phase3b-review2", paste0("cov-", which_lib))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
NS <- 30L; NT <- 400L

for (seed in seeds) {
  t0 <- proc.time()[["elapsed"]]
  set.seed(seed)
  u <- rnorm(NS, 0, 0.35)
  b <- rnorm(NS, 0, 0.20)
  d <- do.call(rbind, lapply(seq_len(NS), function(s) {
    cond <- rep(0:1, length.out = NT)
    x <- ddm_simulate(NT, mu = 0.4 + 0.9 * cond + u[s],
                      bs = 1.4 * exp(b[s]), ndt = 0.25, bias = 0.5)
    x$cond <- cond
    x$s <- factor(s, levels = seq_len(NS))
    x
  }))
  res <- tryCatch({
    fit <- frm(bf(rt | dec(upper) ~ cond + (1 | s),
                  bs ~ 1 + (1 | s), ndt ~ 1, bias = 0.5),
               family = wiener(), data = d)
    ci <- confint(fit)
    V <- vcov(fit)
    list(seed = seed, ci = ci, vn = dimnames(V),
         V = V, conv = fit$opt$convergence, loglik = as.numeric(logLik(fit)),
         mean_u = mean(u), sd_u = sd(u), mean_b = mean(b),
         pdHess = tryCatch(fit$sdr$pdHess, error = function(e) NA),
         elapsed = proc.time()[["elapsed"]] - t0,
         eam_path = find.package("frmtmb.eam"),
         eam_ver = as.character(packageVersion("frmtmb.eam")))
  }, error = function(e) list(seed = seed, error = conditionMessage(e)))
  saveRDS(res, file.path(out_dir, sprintf("seed-%d.rds", seed)))
  cat(which_lib, seed, if (is.null(res$error)) "ok" else res$error,
      round(proc.time()[["elapsed"]] - t0, 1), "\n")
}
