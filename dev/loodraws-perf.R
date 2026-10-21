# Time and memory of loo_predict() at 4000 draws x 1e4 rows, against
# brms 2.23.0. Arm A runs both PSIS-plus-E_loo paths on the SAME
# matrices (frmtmb's log_lik() and posterior_predict() of one draws
# object), interleaved, with a control that repeats frmtmb's path.
# Arm B runs each package's loo_predict() end to end on its own fit of
# the same data. Seeds: data 1, samplers 1, predictions 3.
#   Rscript dev/loodraws-perf.R  -> dev/loodraws-log/perf.txt
.libPaths(c("C:/Users/adf44/source/r/wt-loodraws-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-loodraws/dev/stan-cache")
options(mc.cores = 1)
suppressMessages({library(frmtmb); library(frmtmb.sample); library(brms)})
cat("lib:", dirname(system.file(package = "frmtmb.sample")),
    format(packageVersion("frmtmb.sample")), "brms",
    format(packageVersion("brms")), "loo", format(packageVersion("loo")),
    "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(1)
n <- 1e4
dd <- data.frame(x = rnorm(n))
dd$y <- rnorm(n, 1 + 0.5 * dd$x)
t0 <- proc.time()[["elapsed"]]
ds <- q(frm_sample(frmtmb::bf(y ~ x), family = gaussian(), data = dd, chains = 4,
                   iter = 2000, warmup = 1000, seed = 1, refresh = 0))
cat("frm_sample s:", round(proc.time()[["elapsed"]] - t0, 1), "draws",
    ndraws(ds), "\n")
t0 <- proc.time()[["elapsed"]]
bf_fit <- q(brm(y ~ x, data = dd, chains = 4, iter = 2000, seed = 1,
                refresh = 0, cores = 1))
cat("brm s:", round(proc.time()[["elapsed"]] - t0, 1), "draws",
    ndraws(bf_fit), "\n")

# memory: the peak of R's heap during the call, above what was live
# before it (gc()'s "max used" after a reset, in MB)
measure <- function(f) {
  invisible(gc(reset = TRUE))
  base <- sum(gc()[, 2])
  t <- system.time(v <- f())[["elapsed"]]
  peak <- sum(gc()[, 6])
  list(t = t, mb = peak - base, v = v)
}

t0 <- proc.time()[["elapsed"]]
ll <- log_lik(ds)
t_ll <- proc.time()[["elapsed"]] - t0
set.seed(3)
t0 <- proc.time()[["elapsed"]]
yrep <- posterior_predict(ds)
t_pp <- proc.time()[["elapsed"]] - t0
cat(sprintf("frmtmb log_lik() %.1f s, posterior_predict() %.1f s, %d x %d\n",
            t_ll, t_pp, nrow(ll), ncol(ll)))
fr <- function() {
  ps <- frmtmb.sample:::draws_psis_of(ll, ds)
  frmtmb.sample:::loo_expect_value(yrep, ps, "mean", 0.5)
}
br <- function() {
  # brms's loo_predict() body after its log_lik() and posterior_predict():
  # loo(save_psis = TRUE) with r_eff_log_lik(), then E_loo_value()
  r_eff <- loo::relative_eff(exp(ll), chain_id = rep(1:4, each = 1000))
  lo <- q(loo::loo(ll, r_eff = r_eff, save_psis = TRUE))
  brms:::E_loo_value(yrep, lo$psis_object, type = "mean", probs = 0.5)
}
arms <- list(frmtmb = fr, brms = br, control = fr)
res <- list()
for (round in 1:3) for (a in names(arms)) {
  m <- measure(function() q(arms[[a]]()))
  res[[length(res) + 1L]] <- data.frame(round = round, arm = a, s = m$t,
                                        mb = m$mb)
  if (round == 1 && a == "frmtmb") vf <- m$v
  if (round == 1 && a == "brms") vb <- m$v
}
res <- do.call(rbind, res)
print(res)
cat("ARM A (same matrices) min s:",
    paste(sprintf("%s %.2f", unique(res$arm),
                  tapply(res$s, res$arm, min)[unique(res$arm)]),
          collapse = ", "), "\n")
cat("ARM A peak MB (max over rounds):",
    paste(sprintf("%s %.0f", unique(res$arm),
                  tapply(res$mb, res$arm, max)[unique(res$arm)]),
          collapse = ", "), "\n")
cat("ARM A frmtmb == brms E_loo_value on the same matrices:",
    identical(vf, vb), "\n")

# arm B: each package end to end, one round each, fresh gc
eb <- list()
for (round in 1:2) {
  set.seed(3)
  mf <- measure(function() q(loo_predict(ds)))
  set.seed(3)
  mb <- measure(function() q(loo_predict(bf_fit)))
  eb[[round]] <- data.frame(round = round, arm = c("frmtmb", "brms"),
                            s = c(mf$t, mb$t), mb = c(mf$mb, mb$mb))
}
eb <- do.call(rbind, eb)
print(eb)
cat("ARM B (end to end) min s:",
    sprintf("frmtmb %.1f, brms %.1f", min(eb$s[eb$arm == "frmtmb"]),
            min(eb$s[eb$arm == "brms"])), "\n")
cat("ARM B peak MB:",
    sprintf("frmtmb %.0f, brms %.0f", max(eb$mb[eb$arm == "frmtmb"]),
            max(eb$mb[eb$arm == "brms"])), "\n")
