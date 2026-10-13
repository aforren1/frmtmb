# Reviewer: the default frm_sample(fit) call (init "last.par.best": chain
# 1 at the ML mode, chains 2-4 jittered) on an exact y ~ gp(x), the
# construction of dev/moparity-gp-init.R, over data seeds, with the
# formula route (init "random") beside it.
#   Rscript dev/moparity-rev-gp-init.R base|lane dataseed
args <- commandArgs(TRUE)
arm <- args[1]
ds <- as.integer(args[2])
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
set.seed(ds); n <- 60
d9 <- data.frame(x = round(stats::runif(n, 0, 6), 1))
d9$y <- 0.5 + sin(d9$x) + stats::rnorm(n, 0, 0.3)
fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = d9)
cat("data seed", ds, "ML theta", format(fit$estimates$theta, digits = 4),
    "sigma", format(exp(fit$estimates$betad[1]), digits = 4), "\n")
report <- function(lab, s) {
  sp <- rstan::get_sampler_params(s$stanfit, inc_warmup = FALSE)
  M <- s$draws
  nch <- length(sp); per <- nrow(M) / nch
  cols <- intersect(c("b_Intercept", "sigma", "sdgp_gpx", "lscale_gpx"),
                    colnames(M))
  arr <- array(M[, cols], c(per, nch, length(cols)),
               dimnames = list(NULL, NULL, cols))
  sm <- posterior::summarise_draws(posterior::as_draws_array(arr), "rhat",
                                   "ess_bulk")
  div <- vapply(sp, function(x) sum(x[, "divergent__"]), 0)
  td <- vapply(sp, function(x) sum(x[, "treedepth__"] >= 10), 0)
  cat(sprintf("RES %s %s data seed %d | div per chain %s | max treedepth hits %d | max R-hat %.3f | min bulk ESS %.0f | chain means sigma %s\n",
              arm, lab, ds, paste(div, collapse = ","), sum(td),
              max(sm$rhat), min(sm$ess_bulk),
              paste(vapply(seq_len(nch), function(k) {
                format(mean(M[(k - 1) * per + seq_len(per), "sigma"]),
                       digits = 3)
              }, ""), collapse = ",")))
}
t0 <- proc.time()[[3]]
s1 <- suppressWarnings(suppressMessages(frm_sample(
  fit, chains = 4, iter = 1000, seed = 1, refresh = 0, cores = 1)))
report("fit-default", s1)
s2 <- suppressWarnings(suppressMessages(frm_sample(
  bf(y ~ gp(x)), family = gaussian(), data = d9, chains = 4, iter = 1000,
  seed = 1, refresh = 0, cores = 1)))
report("formula-default", s2)
cat("secs", round(proc.time()[[3]] - t0), "\n")
