# frm_sample(fit) on lane surface's exact y ~ gp(x): per chain, stepsize,
# divergences and means, under init = "last.par.best" and "random".
#   Rscript dev/moparity-gp-init.R
.libPaths(c("C:/Users/adf44/source/r/wt-moparity-lib", "C:/Users/adf44/source/r/rellib-r7", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(5); n <- 60
d9 <- data.frame(x = round(stats::runif(n, 0, 6), 1))
d9$y <- 0.5 + sin(d9$x) + stats::rnorm(n, 0, 0.3)
fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = d9)
cat("ML theta", fit$estimates$theta, "\n")
for (init in c("last.par.best", "random")) {
s <- suppressWarnings(suppressMessages(frm_sample(fit, chains = 4, iter = 1000, seed = 1, refresh = 0, cores = 1, init = init)))
sp <- rstan::get_sampler_params(s$stanfit, inc_warmup = FALSE)
M <- s$draws
nch <- 4; per <- nrow(M) / nch
for (k in 1:4) {
  rows <- (k - 1) * per + seq_len(per)
  cat(sprintf("%s chain %d: stepsize %.4g div %d mean sigma %.4f mean sdgp %.3f\n", init, k, sp[[k]][1, "stepsize__"], sum(sp[[k]][, "divergent__"]), mean(M[rows, "sigma"]), mean(M[rows, "sdgp_gpx"])))
}
}
