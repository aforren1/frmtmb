# Lane surface, item 9: frm_sample(fit) on the gpby review's exact
# y ~ gp(x) fit (60 points, data seed 5), over sampler seeds 1 to 8,
# chains = 1, iter = 600. Does the chain move, on one build?
#
#   Rscript dev/surface-gp-sweep.R lane|base > dev/surface-out/gp-sweep-<arm>.txt
arm <- commandArgs(TRUE)[1]
source("dev/surface-env.R")
surface_env(arm)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "| frmtmb", find.package("frmtmb"), "\n")
options(mc.cores = 1)
set.seed(5)
n <- 60
d9 <- data.frame(x = round(stats::runif(n, 0, 6), 1))
d9$y <- 0.5 + sin(d9$x) + stats::rnorm(n, 0, 0.3)
fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = d9)
for (s in 1:8) {
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 600, refresh = 0, seed = s)))
  sp <- rstan::get_sampler_params(ds$stanfit, inc_warmup = FALSE)[[1]]
  m <- as.matrix(ds)
  cat(sprintf(paste0("seed %d: accept %.3f, stepsize %.3g, divergent %d,",
                     " distinct draws of b_Intercept %d\n"),
              s, mean(sp[, "accept_stat__"]), sp[1, "stepsize__"],
              sum(sp[, "divergent__"]), length(unique(m[, "b_Intercept"]))))
}
