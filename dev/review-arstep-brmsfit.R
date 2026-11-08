# Review of lane arstep: brms 2.23.0 on the lane's step row beside ar()
# (step n 80 seed 22, data as tests/testthat/test-gp-ar-start.R ga_step()),
# where frmtmb's maximum has ar 1.022. brms's own priors, sampled with rstan
# from brms's Stan code and data (no brm() compile cache: the program is
# compiled here or read from dev/stan-cache by its content key, as
# helper-brms.R does, and a new program is saved under review-arstep-out).
# Does brms's posterior sit past the unit root too?
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/review-arstep-lib.R")
rv_libs("lane")
suppressPackageStartupMessages({library(brms); library(rstan)})
cat("brms", as.character(packageVersion("brms")), "rstan", as.character(packageVersion("rstan")), "\n")
n <- 80L; seed <- 22L
set.seed(74000 + seed * 31 + n)
x <- sort(round(runif(n, 0, 10), 4))
d <- data.frame(x = x, t = seq_len(n),
                y = 2 * (x > 5) - 1 + as.numeric(stats::arima.sim(list(ar = 0.4), n, sd = 0.03)))
f <- y ~ gp(x) + ar(time = t)
fit <- frm(f, data = d, family = gaussian())
cat(sprintf("FRM ll %.4f ar %.4f sdgp %.4f lscale %.4f sigma %.4f\n", as.numeric(logLik(fit)),
            fit$estimates$thetaac, exp(fit$estimates$theta[1]), exp(fit$estimates$theta[2]),
            exp(fit$estimates$betad[1])))
code <- brms::stancode(f, data = d)
sdata <- brms::standata(f, data = d)
key <- {
  tf <- tempfile(fileext = ".stan")
  writeLines(c(code, paste0("// rstan ", packageVersion("rstan"))), tf)
  unname(tools::md5sum(tf))
}
p1 <- file.path(RV_WT, "dev/stan-cache", paste0(key, ".rds"))
mydir <- file.path(RV_OUT, "stan-cache"); dir.create(mydir, showWarnings = FALSE)
p2 <- file.path(mydir, paste0(key, ".rds"))
mod <- if (file.exists(p1)) { cat("cache: worktree\n"); readRDS(p1) } else if (file.exists(p2)) {
  cat("cache: review\n"); readRDS(p2)
} else {
  cat("cache: compiling\n")
  m <- rstan::stan_model(model_code = code, save_dso = TRUE); saveRDS(m, p2); m
}
t0 <- proc.time()[[3]]
sf <- rstan::sampling(mod, data = sdata, chains = 2, iter = 1500, seed = 4401,
                      refresh = 0, cores = 1)
cat("sampling secs", proc.time()[[3]] - t0, "\n")
dr <- as.matrix(sf, pars = c("ar", "sdgp_1", "lscale_1", "sigma"))
qs <- apply(dr, 2, quantile, c(0.05, 0.5, 0.95))
print(round(qs, 4))
cat(sprintf("BRMS P(ar > 1) %.3f  P(ar > 0.9) %.3f  draws %d  divergent %d\n",
            mean(dr[, 1] > 1), mean(dr[, 1] > 0.9), nrow(dr),
            sum(rstan::get_divergent_iterations(sf))))
print(rstan::summary(sf, pars = c("ar", "sdgp_1", "lscale_1", "sigma"))$summary[, c("mean", "n_eff", "Rhat")])
cat("DONE\n")
