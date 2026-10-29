# Review of lane meme: the sampled densities on sdme and corme, read
# off draws. sample_prior = "only" keeps the measurement density of the
# noisy values (as brms does), so the noise is made 100 times the
# latent sd: what x says about the hyperparameters is then negligible
# and their draws are their priors. Three ungrouped me() terms give one
# 3 x 3 correlation, whose every pair is 2 Beta(eta + 1/2, eta + 1/2) - 1
# under lkj(eta). sdme ~ lognormal(0, 0.3), so log(sdme) ~ N(0, 0.3).
# Seeds: data 5, sampler 11. Output dev/review-meme-out/prioronly.txt
.libPaths(c("C:/Users/adf44/source/r/wt-meme-lib",
            "C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
set.seed(5)
n <- 3
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n), s = 100,
                y = rnorm(n))
eta <- 2
pr <- set_prior("normal(0, 1)", class = "b") +
  set_prior("normal(0, 1)", class = "Intercept") +
  set_prior("exponential(1)", class = "sigma") +
  set_prior("normal(0, 1)", class = "meanme") +
  set_prior("lognormal(0, 0.3)", class = "sdme") +
  set_prior(paste0("lkj(", eta, ")"), class = "corme")
f <- bf(y ~ me(x, s) + me(z, s) + me(w, s))
zt <- function(v, mean_t, sd_t, q_t, lab) {
  m <- matrix(v, ncol = 4)
  q <- stats::quantile(v, c(0.05, 0.5, 0.95), names = FALSE)
  z <- c(mean = (mean(v) - mean_t) / posterior::mcse_mean(m),
         sd = (stats::sd(v) - sd_t) / posterior::mcse_sd(m),
         q05 = (q[1] - q_t[1]) / posterior::mcse_quantile(m, 0.05),
         q50 = (q[2] - q_t[2]) / posterior::mcse_quantile(m, 0.5),
         q95 = (q[3] - q_t[3]) / posterior::mcse_quantile(m, 0.95))
  cat(sprintf("%-22s ess %6.0f  z: %s\n", lab, posterior::ess_bulk(m),
              paste(sprintf("%s %6.2f", names(z), z), collapse = "  ")))
  z
}
a <- eta + 1 / 2
beta_mom <- function(a) {
  c(mean = 0, sd = 2 * sqrt(a * a / ((2 * a)^2 * (2 * a + 1))))
}
bq <- 2 * stats::qbeta(c(0.05, 0.5, 0.95), a, a) - 1
lq <- stats::qnorm(c(0.05, 0.5, 0.95), 0, 0.3)
allz <- c()
for (sp in c("only", "yes")) {
  cat("==== sample_prior =", sp, "\n")
  t0 <- proc.time()[["elapsed"]]
  ds <- suppressWarnings(suppressMessages(frm_sample(
    f, data = d, family = gaussian(), prior = pr, sample_prior = sp,
    chains = 4, iter = 3000, warmup = 1000, seed = 11, refresh = 0,
    control = list(adapt_delta = 0.995))))
  cat("secs", round(proc.time()[["elapsed"]] - t0), "\n")
  m <- as.matrix(ds)
  cat("divergent:", tryCatch(sum(rstan::get_divergent_iterations(ds$stanfit)),
                             error = function(e) NA), "\n")
  if (sp == "only") {
    for (cn in grep("^corme__", colnames(m), value = TRUE)) {
      allz <- c(allz, zt(m[, cn], 0, beta_mom(a)[["sd"]], bq, cn))
    }
    for (cn in grep("^sdme_", colnames(m), value = TRUE)) {
      allz <- c(allz, zt(log(m[, cn]), 0, 0.3, lq, paste0("log ", cn)))
    }
    for (cn in grep("^meanme_", colnames(m), value = TRUE)) {
      allz <- c(allz, zt(m[, cn], 0, 1, stats::qnorm(c(0.05, 0.5, 0.95)), cn))
    }
  } else {
    pd <- prior_draws(ds)
    cat("prior_draws:", names(pd), "\n")
    # independent draws, so a KS test is valid here
    cat(sprintf("KS prior_Lme_1 vs 2 Beta(%.1f, %.1f) - 1: p = %.3g\n", a, a,
                suppressWarnings(stats::ks.test((pd$Lme_1 + 1) / 2, "pbeta",
                                                a, a)$p.value)))
    cat(sprintf("KS prior_sdme_1 vs lognormal(0, 0.3): p = %.3g\n",
                suppressWarnings(stats::ks.test(pd$sdme_1, "plnorm", 0,
                                                0.3)$p.value)))
    cat(sprintf("KS prior_meanme_1 vs N(0, 1): p = %.3g\n",
                suppressWarnings(stats::ks.test(pd$meanme_1, "pnorm")$p.value)))
  }
}
cat(sprintf("RESULT max |z| %.2f over %d statistics\n", max(abs(allz)),
            length(allz)))
cat("DONE\n")
