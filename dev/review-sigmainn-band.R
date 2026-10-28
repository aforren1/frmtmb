# Short check: conditional_effects(method = "posterior_predict") band
# under ar(cov = TRUE), brms 2.23.0 against frmtmb (ML) and frmtmb.sample
# (draws), on the data of dev/review-sigmainn-readers.R (seed 20261028).
# Usage: Rscript dev/review-sigmainn-band.R brms|lane|r14|summary
arm <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sigmainn"
out <- file.path(wt, "dev", "review-sigmainn-out")
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-sigmainn-lib",
            if (arm != "brms") "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
set.seed(20261028)
ng <- 20; k <- 8
d <- data.frame(t = rep(seq_len(k), ng), g = factor(rep(seq_len(ng), each = k)))
d$x <- rnorm(nrow(d)); d$x2 <- rnorm(nrow(d))
d$y <- 1 + 0.5 * d$x + unlist(lapply(seq_len(ng), function(i)
  as.numeric(arima.sim(list(ar = 0.8), k, sd = 0.5))))
hw <- function(ce) (ce$upper__ - ce$lower__) / (2 * qnorm(0.975))
if (arm == "summary") {
  B <- readRDS(file.path(out, "band-brms.rds"))
  for (a in c("lane", "r14")) {
    L <- readRDS(file.path(out, paste0("band-", a, ".rds")))
    cat(sprintf("%-4s ML band/sigma(fit) %.4f (band %.4f, sigma %.4f) | draws band %.4f vs brms band %.4f: ratio %.4f (range of pointwise ratios %.3f to %.3f)\n",
                a, mean(L$ml_hw) / L$sigma, mean(L$ml_hw), L$sigma, mean(L$dr_hw),
                mean(B$hw), mean(L$dr_hw) / mean(B$hw), min(L$dr_hw / B$hw), max(L$dr_hw / B$hw)))
  }
  cat(sprintf("brms band %.4f, posterior mean sigma %.4f, band/sigma %.4f\n",
              mean(B$hw), B$sigma, mean(B$hw) / B$sigma))
  quit(save = "no")
}
if (arm == "brms") {
  suppressPackageStartupMessages(library(brms))
  bf1 <- brm(y ~ x + ar(time = t, gr = g, cov = TRUE), data = d, chains = 4,
             iter = 2000, seed = 1, refresh = 0, backend = "rstan")
  set.seed(2)
  ce <- conditional_effects(bf1, effects = "x", method = "posterior_predict")[[1]]
  saveRDS(list(hw = hw(ce), sigma = mean(as.matrix(bf1, variable = "sigma"))),
          file.path(out, "band-brms.rds"))
  quit(save = "no")
}
suppressPackageStartupMessages({ library(frmtmb); library(frmtmb.sample) })
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
fit <- frm(y ~ x + ar(time = t, gr = g, cov = TRUE), data = d)
set.seed(2)
ce <- conditional_effects(fit, effects = "x", method = "posterior_predict", ndraws = 4000)[[1]]
ds <- frm_sample(y ~ x + ar(time = t, gr = g, cov = TRUE), data = d, chains = 4,
                 iter = 2000, seed = 1, refresh = 0)
set.seed(2)
ced <- conditional_effects(ds, effects = "x", method = "posterior_predict")[[1]]
saveRDS(list(ml_hw = hw(ce), sigma = sigma(fit), dr_hw = hw(ced)),
        file.path(out, paste0("band-", arm, ".rds")))
