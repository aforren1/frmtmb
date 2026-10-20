# Lane sampleprior: what frmtmb's shifted_lognormal objective does at a
# non-decision time past the smallest response, with and without brms's
# uniform(0, min_Y) default (data seed 20261020).
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
set.seed(20261020)
d <- data.frame(x = stats::rnorm(60))
d$y <- 0.2 + exp(stats::rnorm(60, 0.3 * d$x))
fit <- frm(bf(y ~ x), family = shifted_lognormal(), data = d)
rp <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, NULL))
fit$prior <- rp$effective
obj <- frmtmb.sample:::prior_augmented_obj(fit, rp$ri$entries)
p0 <- unlist(fit$estimates[c("beta", "betad")])
j <- grep("ndt", names(p0))
stopifnot(length(j) == 1L)
for (ndt in c(0.5, 0.99, 1.01, 2) * min(d$y)) {
  p <- p0
  p[j] <- log(ndt)
  cat(sprintf("ndt = %.4f x min(y): likelihood %s, with the default %s\n",
              ndt / min(d$y), format(fit$obj$fn(p)), format(obj$fn(p))))
}
