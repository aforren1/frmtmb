# posterior_average()'s stacking weights before and after: loo objects
# (the 0.20.0 route) against log_lik() matrices with r_eff (the new
# route), on the same draws. Log: dev/loodraws-log/probe-pa.txt
.libPaths(c("C:/Users/adf44/source/r/wt-loodraws-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-loodraws/dev/stan-cache")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(1)
dd <- data.frame(x = rnorm(60))
dd$y <- rnorm(60, 1 + 0.5 * dd$x)
d1 <- q(frm_sample(bf(y ~ x), family = gaussian(), data = dd, chains = 2,
                   iter = 600, seed = 1, refresh = 0))
d2 <- q(frm_sample(bf(y ~ 1), family = gaussian(), data = dd, chains = 2,
                   iter = 600, seed = 2, refresh = 0))
for (m in c("stacking", "pseudobma")) {
  set.seed(5)
  old <- as.numeric(q(loo::loo_model_weights(list(loo(d1), loo(d2)),
                                             method = m)))
  set.seed(5)
  new <- as.numeric(q(frmtmb.sample:::draws_model_weights(
    list(d1, d2), list(method = m), NULL, NULL, "x")))
  cat(m, "old", format(old, digits = 17), "new", format(new, digits = 17),
      "identical", identical(old, new), "\n")
}
