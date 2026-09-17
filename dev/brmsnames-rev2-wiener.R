## Reviewer recheck, MAJOR 3 on an extension family: wiener()'s natural
## columns (bs, ndt through a closure-bound scaled logit, bias). The
## readers that hand a draw to the model must give base's numbers.
##   Rscript dev/brmsnames-rev2-wiener.R base|lane
## Data seed 77 (the eam sampling test's), sampler seed 3.
arm <- commandArgs(trailingOnly = TRUE)[1L]
libs <- list(
  base = c("C:/Users/adf44/source/r/rellib-r3"),
  lane = c("C:/Users/adf44/source/r/brmsnames-lib",
           "C:/Users/adf44/source/r/rellib-r3"))
.libPaths(c(libs[[arm]], "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample)); q(library(frmtmb.eam))
set.seed(77)
dat <- ddm_simulate(250, mu = 0.9, bs = 1.4, ndt = 0.25)
fit <- q(frm(bf(rt | vint(upper) ~ 1), family = wiener(), data = dat))
set.seed(3)
ds <- q(frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 3))
cat("column means:\n"); print(colMeans(ds$draws), digits = 10)
ll <- q(log_lik(ds))
cat(sprintf("sum(log_lik(ds)) %.12f\n", sum(ll)))
ep <- q(posterior_epred(ds))
cat(sprintf("sum(posterior_epred(ds)) %.12f\n", sum(ep)))
set.seed(4)
pp <- q(posterior_predict(ds, ndraws = 20))
cat(sprintf("sum(posterior_predict(ds, ndraws = 20)) %.12f\n", sum(pp)))
lo <- q(loo(ds))$estimates
print(lo, digits = 12)
