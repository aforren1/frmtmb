# frm_sample() on gp(x, by = f): names and scales of the draws.
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-gpby/dev/stan-cache")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
set.seed(1)
n <- 60
dd <- data.frame(x = runif(n, 0, 5),
                 f = factor(sample(c("a", "b"), n, TRUE)))
dd$y <- ifelse(dd$f == "a", sin(dd$x), cos(dd$x)) + rnorm(n, 0, 0.3)
fit <- frm(bf(y ~ gp(x, by = f)), data = dd)
print(variables(fit))
ds <- frm_sample(fit, chains = 1, iter = 400, warmup = 200, seed = 1,
                 refresh = 0)
v <- variables(ds)
print(v[!grepl("^b\\[", v)])
print(summary(ds))
nd <- data.frame(x = c(1, 6), f = factor(c("a", "b")))
ep <- posterior_epred(ds, newdata = nd)
print(apply(ep, 2, quantile, c(0.05, 0.5, 0.95)))
