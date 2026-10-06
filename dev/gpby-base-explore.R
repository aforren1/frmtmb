# How 0.67.0 names and reports a gp() fit, so the by = names can match.
lib <- Sys.getenv("GPBY_LIB", "")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
set.seed(1)
n <- 40
dd <- data.frame(x = runif(n, 0, 5), z = runif(n, 0, 3),
                 f = factor(sample(c("a", "b", "c"), n, TRUE)),
                 w = runif(n, 0.5, 2))
dd$y <- sin(dd$x) + rnorm(n, 0, 0.3)
fit <- frm(bf(y ~ gp(x)), data = dd)
print(variables(fit))
s <- summary(fit)
print(names(s))
print(s$gp)
print(VarCorr(fit))
print(fixef(fit))
print(default_prior(bf(y ~ gp(x)), data = dd))
r <- try(frm(bf(y ~ gp(x, by = f)), data = dd))
r <- try(frm(bf(y ~ gp(x, z)), data = dd)); print(VarCorr(r))
ds <- frm_sample(bf(y ~ gp(x)), data = dd, chains = 1, iter = 200,
                 warmup = 100, seed = 1, refresh = 0)
print(head(variables(ds), 20))
print(summary(ds))
