.libPaths(c("C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb); library(frmtmb.eam)
set.seed(5)
d <- do.call(rbind, lapply(1:6, function(s) {
  x <- ddm_simulate(100, mu = 0.8 + rnorm(1, 0, 0.3), bs = 1.4, ndt = 0.25)
  x$s <- factor(s); x$cond <- rep(0:1, 50); x
}))
d$cens <- ifelse(d$rt > 1, "right", "none"); d$rt <- pmin(d$rt, 1)
t0 <- proc.time()
fit <- frm(bf(rt | dec(upper) + cens(cens) ~ cond + (1 | s), bs ~ 1,
              ndt ~ 1, bias = 0.5), family = wiener(), data = d)
print(proc.time() - t0)
str(fixef_by_dpar(fit))
print(names(fit$sdr$par.fixed))
print(fit$sdr$par.fixed)
print(confint(fit))
print(VarCorr(fit))
