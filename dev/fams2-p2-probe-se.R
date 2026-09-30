# Punch 2, n3: how a fit-end check can read the standard errors of
# kappa's coefficients (names, and the component and index route).
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(31)
n <- 400
x <- rnorm(n)
d <- data.frame(y = rbeta(n, 2, 5), x = x)
f <- suppressWarnings(frm(bf(y ~ 1, kappa ~ x), family = xbeta(), data = d))
lp <- f$frame$linpreds[[frmtmb:::linpred_key("y", "kappa")]]
str(lp[c("par", "idx")])
sdr <- frmtmb:::sdr_of(f)
print(rownames(sdr$cov.fixed))
print(frmtmb:::outer_par_names(f))
print(names(f$opt$par))
print(sqrt(diag(sdr$cov.fixed)))
print(range(frmtmb:::eval_dpars(f)[["y"]][["kappa"]]))
