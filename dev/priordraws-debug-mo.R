# Why the mo() simplex has no prior draw statement: the names the plan
# reads.
.libPaths(c("C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(1)
d <- data.frame(x = rnorm(60), m = sample(1:4, 60, TRUE))
d$y <- 1 + 0.5 * d$x + rnorm(60)
fit <- frm(bf(y ~ mo(m)), data = d)
nc <- frmtmb.sample:::draws_natural_cols(fit)
str(nc$ordinal, max.level = 2)
print(nc$extra)
print(nc$derived$names)
tab <- validate_prior(NULL, fit, route = "fit")
print(as.data.frame(tab))
