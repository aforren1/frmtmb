# Lane adhess: how far the covariance of test-se-check.R's "healthy fit"
# (ls ~ mo(income), se_mo_data(7), ls / 16) moves from 0.71.0's
# (optimHess() differences) to the inverse of the AD Hessian.
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-common.R")
set.seed(7)
lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
d <- data.frame(income, ls)
d$age <- rnorm(100, mean = 40, sd = 10)
d$ls <- d$ls / 16
fit <- frm(ls ~ mo(income), data = d)
old <- RTMB::sdreport(fit$obj)$cov.fixed
H <- fit$obj$he(ns$sdr_outer_point(fit$obj))
new <- solve((H + t(H)) / 2)
print(cov_err(old, new))
