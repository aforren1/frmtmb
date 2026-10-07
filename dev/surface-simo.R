# Lane surface: the draws names of a mo() model (is there a simo_ column
# for plot(variable = "simo", regex = TRUE) to find?).
source("dev/surface-env.R")
surface_env("lane")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(1)
income <- factor(sample(c("below_20", "20_to_40", "40_to_100", "greater_100"), 100, TRUE),
                 levels = c("below_20", "20_to_40", "40_to_100", "greater_100"), ordered = TRUE)
ls <- c(30, 60, 70, 75)[income] + rnorm(100, 0, 10)
dat <- data.frame(income, ls)
fit1 <- frm(ls ~ mo(income), data = dat)
ds <- suppressWarnings(frm_sample(fit1, chains = 1, iter = 300, refresh = 0, seed = 1))
print(variables(ds))
print(grep("simo", variables(ds), value = TRUE))
