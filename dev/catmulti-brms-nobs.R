# brms 2.23.0: nobs() and the summary's observation count of a subset()
# model, for the stretch item of lane catmulti; and drop1() /
# marginaleffects on a multinomial brms fit.
# Output: dev/catmulti-out/brms-nobs.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
set.seed(3)
n <- 300
dd <- data.frame(x = rnorm(n))
dd$y <- 1 + dd$x + rnorm(n)
dd$keep <- c(FALSE, rep(TRUE, n - 1))
fit <- brm(bf(y | subset(keep) ~ x), data = dd, chains = 1, iter = 400,
           refresh = 0, backend = "rstan", seed = 1)
cat("nobs:", nobs(fit), "\n")
s <- summary(fit)
cat("summary nobs:", s$nobs, "\n")
print(s)
cat("log_lik columns:", ncol(log_lik(fit)), "\n")
