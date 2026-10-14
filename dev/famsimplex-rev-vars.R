# Reviewer: variables(), hypothesis(), confint() on ML fits; CE dpar se.
# Output: dev/famsimplex-rev-vars.log
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-sim.R")
set.seed(1)
dg <- data.frame(x = rnorm(200), z = rnorm(200))
dg$y <- 1 + dg$x + rnorm(200, 0, exp(0.2 * dg$z))
fg <- frm(bf(y ~ x), family = gaussian(), data = dg)
cat("gaussian variables:", variables(fg), "\n")
fg2 <- frm(bf(y ~ x, sigma ~ z), family = gaussian(), data = dg)
ce <- conditional_effects(fg2, "z", dpar = "sigma")[[1]]
cat("gaussian CE dpar sigma: est", ce$estimate__[1], "se", ce$se__[1],
    "lo", ce$lower__[1], "hi", ce$upper__[1], "\n")
dl <- fs_data_logistic_normal(seed = 13)
fl <- frm(bf(Y ~ x), family = logistic_normal(), data = dl)
cat("ln variables:", variables(fl), "\n")
print(hypothesis(fl, "lncor__b__c > 0", class = NULL))
print(hypothesis(fl, "sigmac > 0.5", class = NULL))
print(confint(fl))
print(summary(fl)$spec_pars)
# K = 4: three correlations, names and order against brms's loop
set.seed(6)
d6 <- data.frame(x = rnorm(300))
M <- cbind(0.2 + 0.4 * d6$x, -0.1 * d6$x, 0.5 + 0 * d6$x)
C <- matrix(c(1, 0.3, -0.2, 0.3, 1, 0.4, -0.2, 0.4, 1), 3)
d6$Y <- fs_rlogistic_normal(M, diag(c(0.5, 0.7, 0.9)) %*% C %*%
                              diag(c(0.5, 0.7, 0.9)))
colnames(d6$Y) <- c("a", "b", "c", "d")
f6 <- frm(bf(Y ~ x), family = logistic_normal(), data = d6)
print(summary(f6)$spec_pars)
cat("empirical alr cor (truth b,c 0.3; b,d -0.2; c,d 0.4):\n")
X <- log(d6$Y[, -1]) - log(d6$Y[, 1])
R <- resid(lm(X ~ d6$x))
print(round(cor(R), 4))
f6c <- frm(bf(Y ~ x), family = logistic_normal(refcat = "c"), data = d6)
print(summary(f6c)$spec_pars)
cat(sprintf("logLik ref a %.10f ref c %.10f\n", logLik(f6), logLik(f6c)))
