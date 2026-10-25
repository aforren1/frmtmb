# Smoke test of fcor() in core: fit, logLik against a dense MVN at the
# estimates, and the post-fit methods. Seed 1.
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
set.seed(1)
n <- 40
V <- 0.6^abs(outer(1:n, 1:n, "-"))
d <- data.frame(x = rnorm(n), g = factor(rep(1:8, each = 5)))
L <- t(chol(V))
d$y <- 1 + 0.5 * d$x + as.vector(L %*% rnorm(n)) * 0.8
dense <- function(y, mu, S) {
  R <- chol(S)
  z <- backsolve(R, y - mu, transpose = TRUE)
  -0.5 * n * log(2 * pi) - sum(log(diag(R))) - 0.5 * sum(z^2)
}
f1 <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V))
print(f1)
mu <- fitted(f1)[, "Estimate"]
s <- sigma(f1)
cat("logLik", format(as.numeric(logLik(f1)), digits = 12),
    "dense", format(dense(d$y, mu, s^2 * V), digits = 12), "\n")
# a covariance M: sigma scales it
f2 <- frm(y ~ x + fcor(V2), data = d, data2 = list(V2 = 4 * V))
cat("cov M: sigma", sigma(f2), "vs", sigma(f1) / 2, " logLik",
    as.numeric(logLik(f2)), as.numeric(logLik(f1)), "\n")
# predicted sigma
f3 <- frm(bf(y ~ x + fcor(V), sigma ~ x), data = d, data2 = list(V = V))
dp <- frm_linpred(f3, dpar = "sigma", type = "response")
mu3 <- fitted(f3)[, "Estimate"]
cat("het logLik", as.numeric(logLik(f3)), "dense",
    dense(d$y, mu3, diag(dp) %*% V %*% diag(dp)), "\n")
# student
f4 <- frm(y ~ x + fcor(V), family = student(), data = d,
          data2 = list(V = V))
print(summary(f4))
# random effect
f5 <- frm(y ~ x + (1 | g) + fcor(V), data = d, data2 = list(V = V))
print(VarCorr(f5))
print(head(residuals(f1, type = "pearson")))
print(dim(autocor_matrix(f1)))
print(head(simulate(f1, nsim = 2, seed = 3)))
print(head(predict(f1)))
r <- try(fitted(f1, newdata = d[1:5, ]))
r <- try(simulate(f1, newdata = d[1:5, ]))
print(head(frm_linpred(f1, newdata = d[1:5, ])))
print(confint_varcorr(f1))
print(variables(f1))
print(frm_compat("fcor()", "gaussian"))
