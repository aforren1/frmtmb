# Smoke test of gp(x, by = ) in the lane build.
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(1)
n <- 90
dd <- data.frame(x = runif(n, 0, 5), z = runif(n, 0, 3),
                 f = factor(sample(c("a", "b", "c"), n, TRUE)),
                 w = runif(n, 0.5, 2))
dd$y <- ifelse(dd$f == "a", sin(dd$x), ifelse(dd$f == "b", cos(dd$x),
                                               0.3 * dd$x)) +
  rnorm(n, 0, 0.3)
f1 <- frm(bf(y ~ gp(x, by = f)), data = dd)
print(logLik(f1))
print(variables(f1))
print(summary(f1)$gp)
nd <- data.frame(x = c(0.5, 2.5, 6), f = factor(c("a", "b", "c")))
print(frm_linpred(f1, newdata = nd, se.fit = TRUE))
lb <- frm_lp_basis(f1, newdata = nd, extra_cov = TRUE)
print(lb$extra_var); print(lb$extra_cov)
r <- try(frm_linpred(f1, newdata = data.frame(x = 1, f = "d")))
f2 <- frm(bf(y ~ gp(x, by = f, k = 10)), data = dd)
print(logLik(f2)); print(summary(f2)$gp)
f3 <- frm(bf(y ~ gp(x, by = w)), data = dd)
print(logLik(f3)); print(variables(f3))
f4 <- frm(bf(y ~ gp(x, z, by = f, iso = FALSE)), data = dd)
print(summary(f4)$gp)
f5 <- frm(bf(y ~ gp(x, by = f, cmc = FALSE)), data = dd)
print(variables(f5))
f6 <- frm(bf(y ~ gp(x, by = f, gr = FALSE)), data = dd)
print(c(logLik(f1), logLik(f6)))
