.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1)
dd <- data.frame(x = rnorm(100))
G <- matrix(rgamma(400, 2), 100); dd$Y <- G / rowSums(G)
colnames(dd$Y) <- letters[1:4]
f <- dirichlet(refcat = "c")
fr <- f[["defer"]](bf(Y ~ x)$formula, dd)
print(fr$dpars)
r <- tryCatch(frm(bf(Y ~ x, mua ~ 1), family = dirichlet(refcat = "c"), data = dd), error = function(e) {print(sys.calls()); conditionMessage(e)})
print(r)
r2 <- tryCatch(frm(bf(Y ~ x, mub ~ 1), family = dirichlet(refcat = "c"), data = dd), error = function(e) conditionMessage(e))
print(r2)
