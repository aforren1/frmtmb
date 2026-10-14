# The GEV mean and variance against the same formulas in 256-bit
# arithmetic, which separates the formula's rounding from the
# quadrature tolerance in dev/famuni-out/moments.txt.
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib", "C:/Users/adf44/source/r/rellib-r8", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(Rmpfr)})
P <- unclass(gen_extreme_value())$post
xi <- c(-0.6, -0.2, -1e-5, -1e-6, -1e-8, 1e-8, 1e-6, 1e-5, 0.1, 0.4, 0.8)
m <- P$mean_fn(list(mu = rep(0, length(xi)), sigma = rep(1, length(xi)), xi = xi), list())
x <- mpfr(xi, 256)
rm_ <- as.numeric((gamma(1 - x) - 1) / x)
cat("mean: max rel err vs 256-bit", format(max(abs(m - rm_) / abs(rm_)), digits = 3), "over", length(xi), "xi\n")
ok <- xi < 0.5
v <- P$var_fn(list(mu = rep(0, sum(ok)), sigma = rep(1, sum(ok)), xi = xi[ok]), list())
rv <- as.numeric((gamma(1 - 2 * x[ok]) - gamma(1 - x[ok])^2) / x[ok]^2)
cat("var: max rel err vs 256-bit", format(max(abs(v - rv) / abs(rv)), digits = 3), "over", sum(ok), "xi\n")
