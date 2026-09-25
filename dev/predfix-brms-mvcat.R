# brms 2.23.0: fitted() on a multivariate model with a categorical
# response (punch round 1, minor 5). Same data as test-predfix.R.
#   R_MAKEVARS_USER=... Rscript dev/predfix-brms-mvcat.R > dev/predfix-log/brms-mvcat.txt
.libPaths(c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
set.seed(4)
n <- 120
d <- data.frame(x = rnorm(n))
d$y1 <- 1 + 0.5 * d$x + rnorm(n)
d$o <- cut(d$x + rnorm(n), c(-Inf, -0.5, 0.5, Inf), labels = c("a", "b", "c"))
# mvbf() with a family per response: `bf() + gaussian() + bf() +
# categorical()` gives the LAST family to both responses
b <- brm(mvbf(bf(y1 ~ x, family = gaussian()),
              bf(o ~ x, family = categorical())), data = d,
         chains = 1, iter = 600, refresh = 0, seed = 1,
         backend = "cmdstanr")
# rstan's stanc.js fails on this model without V8 ("parser failed
# badly"), so it is compiled by CmdStan 2.39.0 instead
f <- fitted(b)
cat("dim", dim(f), "\n")
print(dimnames(f)[2:3])
f2 <- fitted(b, resp = c("o", "y1"))
cat("resp = c('o', 'y1'): dim", dim(f2), "\n")
print(dimnames(f2)[[3]])
f3 <- fitted(b, resp = "o")
cat("resp = 'o': dim", dim(f3), "\n")
print(dimnames(f3)[[3]])
print(round(f[1:2, , ], 4))
