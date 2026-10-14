# Does brms 2.23.0 let the FIRST predicted category take its own formula?
# Output: dev/famsimplex-brms-first.log
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(1)
dd <- data.frame(x = rnorm(60), z = rnorm(60))
G <- matrix(rgamma(180, 2), 60); dd$Y <- G / rowSums(G)
colnames(dd$Y) <- c("a", "b", "c")
dd$f <- factor(sample(c("a", "b", "c"), 60, TRUE))
grab <- function(code, pat) grep(pat, strsplit(code, "\n")[[1]], value = TRUE)
c1 <- stancode(bf(Y ~ x, mub ~ z), data = dd, family = dirichlet())
print(grab(c1, "X_mub|X_muc"))
c2 <- tryCatch(stancode(bf(f ~ x, mub ~ z), data = dd, family = categorical()),
               error = function(e) paste("ERROR:", conditionMessage(e)))
print(grab(c2, "X_mub|X_muc"))
c3 <- stancode(bf(Y ~ x, mua ~ z), data = dd,
               family = dirichlet(refcat = "c"))
print(grab(c3, "X_mua|X_mub"))
print(brms::default_prior(bf(Y ~ x, mub ~ z), data = dd,
                          family = dirichlet())[, c("class", "coef", "dpar")])
# the mixture spelling the vigport ledger records
set.seed(3)
dm <- data.frame(x1 = rnorm(60), y = rnorm(60))
cm <- stancode(bf(y ~ 1, mu1 ~ x1, mu2 ~ x1), data = dm,
               family = mixture(gaussian, gaussian))
print(grab(cm, "X_mu1|X_mu2"))
