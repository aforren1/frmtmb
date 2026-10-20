# Review of lane illse2, finding 7: which SEs the poisson GLMM with a
# covariate on scale 500 loses, and what glmer reports.
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
set.seed(13)
n <- 300
d <- data.frame(x = rnorm(n, 500, 100), g = factor(rep(1:20, length.out = n)))
d$y <- rpois(n, exp(0.3 + 0.004 * (d$x - 500) + rnorm(20, 0, 0.3)[d$g]))
r <- cap(frm(bf(y ~ x + (1 | g)) + poisson(), data = d))
f <- r$v
cat("lost:", lost_of(f), "\n")
for (w in r$w) cat("W:", substr(w, 1, 250), "\n")
print(suppressWarnings(fixef(f)))
l <- lme4::glmer(y ~ x + (1 | g), data = d, family = poisson)
print(summary(l)$coefficients)
