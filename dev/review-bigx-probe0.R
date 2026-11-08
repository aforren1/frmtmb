source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-common.R")
set.seed(1); d <- data.frame(z = rnorm(100), g = factor(rep(1:10, 10)))
d$y <- rpois(100, exp(0.5 + 0.4 * d$z)); d$w <- d$z + 1e6
f <- frm(bf(y ~ w + (1 | g)), family = poisson(), data = d)
str(f$estimates[1:3]); print(names(f$opt$par)); print(head(fixef(f)))
print(class(vcov(f))); print(dimnames(vcov(f)))
print(exists("fixef")); print(ls("package:frmtmb")[1:400])
