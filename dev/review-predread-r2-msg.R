# Short check of fix round 2: the whole text of the warning on the build with the fixes.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7411)
d <- data.frame(x = rnorm(120), xo = rnorm(120), sx = 0.4, g = factor(rep(1:12, 10)))
d$y <- 1 + d$x + rep(rnorm(12), 10) + rnorm(120)
fit <- quiet(frm(bf(y ~ me(xo, sx) + x + (1 | g)), data = d))
m <- tryCatch(marginaleffects::get_vcov(fit), warning = function(w) conditionMessage(w))
cat(strwrap(m, 80), sep = "\n")
