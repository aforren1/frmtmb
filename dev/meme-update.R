# Lane meme: update() carries a me() prior and drops it on a model
# without me(), as brms's update() does. Output: dev/meme-out/update.txt
.libPaths(c("C:/Users/adf44/source/r/wt-meme-lib", "C:/Users/adf44/source/r/rellib-r15", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(4)
d <- data.frame(x = rnorm(50), z = rnorm(50), sx = 0.3, sz = 0.3)
d$y <- 1 + d$x - d$z + rnorm(50)
fit <- frm(y ~ me(x, sx) + me(z, sz), data = d,
           prior = set_prior("normal(0, 1)", class = "meanme") +
             set_prior("lkj(2)", class = "corme"))
print(prior_summary(fit))
f2 <- update(fit, . ~ . - me(z, sz))
print(prior_summary(f2))
f3 <- tryCatch(update(fit, . ~ x), error = function(e) conditionMessage(e))
print(if (is.character(f3)) f3 else prior_summary(f3))
print(anova(fit, update(fit, prior = set_prior("normal(0, 1)", class = "meanme") + set_prior("lkj(2)", class = "corme"))))
