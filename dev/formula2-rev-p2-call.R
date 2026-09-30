# Reviewer, punch round 2: eval(fit$call) and update() refit the same model.
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(6)
d <- data.frame(y = rpois(120, 3), x1 = rnorm(120), x2 = rnorm(120),
                g = factor(rep(1:6, each = 20)))
fit <- frm(y ~ ., data = d, family = poisson())
print(fit$call)
d$late <- rnorm(120)
cat("eval(call) logLik identical:", identical(logLik(eval(fit$call)),
                                               logLik(fit)), "\n")
cat("update() logLik identical:", identical(logLik(update(fit)),
                                             logLik(fit)), "\n")
fb <- frm(bf(y ~ ., sigma ~ x1), data = transform(d, y = y + rnorm(120)))
cat("bf() call formula class:", class(fb$call$formula), "\n")
