# Reviewer: where a cumulative fit keeps its thresholds, and their sign.
.libPaths(c("C:/Users/adf44/source/r/wt-fuzzfix-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
n <- 300; x <- rnorm(n)
y <- factor(cut(0.5 * x + rlogis(n), c(-Inf, -0.5, 0.7, Inf), labels = FALSE),
            ordered = TRUE)
d <- data.frame(x = x, y = y)
f <- frm(bf(y ~ x), family = cumulative(), data = d)
str(f$estimates)
print(rownames(vcov(f)))
o <- MASS::polr(y ~ x, data = d)
print(o$zeta); print(coef(o)); print(c(logLik(o), logLik(f)))
