# drop_unused_levels = FALSE on a grouping factor with an unused level.
.libPaths(c("C:/Users/adf44/source/r/wt-formrobust-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(7)
d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 10),
                           levels = c("a", "b", "c", "d")), x = rnorm(30))
d$y <- rnorm(30, as.numeric(d$g))
f1 <- frm(bf(y ~ x + (1 | g)), data = d)
f0 <- frm(bf(y ~ x + (1 | g)), data = d, drop_unused_levels = FALSE)
print(rownames(ranef(f1)$g)); print(rownames(ranef(f0)$g))
cat("logLik rel diff:", format(abs(as.numeric(logLik(f0) - logLik(f1))) /
                               abs(as.numeric(logLik(f1))), digits = 3), "\n")
print(ranef(f0)$g)
