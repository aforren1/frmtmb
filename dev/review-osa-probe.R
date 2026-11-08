# Review lane osa: probe how a fit exposes its estimates (scratch).
.libPaths(c("C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
set.seed(5)
d <- expand.grid(t = 1:10, g = factor(1:20))
d$x <- rnorm(nrow(d))
d$y <- rpois(nrow(d), exp(0.3 + 0.4 * d$x + rnorm(20, 0, .8)[d$g]))
f1 <- frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d)
str(f1$estimates)
print(summary(f1))
f2 <- frm(y ~ x + (1 | g), family = zero_inflated_poisson(), data = d)
str(f2$estimates)
str(lapply(frmtmb:::eval_dpars(f2)[[1]], function(v) head(as.numeric(v))))
print(fixef(f2))
