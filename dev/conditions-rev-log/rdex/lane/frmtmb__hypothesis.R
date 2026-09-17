.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: hypothesis
### Title: Hypothesis tests on parameter expressions
### Aliases: hypothesis hypothesis.frmtmb_fit hypothesis.frmtmb_multiple

### ** Examples

set.seed(4)
dd <- data.frame(x1 = rnorm(120), x2 = rnorm(120),
                 g = factor(rep(1:10, 12)))
dd$y <- rnorm(120, 1 + 0.6 * dd$x1 + 0.4 * dd$x2 +
                rnorm(10, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x1 + x2 + (1 | g)) + gaussian(), data = dd)
h <- hypothesis(fit, c("x1 - x2 = 0", "exp(Intercept) = 1"))
h
h$hypothesis$Est.Error
attr(h, "test")$p
# brms's directional form: one-sided p, and a 90% interval
hypothesis(fit, "x1 > x2")
# class/group name the natural-scale random-effect summaries
hypothesis(fit, "Intercept > 0", class = "sd", group = "g")
# variance-component expressions need class = NULL, as in brms: an
# ICC with bootstrap intervals
hypothesis(fit, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2) = 0",
           class = NULL, method = "boot", nsim = 20, seed = 1)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
