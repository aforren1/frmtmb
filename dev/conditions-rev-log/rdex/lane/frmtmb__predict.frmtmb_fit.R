.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: predict.frmtmb_fit
### Title: Predictions from a frmtmb fit
### Aliases: predict.frmtmb_fit

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)

# the link scale by default; "response" is what fitted() returns
head(predict(fit))
max(abs(predict(fit, type = "response") - fitted(fit)))

# re_formula = NA drops the random effects: the population prediction
nd <- data.frame(x = c(-1, 0, 1), g = factor(1, levels = levels(dd$g)))
predict(fit, newdata = nd, re_formula = NA, type = "response")

# delta-method standard errors, on whichever scale was asked for
p <- predict(fit, newdata = nd, se.fit = TRUE)
cbind(fit = p$fit, se = p$se.fit)

# a level the fit never saw errors unless it is allowed explicitly,
# in which case it is predicted at the population level
nd_new <- data.frame(x = 0, g = factor("new"))
try(predict(fit, newdata = nd_new))
predict(fit, newdata = nd_new, allow_new_levels = TRUE)

# a distributional parameter instead of the mean
fit2 <- frm(bf(y ~ x, sigma ~ x) + gaussian(), data = dd)
head(predict(fit2, dpar = "sigma", type = "response"))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
