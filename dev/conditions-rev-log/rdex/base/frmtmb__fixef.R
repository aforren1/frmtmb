.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: fixef.frmtmb_fit
### Title: Extract fixed effects
### Aliases: fixef.frmtmb_fit fixef

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)

# one entry per distributional parameter, each on its link scale
fit <- frm(bf(y ~ x + (1 | g), sigma ~ x) + gaussian(), data = dd)
fixef(fit)
exp(fixef(fit)$sigma[["(Intercept)"]])   # sigma is modeled on the log

# flatten to the vector vcov() and confint() name their rows by
fixef(fit, flatten = TRUE)
all(names(fixef(fit, flatten = TRUE)) %in% rownames(confint(fit)))

# so a standard error goes with its coefficient by name
cf <- fixef(fit, flatten = TRUE)
cf / sqrt(diag(vcov(fit))[names(cf)])



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
