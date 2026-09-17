.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: diagnose
### Title: Convergence diagnostics for a frmtmb fit
### Aliases: diagnose

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
diagnose(fit)

# a random effect the data cannot support collapses to the boundary,
# which is a valid fit but a warning about the model
dd$h <- factor(rep(1:10, each = 10))
fit_s <- frm(bf(y ~ x + (1 | g) + (1 | h)) + gaussian(), data = dd)
d <- diagnose(fit_s, quiet = TRUE)
d$singular

# a predictor scaled far from one slows the optimizer down; the
# remedy is frmtmb_control(autoscale = TRUE)
dd$xbig <- dd$x * 1e5
diagnose(frm(bf(xbig ~ 1) + gaussian(), data = dd), quiet = TRUE)$scale



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
