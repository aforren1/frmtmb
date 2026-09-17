.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: par_template
### Title: Parameter names and starting values
### Aliases: par_template par_template.frmtmb_fit par_template.default

### ** Examples

dd <- data.frame(x = rnorm(30), g = factor(rep(1:5, 6)))
dd$y <- 1 + 2 * dd$x + rnorm(30)

# before fitting: the names and the cold starting values
par_template(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# after fitting: the same layout, holding the estimates
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), dd)
par_template(fit)

# discover, edit, fit
st <- par_template(bf(y ~ x) + gaussian(), data = dd)
st$beta["x"] <- 2
frm(bf(y ~ x) + gaussian(), dd, start = st)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
