.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.eam))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: lba
### Title: The linear ballistic accumulator
### Aliases: lba

### ** Examples

set.seed(1)
dat <- lba_simulate(400, v = c(2.4, 1.6, 1.0), A = 0.5, k = 0.4,
                    ndt = 0.2)
fit <- frm(bf(rt | vint(choice) ~ 1), family = lba(3), data = dat)
fixef(fit)




'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
