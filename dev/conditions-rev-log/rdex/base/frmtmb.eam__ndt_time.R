.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.eam))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: ndt_time
### Title: The non-decision time, in the response's own units
### Aliases: ndt_time

### ** Examples

set.seed(1)
d <- ddm_simulate(300, mu = 1.2, bs = 1.5, ndt = 0.25)
fit <- frm(bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5),
           family = wiener(), data = d)
# with no ndt_group() the two agree: `ndt` is already a time
head(predict(fit, dpar = "ndt", type = "response"), 3)
head(ndt_time(fit), 3)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
