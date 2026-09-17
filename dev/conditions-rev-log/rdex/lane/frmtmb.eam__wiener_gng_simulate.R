.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.eam))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: wiener_gng_simulate
### Title: Simulate from a go/no-go diffusion model
### Aliases: wiener_gng_simulate

### ** Examples

set.seed(1)
dat <- wiener_gng_simulate(1000, mu = 1.0, bs = 1.4, ndt = 0.25,
                           deadline = 1.5)
mean(dat$responded)
summary(dat$rt[dat$responded == 1])




'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
