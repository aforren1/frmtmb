.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.eam))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: rdm_simulate
### Title: Simulate from a racing diffusion model
### Aliases: rdm_simulate

### ** Examples

set.seed(1)
dat <- rdm_simulate(500, v = c(3.0, 2.0, 1.2), A = 0.5, k = 0.5,
                    ndt = 0.2)
table(dat$choice)
tapply(dat$rt, dat$choice, mean)




'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
