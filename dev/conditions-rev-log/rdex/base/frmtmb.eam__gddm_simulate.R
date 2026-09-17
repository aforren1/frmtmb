.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.eam))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: gddm_simulate
### Title: Simulate from a generalized drift-diffusion model
### Aliases: gddm_simulate

### ** Examples

set.seed(2)
head(gddm_simulate(20, mu = 1.5, bs = 2, ndt = 0.2,
                   control = gddm_control(t_max = 2, dt = 0.02,
                                          ny = 101)))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
