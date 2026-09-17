.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.eam))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: ddm_simulate
### Title: Simulate a drift-diffusion data set
### Aliases: ddm_simulate

### ** Examples

set.seed(1)
cond <- rep(c(0, 1), each = 100)
dat <- ddm_simulate(200, mu = 0.2 + 1.1 * cond, bs = 1.4,
                    ndt = 0.3, bias = 0.5)
dat$cond <- factor(cond)
str(dat)

# Ratcliff's full model, at values the literature uses
full <- ddm_simulate(200, mu = 1.2, bs = 1.5, ndt = 0.3,
                     sv = 1.0, sz = 0.2, st = 0.1)




'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
