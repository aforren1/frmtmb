.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.eam))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: lba_simulate
### Title: Simulate from a linear ballistic accumulator
### Aliases: lba_simulate

### ** Examples

set.seed(1)
dat <- lba_simulate(500, v = c(2.5, 1.5, 1.0), A = 0.5, k = 0.4,
                    ndt = 0.2)
table(dat$choice)
tapply(dat$rt, dat$choice, mean)




'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
