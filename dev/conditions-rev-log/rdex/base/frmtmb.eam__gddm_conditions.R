.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.eam))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: gddm_conditions
### Title: Build a condition index for 'gddm()'
### Aliases: gddm_conditions

### ** Examples

d <- data.frame(coh = c(0, 0, 0.5, 0.5), block = c(1, 2, 1, 2))
gddm_conditions(d, coh, block)
gddm_conditions(d, ~ coh)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
