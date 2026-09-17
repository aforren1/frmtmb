.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.learn))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: igt_orl
### Title: Outcome-representation learning for the Iowa gambling task
### Aliases: igt_orl

### ** Examples

d <- frm_task_design("igt", n_subject = 6, n_trial = 60, seed = 8)
d$choice <- frm_task_simulate(
  igt_orl(subject = id, trial = trial), d,
  pars = list(Arew = 0.3, Apun = 0.1, k = 0.5, betaF = 1,
              betaP = 1), seed = 8)[[1]]$choice
table(d$choice)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
