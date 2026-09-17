.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.learn))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: ts_par7
### Title: The two-step task: a model-based and model-free hybrid
### Aliases: ts_par7

### ** Examples

d <- frm_task_design("twostep", n_subject = 5, n_trial = 40, seed = 8)
d <- frm_task_simulate(
  ts_par7(subject = id, trial = trial), d,
  pars = list(alpha1 = 0.4, tau1 = 3, alpha2 = 0.4, tau2 = 3,
              lambda = 0.6, w = 0.5, pers = 0.2), seed = 8)[[1]]
head(d[, c("id", "trial", "choice", "state2", "choice2")])



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
