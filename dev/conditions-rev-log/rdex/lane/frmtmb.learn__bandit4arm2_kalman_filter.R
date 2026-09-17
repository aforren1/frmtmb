.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.learn))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: bandit4arm2_kalman_filter
### Title: A Kalman filter over the restless four-armed bandit
### Aliases: bandit4arm2_kalman_filter

### ** Examples

d <- frm_task_design("bandit4arm_restless", n_subject = 5, n_trial = 40,
                     seed = 6)
d$choice <- frm_task_simulate(
  bandit4arm2_kalman_filter(subject = id, trial = trial), d,
  pars = list(lambda = 0.98, center = 50, tau = 0.15, mu0 = 50,
              sigma0 = 10, sigmaD = 3), seed = 6)[[1]]$choice
table(d$choice)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
