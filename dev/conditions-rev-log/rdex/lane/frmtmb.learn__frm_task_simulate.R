.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.learn))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_task_simulate
### Title: Draw datasets from a learning family's generative process
### Aliases: frm_task_simulate

### ** Examples

d <- frm_task_design("bandit2arm", n_subject = 5, n_trial = 30,
                     seed = 1)
sims <- frm_task_simulate(bandit2arm_delta(subject = id, trial = trial),
                          d, pars = list(alpha = 0.4, tau = 3),
                          nsim = 2, seed = 1)
table(sims[[1]]$choice)

# a learning rate that changes within a subject, one value per row
r <- frm_task_design("reversal", n_subject = 5, n_trial = 30, seed = 2)
a <- ifelse(r$after_reversal == "after", 0.6, 0.2)
rs <- frm_task_simulate(bandit2arm_delta(subject = id, trial = trial),
                        r, pars = list(alpha = a, tau = 3), seed = 2)
table(rs[[1]]$choice)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
