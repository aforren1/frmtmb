.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.learn))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: igt_pvl_delta
### Title: Prospect-valence learning with a delta rule, for the Iowa
###   gambling task
### Aliases: igt_pvl_delta

### ** Examples

d <- frm_task_design("igt", n_subject = 6, n_trial = 50, seed = 7)
d$choice <- frm_task_simulate(
  igt_pvl_delta(subject = id, trial = trial), d,
  pars = list(alpha = 0.3, shape = 0.4, lambda = 1.5, tau = 1),
  seed = 7)[[1]]$choice
# the good decks are 3 and 4
table(d$choice)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
