.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: prior
### Title: Set up priors with brms's quoting spelling
### Aliases: prior prior_ prior_string

### ** Examples

# the brms nonlinear vignette's spelling
prior(normal(5000, 1000), nlpar = "ult")

# combine with c() or `+`, as with set_prior()
c(prior(normal(1, 2), nlpar = "omega"),
  prior(normal(45, 10), nlpar = "theta"))

# the programmatic spellings
prior_(~normal(0, 10), class = ~b)
prior_string(paste0("normal(0, ", 2 * 5, ")"), class = "b")



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
