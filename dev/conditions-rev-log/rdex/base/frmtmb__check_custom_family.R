.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: check_custom_family
### Title: Check a custom family's log-density for AD safety
### Aliases: check_custom_family

### ** Examples

set.seed(1)
y <- rpois(50, 3)

# a hand-written poisson: check it before fitting anything with it
ok <- custom_family(
  "my_poisson", dpars = "mu", links = list(mu = "log"),
  lpdf = function(y, dpars, aterms) {
    y * log(dpars$mu) - dpars$mu - lgamma(y + 1)
  },
  type = "discrete"
)
check_custom_family(ok, y = y, dpars = list(mu = rep(2.5, 50)))

# base matrix() strips the advector class, so the tape sees constants
# and the gradient is silently wrong. The check catches it.
bad <- custom_family(
  "bad", dpars = "mu", links = list(mu = "log"),
  lpdf = function(y, dpars, aterms) {
    m <- matrix(dpars$mu, ncol = 1)
    y * log(m[, 1]) - m[, 1] - lgamma(y + 1)
  },
  type = "discrete"
)
try(check_custom_family(bad, y = y, dpars = list(mu = rep(2.5, 50))))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
