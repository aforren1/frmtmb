.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb_family
### Title: Define a model family
### Aliases: frmtmb_family custom_family

### ** Examples

# a custom family is a plain R log-density over taped parameters
dd <- data.frame(y = rbinom(100, 5, 0.4),
                 size = 5, x = rnorm(100))
fam <- custom_family(
  "vbinom", dpars = "mu", links = list(mu = "logit"),
  lpdf = function(y, dpars, aterms) {
    RTMB::dbinom(y, aterms$vint1, dpars$mu, log = TRUE)
  },
  # the density indexes vint1, so a model without it is refused
  # rather than fitted against a zero-length log-likelihood
  required_aterms = "vint1",
  type = "discrete"
)
fit <- frm(bf(y | vint(size) ~ x) + fam, data = dd)
fixef(fit)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
