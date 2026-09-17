.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: mixture
### Title: Finite mixture families
### Aliases: mixture

### ** Examples

# two well-separated gaussian components
set.seed(3)
dd <- data.frame(y = c(rnorm(80, 0, 1), rnorm(80, 5, 1)),
                 x = rnorm(160))
fit <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian()), data = dd)
# one mu and sigma per component, plus the mixing weight theta1
fixef(fit)
# posterior class probability per observation
head(mixture_probs(fit))

# the mixing weight can take its own predictor
frm(bf(y ~ 1, theta1 ~ x) + mixture(gaussian(), gaussian()), data = dd)

## No test: 
# latent classes: every observation of a group shares one class
set.seed(4)
n_g <- 40
cls <- rep(c(1, 2), each = n_g / 2)
dg <- data.frame(g = factor(rep(seq_len(n_g), each = 5)))
dg$y <- rnorm(nrow(dg), c(0, 4)[cls[as.integer(dg$g)]], 1)
fg <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian(), groups = ~g),
          data = dg)
head(mixture_probs(fg))   # one row per group, not per observation
## End(No test)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
