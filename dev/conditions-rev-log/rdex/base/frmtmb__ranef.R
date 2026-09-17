.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: ranef.frmtmb_fit
### Title: Extract random-effect modes
### Aliases: ranef.frmtmb_fit ranef

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# one matrix per random-effect term, levels by coefficients
ranef(fit)
ranef(fit)$g[1:3, ]                 # keyed by the grouping factor
attr(ranef(fit)$g, "term")          # the block it came from
ranef(fit)[["1 | g"]][1:3, ]        # or by that label

# condVar adds the conditional SDs a caterpillar plot needs
re <- as.data.frame(ranef(fit, condVar = TRUE))
head(re)
with(re[order(re$condval), ],
     plot(condval, seq_along(condval), pch = 16,
          xlim = range(condval - 2 * condsd, condval + 2 * condsd),
          xlab = "conditional mode", ylab = "group"))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
