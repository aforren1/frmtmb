.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: conditional_effects
### Title: Conditional effects of predictors
### Aliases: conditional_effects conditional_effects.frmtmb_fit

### ** Examples

set.seed(5)
dd <- data.frame(x = rnorm(120), f = factor(rep(c("a", "b"), 60)))
dd$y <- rnorm(120, 1 + 0.5 * dd$x + (dd$f == "b"), 1)
fit <- frm(bf(y ~ x * f), family = gaussian(), data = dd)
ce <- conditional_effects(fit, effects = c("x", "x:f"))
plot(ce, ask = FALSE)
# prediction intervals instead of epred bands
ce_p <- conditional_effects(fit, effects = "x", method = "predict")
## No test: 
# a likelihood-profile band: asymmetric, and no quadratic assumption
ce_pr <- conditional_effects(fit, effects = "x", band = "profile",
                             resolution = 20, profile_points = 5)

# a bootstrap band, reused for a second effect without refitting
ce_b <- conditional_effects(fit, effects = "x", band = "boot",
                            boot = 25, seed = 1)
ce_b2 <- conditional_effects(fit, effects = "x", band = "boot",
                             boot = attr(ce_b, "boot"))
## End(No test)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
