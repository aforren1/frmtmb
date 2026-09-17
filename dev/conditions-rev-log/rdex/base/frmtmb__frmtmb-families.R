.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb-families
### Title: Additional response families
### Aliases: frmtmb-families student lognormal negbinomial nbinom1 Beta
###   tweedie compois zero_inflated_poisson zero_inflated_negbinomial
###   hurdle_poisson multinomial cumulative beta_binomial skew_normal
###   exgaussian bernoulli geometric exponential weibull shifted_lognormal
###   hurdle_gamma hurdle_lognormal zero_inflated_binomial
###   zero_inflated_beta asym_laplace zero_inflated_asym_laplace huber
###   sratio cratio acat von_mises categorical cox

### ** Examples

set.seed(4)
n <- 120
dd <- data.frame(x = rnorm(n))

# heavier tails than gaussian(), with an estimated df
dd$y <- 1 + 0.8 * dd$x + rt(n, df = 4)
fixef(frm(bf(y ~ x) + student(), data = dd))

# counts with more spread than poisson() allows
dd$cnt <- rnbinom(n, mu = exp(0.5 + 0.4 * dd$x), size = 2)
fit <- frm(bf(cnt ~ x) + negbinomial(), data = dd)
fixef(fit)$mu

# a zero-inflated count: the zi dpar gets its own predictor
dd$zi <- ifelse(runif(n) < 0.3, 0, dd$cnt)
frm(bf(zi ~ x, zi ~ 1) + zero_inflated_poisson(), data = dd)

# an ordered response: level order is the category order
dd$grade <- cut(1 + 0.8 * dd$x + rlogis(n), 3,
                labels = c("low", "mid", "high"), ordered_result = TRUE)
frm(bf(grade ~ x) + cumulative(), data = dd)

# a proportion in (0, 1)
dd$p <- plogis(0.2 + 0.6 * dd$x + rnorm(n, 0, 0.3))
frm(bf(p ~ x) + Beta(), data = dd)

# bounded influence: a few wild points barely move the slope
dd$rob <- 1 + 0.8 * dd$x + rnorm(n)
dd$rob[1:5] <- dd$rob[1:5] + 30
fixef(frm(bf(rob ~ x), family = huber(), data = dd))$mu
fixef(frm(bf(rob ~ x), family = gaussian(), data = dd))$mu

# an unordered factor: one predictor per non-reference category,
# named after the level it belongs to
dd$pick <- factor(sample(c("ale", "stout", "lager"), n, TRUE))
cat_fit <- frm(bf(pick ~ x), family = categorical(), data = dd)
fixef(cat_fit)                     # mulager and mustout; ale is the
                                   # reference
head(fitted(cat_fit))              # n x K category probabilities

# one category may take its own predictor
dd$w <- rnorm(n)
frm(bf(pick ~ x, mustout ~ w), family = categorical(), data = dd)

# an angle: mu is the mean direction, kappa the concentration
dd$angle <- atan2(sin(0.5 + dd$x), cos(0.5 + dd$x))
vm_fit <- frm(bf(angle ~ x), family = von_mises(), data = dd)
head(fitted(vm_fit))               # the mean direction, in radians

# proportional hazards with a spline baseline; (1 | g) is a frailty
dd$time <- rexp(n, exp(-0.5 + 0.7 * dd$x))
dd$out <- rbinom(n, 1, 0.3)        # 1 = right censored
cox_fit <- frm(bf(time | cens(out) ~ x), family = cox(), data = dd)
fixef(cox_fit)$mu                  # log hazard ratios
cox_baseline(cox_fit)              # the baseline hazard weights



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
