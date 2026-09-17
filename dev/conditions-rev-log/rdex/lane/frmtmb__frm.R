.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm
### Title: Fit a model
### Aliases: frm

### ** Examples

## Not run: 
##D data(sleepstudy, package = "lme4")
##D fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
##D               data = sleepstudy)
##D summary(fit)
## End(Not run)
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
summary(fit)
fixef(fit)
VarCorr(fit)

# distributional regression: model sigma too
fit2 <- frm(bf(y ~ x + (1 | g), sigma ~ x) + gaussian(), data = dd)
anova(fit, fit2)

# a nonlinear model: discover the names, edit, fit. A nonlinear body
# is undefined at the zero start, so `start` is not optional here.
nd <- data.frame(t = rep(seq(0, 10, length.out = 25), 4))
nd$y <- 8 * (1 - exp(-nd$t / 3)) + rnorm(nrow(nd), 0, 0.3)
nf <- bf(y ~ asym * (1 - exp(-t / lrc)), asym ~ 1, lrc ~ 1,
         nl = TRUE) + gaussian()
st <- par_template(nf, data = nd)
st
st$beta[["asym_(Intercept)"]] <- 5
st$beta[["lrc_(Intercept)"]] <- 1
fixef(frm(nf, data = nd, start = st))

# or let a located prior place the same starts, brms-style
fixef(frm(nf, data = nd,
          prior = prior(normal(5, 5), nlpar = "asym") +
            prior(normal(1, 5), nlpar = "lrc")))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
