.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: set_prior
### Title: Set up priors brms-style
### Aliases: set_prior

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(100), z = rnorm(100),
                 g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)

# `+` combines specifications; the class-wide one goes first so the
# coefficient-specific one can override it
pr <- set_prior("normal(0, 1)", class = "b") +
  set_prior("normal(0, 0.2)", class = "b", coef = "z") +
  set_prior("exponential(1)", class = "sd", group = "g")
pr

# the priors penalize the likelihood: the fit is a MAP estimate
fit <- frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd, prior = pr)
fixef(fit)$mu
# the tight prior on z shrinks it toward zero
fixef(frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd))$mu

# an empty distribution string sets a hard bound only
set_prior("", class = "b", coef = "x", lb = 0)

# a distributional parameter's own class is a density on the
# parameter itself, on its own scale, where the model gives that
# parameter no predictor. This model does not, so this is the
# spelling it offers, and get_prior() lists it
fit_s <- frm(bf(y ~ x + z) + gaussian(), data = dd,
             prior = set_prior("student_t(3, 0, 2.5)",
                               class = "sigma"))
sigma(fit_s)
# with sigma ~ 1 the model has a log-scale intercept instead, and
# that is the slot the prior addresses. brms draws the same line
set_prior("student_t(3, 0, 2.5)", class = "Intercept",
          dpar = "sigma")

# bounds address a nonlinear parameter the way a distribution does,
# so a guessing rate is held in [0, 1]
set_prior("", nlpar = "guess", lb = 0, ub = 1)

# the residual-correlation classes are brms's own names
set_prior("normal(0, 0.5)", class = "ar")
set_prior("lkj(2)", class = "rescor")
# and one internal covariance parameter, by its name
set_prior("", class = "theta", coef = "thetaac_1", lb = -2, ub = 2)

# class "cor" addresses a correlated block as a whole, brms's
# spelling; eta > 1 pulls the correlation toward zero
dd$z <- rnorm(100)
dd$y2 <- dd$y + rnorm(10, 0, 0.6)[dd$g] * dd$z
fitc <- frm(bf(y2 ~ x + z + (z | g)) + gaussian(), data = dd,
            prior = set_prior("lkj(4)", class = "cor"))
VarCorr(fitc)

# get_prior() shows which rows a design offers
get_prior(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd)

# prior() quotes its first argument, brms's spelling, and reaches
# the same machinery
prior(normal(0, 1), class = "b")



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
