.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb-multimembership
### Title: Multi-membership random effects
### Aliases: frmtmb-multimembership

### ** Examples

set.seed(1)
n <- 200
d <- data.frame(
  x = rnorm(n),
  school1 = factor(sample(letters[1:8], n, TRUE)),
  school2 = factor(sample(letters[5:12], n, TRUE)),
  share1 = runif(n, 0.5, 1)
)
d$share2 <- 1 - d$share1
u <- rnorm(12, 0, 0.8)
names(u) <- letters[1:12]
d$y <- 1 + 0.5 * d$x +
  0.5 * u[as.character(d$school1)] +
  0.5 * u[as.character(d$school2)] + rnorm(n, 0, 0.5)

# equal membership: each pupil is half of each school
fit <- frm(bf(y ~ x + (1 | mm(school1, school2))) + gaussian(),
           data = d)
summary(fit)
# one coefficient per pooled school level
ranef(fit)

# the time each pupil spent in each school, as proportions
frm(bf(y ~ x + (1 | mm(school1, school2,
                       weights = cbind(share1, share2)))) + gaussian(),
    data = d)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
